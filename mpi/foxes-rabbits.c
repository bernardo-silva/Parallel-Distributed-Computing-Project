#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
#include <mpi.h>

#include "environment.h"
#include "mpi_utils.h"

#define EMPTY  0
#define ROCK   1
#define RABBIT 2
#define FOX    3

#define RED 0
#define BLACK 1

#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

MPI_Datatype MPI_Entity;

int choose_next_pos(struct Environment* env, int (*criteria)(), int i, int j, int* k, int* l){
    int up = i>0?-1:0;
    int down = i<env->row_block_size_ghost-1?1:0;

    int left = j>0?-1:0;
    int right = j<env->column_block_size_ghost-1?1:0;

    int possible = 0;
    int moves[4][2] = {{0,0},{0,0},{0,0},{0,0}};

    if(up) if(criteria(env, i+up, j)){
        moves[possible][0] = i+up;
        moves[possible][1] = j;
        possible++;
    }
    if(right) if(criteria(env, i, j+right)){
        moves[possible][0] = i;
        moves[possible][1] = j+right;
        possible++;
    }
    if(down) if(criteria(env, i+down, j)){
        moves[possible][0] = i+down;
        moves[possible][1] = j;
        possible++;
    }
    if(left) if(criteria(env, i, j+left)){
        moves[possible][0] = i;
        moves[possible][1] = j+left;
        possible++;
    }
    if(!possible) return 0;


    int choice;
    choice = ((i+env->row_low-env->is_not_top) * env->N +
        j+env->column_low-env->is_not_left) % possible;

    *k = moves[choice][0];
    *l = moves[choice][1];

    return 1;
}

void check_conflict(struct Environment* env, Entity* on_board, Entity* new){
    //Do nothing if on_board was rock or empty (rock should never happen)
    if(on_board->type == EMPTY){
        *on_board = *new;
        return;
    }

    // Check if they are the same type
    if(on_board->type == new->type){
        //Both Foxes
        if(new->type == FOX){
            if(!new->starve || !on_board->starve) on_board->starve = 0;
            if(new->age == on_board->age){
                on_board->starve = (new->starve < on_board->starve)?new->starve:on_board->starve;
                return;
            }
        }
        on_board->age = (new->age > on_board->age)?new->age:on_board->age;
        return;
    }
    //If not, pick the fox
    if(new->type == FOX){
        *on_board = *new;
    }
    on_board->starve = 0;
}

void move_entity(Environment* env, int i, int j){
    int k = 0 , l = 0;

    Entity *ent = &env->board[i][j];
    Entity new = env->board[i][j];

    int type = ent->type;
    int breed_age = 0;

    if (type == EMPTY || type == ROCK || ent->moved) return;

    switch (type) {
        case RABBIT:
            if(!(new.moved = choose_next_pos(env, position_empty, i, j, &k, &l)))
                return;
            breed_age = env->rabbits_breeding;
            break;

        case FOX:
            if(!(new.moved = choose_next_pos(env, position_rabbit, i, j, &k, &l))){
                if(!(new.moved = choose_next_pos(env, position_empty, i, j, &k, &l))){
                    return;
                }
            }
            else new.starve = 0;
            breed_age = env->foxes_breeding;
            break;
    }
    if(ent->age >= breed_age){
        new.age = 0;
        env->temp_board[i][j] = (Entity) {.type = type, .age=0, .starve=0, .moved=0};
    }
    else{
        //Otherwise, leave empty entity behind
        env->temp_board[i][j] = (Entity) {.type = EMPTY, .age=0, .starve=0, .moved=0};
    }

    check_conflict(env, &env->temp_board[k][l], &new);
}

void merge_rows(Environment* env, int row_index, Entity* row){
    for(int i=env->is_not_left; i<env->column_block_size+env->is_not_left; i++){
        if(row[i].type != EMPTY)
            check_conflict(env, &env->temp_board[row_index][i], &row[i]);
    }
}

void merge_columns(Environment* env, int col_index, Entity* col){
    for(int i=env->is_not_top; i<env->row_block_size+env->is_not_top; i++){
        if(col[i].type != EMPTY)
            check_conflict(env, &env->temp_board[i][col_index], &col[i]);
    }
}

void reset_generation(Environment* env){
    for(int i=0; i<env->row_block_size_ghost; i++){
        for(int j=0; j<env->column_block_size_ghost; j++){
            switch (env->temp_board[i][j].type) {
                case EMPTY:
                    env->board[i][j] = env->temp_board[i][j];
                    continue;
                case ROCK:
                    continue;
                case FOX:
                    if(kill_fox(env, i, j)) continue;
                    env->temp_board[i][j].starve++;
                    break;
                default:
                    break;
            }
            env->temp_board[i][j].moved = 0;
            env->temp_board[i][j].age++;
            env->board[i][j] = env->temp_board[i][j];
        }
    }
}

void update_generation(Environment* env, int color, 
                       Entity* row_below, Entity* row_above,
                       Entity* column_left, Entity* column_right,
                       Entity* column_left_send, Entity* column_right_send){

    MPI_Request requests1[8] = {[0 ... 7] = MPI_REQUEST_NULL};
    MPI_Request requests2[8] = {[0 ... 7] = MPI_REQUEST_NULL};
    int index;

    empty_ghost_lines(env);
    
    // MPI recieve ghost line from neighbor processes
    recieve_from_neigs(env, row_above, UP, row_below, DOWN, column_left,
                       LEFT, column_right, RIGHT, requests1, 0, 1, 2, 3);

    // Update block (without updating ghost)
    for(int i=env->is_not_top; i<env->row_block_size + env->is_not_top; i++){
        int start = (i+color+env->row_low-env->is_not_top)%2; 
        for(int j=start; j<env->column_block_size+env->is_not_left; j+= 2){
            move_entity(env, i, j);
        }
    }

    //Convert columns to rows to be able to send
    columns_to_rows(env, column_left_send, column_right_send);

    send_to_neigs(env, 
                  env->temp_board[0], DOWN,
                  env->temp_board[env->row_block_size_ghost - 1], UP, 
                  column_left_send, RIGHT,
                  column_right_send, LEFT,
                  requests1, 4, 5, 6, 7);

    recieve_from_neigs(env,
                       env->temp_board[0], 21,
                       env->temp_board[env->row_block_size_ghost-1], 20,
                       column_left_send, 23,
                       column_right_send, 22,
                       requests2, 4, 5, 6, 7);

    for(int k=0; k<env->n_neigs; k++){
        MPI_Waitany(4, requests1, &index, MPI_STATUS_IGNORE);
        switch (index) {
            case UP:{
                merge_rows(env, 1, row_above);
                MPI_Isend(env->temp_board[1], env->column_block_size_ghost,
                          MPI_Entity, env->neigs_ids[UP], 20, env->cart_comm, requests2);
                break;
            }
            case DOWN:{
                merge_rows(env, env->row_block_size_ghost - 2, row_below);
                MPI_Isend(env->temp_board[env->row_block_size_ghost - 2],
                          env->column_block_size_ghost, MPI_Entity,
                          env->neigs_ids[DOWN], 21, env->cart_comm, requests2 + 1);
                break;
            }
            case LEFT:{
                merge_columns(env, 1, column_left);
                column_to_row(env, 1, column_left);
                MPI_Isend(column_left, env->row_block_size_ghost,
                          MPI_Entity, env->neigs_ids[LEFT], 22, env->cart_comm,
                          requests2 + 2);
                break;
            }
            case RIGHT:{
                merge_columns(env, env->column_block_size_ghost-2, column_right);
                column_to_row(env, env->column_block_size_ghost-2, column_right);
                MPI_Isend(column_right, env->row_block_size_ghost,
                          MPI_Entity, env->neigs_ids[RIGHT], 23,
                          env->cart_comm, requests2 + 3);
                break;
            }
            case MPI_UNDEFINED:{
                printf("%d: Error undefined request\n",env->id);fflush(stdout);
            }
        }
    }
    // MPI_Waitall(4, requests1 + 4, MPI_STATUS_IGNORE);
    
    // printf("%d === %d %d %d %d\n", env->id, requests1[4] == MPI_REQUEST_NULL,
    //        requests1[5] == MPI_REQUEST_NULL, requests1[6] == MPI_REQUEST_NULL,
    //        requests1[7] == MPI_REQUEST_NULL);fflush(stdout);
    // MPI_Waitall(8, requests1, MPI_STATUS_IGNORE);

    MPI_Waitall(8, requests2, MPI_STATUS_IGNORE);

    rows_to_columns(env, column_left_send, column_right_send);
}

void run_simulation(struct Environment* env){
    // Board -> Static
    // Temp - board -> Edit
    Entity* row_below = malloc(env->column_block_size_ghost * sizeof(Entity));
    Entity* row_above = malloc(env->column_block_size_ghost * sizeof(Entity));

    Entity* column_left  = malloc(env->row_block_size_ghost * sizeof(Entity));
    Entity* column_right = malloc(env->row_block_size_ghost * sizeof(Entity));
    Entity* column_left_send  = malloc(env->row_block_size_ghost * sizeof(Entity));
    Entity* column_right_send = malloc(env->row_block_size_ghost * sizeof(Entity));

    // Save current board and increase ages
    for(int i=0; i<env->row_block_size_ghost; i++){
        for(int j=0; j<env->column_block_size_ghost; j++){
            env->board[i][j].age++;
            env->temp_board[i][j].age++;
            env->board[i][j].starve++;
            env->temp_board[i][j].starve++;
        }
    }

    // printf("%d starting to update\n",env->id);fflush(stdout);
    for(int gen = 0; gen<env->generations; ++gen){
        // Update red
        update_generation(env, RED, row_below, row_above,
                          column_left, column_right,
                          column_left_send, column_right_send);
        // printf("%d updated red\n",env->id);fflush(stdout);

        for(int i=0; i<env->row_block_size_ghost; i++){
            for(int j=0; j<env->column_block_size_ghost; j++)
                env->board[i][j] = env->temp_board[i][j];
        }

        MPI_Barrier(env->cart_comm);

        // Update black
        update_generation(env, BLACK, row_below, row_above,
                          column_left, column_right,
                          column_left_send, column_right_send);

        MPI_Barrier(env->cart_comm);//WHYYYY????

        reset_generation(env);
    }
    free(row_below);
    free(row_above);
    free(column_left);
    free(column_right);
    free(column_left_send);
    free(column_right_send);
}


int main (int argc, char *argv[])
{
    MPI_Comm cart_comm;
    int id, p;
    int grid_dims[2];
    grid_dims[0] = grid_dims[1] = 0;
    int periodic[2];
    periodic[0] = periodic[1] = 0;
    // MPI_Datatype MPI_Entity;

    MPI_Init (&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    MPI_Dims_create(p, 2, grid_dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, grid_dims, periodic, 1, &cart_comm);

    if (argc != 11) {
        printf("Incorrect number of parameters\n");
        return -1;
    }

    Environment env;
    generate_world(&env, argv, id, p, grid_dims, cart_comm);
    // printf("%d generated world!\n",id); print_board(&env);fflush(stdout); printf("%d is not %d %d %d %d\n",env.id, 
    //        env.is_not_top,
    //        env.is_not_bottom,
    //        env.is_not_left,
    //        env.is_not_right);
    // printf("%d has neigs %d %d %d %d\n",env.id, env.neigs_ids[0],
    //        env.neigs_ids[1], env.neigs_ids[2], env.neigs_ids[3]);
    printf("%d has row low %d and column low %d \n",env.id, env.row_low,
           env.column_low);
    fflush(stdout);


    MPI_Barrier(env.cart_comm);

    double exec_time = -MPI_Wtime();

    run_simulation(&env);

    MPI_Barrier(env.cart_comm);
    exec_time  += MPI_Wtime();

    // printf("%d final world!\n",id);
    // print_board(&env);fflush(stdout);
    int rocks, rabbits, foxes;
    int t_rocks=0, t_rabbits=0, t_foxes=0;
    
    get_results(&env, &rocks, &rabbits, &foxes);
    printf("%d: %d %d %d\n",env.id, rocks, rabbits, foxes);
    MPI_Barrier(env.cart_comm);

    MPI_Reduce(&rocks, &t_rocks, 1, MPI_INT, MPI_SUM, 0, env.cart_comm);
    MPI_Reduce(&rabbits, &t_rabbits, 1, MPI_INT, MPI_SUM, 0, env.cart_comm);
    MPI_Reduce(&foxes, &t_foxes, 1, MPI_INT, MPI_SUM, 0, env.cart_comm);

    if(!env.id){
        printf("%d %d %d\n", t_rocks, t_rabbits, t_foxes);
        fprintf(stderr, "%.8fs\n", exec_time);
    }

    printf("%d generated world!\n",id);
    print_board(&env);fflush(stdout);
    //Free memory
    for (int i = 0; i < env.row_block_size_ghost; i++) {
        free(env.board[i]);
        free(env.temp_board[i]);
    }
    free(env.board);
    free(env.temp_board);

    MPI_Type_free(&MPI_Entity);
    MPI_Finalize();
    return 0;
}
