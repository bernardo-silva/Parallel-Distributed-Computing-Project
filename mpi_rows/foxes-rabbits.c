#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
#include "environment.h"

#include <mpi.h>

#define EMPTY  0
#define ROCK   1
#define RABBIT 2
#define FOX    3

#define RED 0
#define BLACK 1

MPI_Datatype MPI_Entity;

#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n) - 1)
#define BLOCK_SIZE(id,p,n) \
  (BLOCK_HIGH(id,p,n) - BLOCK_LOW(id,p,n) + 1)
#define BLOCK_OWNER(index,p,n) \
(((p)*((index)+1)-1)/(n)) 


int choose_next_pos(struct Environment* env, int (*criteria)(), int i, int j, int* k, int* l){
    int left = j>0?-1:0;
    int right = j<env->N-1?1:0;

    int up = i>0?-1:0;
    int down = i<env->block_size_ghost-1?1:0;

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

    int choice = ((i+env->row_low-env->is_not_top) * env->N + j)%possible;

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
    for(int i=0; i<env->N; i++){
        if(row[i].type != EMPTY)
            check_conflict(env, &env->temp_board[row_index][i], &row[i]);
    }
}

void reset_generation(Environment* env){
    for(int i=0; i<env->block_size_ghost; i++){
        for(int j=0; j<env->N; j++){
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

void recieve_above_below(Environment* env, Entity* above, int tag_above, Entity*
                      below, int tag_below, MPI_Request* requests, int
                      request_index_above, int request_index_below){

    // MPI request ghost line from process above
    if(env->is_not_top)
        MPI_Irecv(above, env->N, MPI_Entity, env->id - 1, tag_above,
                  MPI_COMM_WORLD, requests + request_index_above);

    // MPI request ghost line from process below
    if(env->is_not_bottom)
        MPI_Irecv(below, env->N, MPI_Entity, env->id + 1, tag_below,
                  MPI_COMM_WORLD, requests + request_index_below);
}

void send_above_below(Environment* env, Entity* above, int tag_above, Entity*
                      below, int tag_below, MPI_Request* requests, int
                      request_index_above, int request_index_below){
    if(env->is_not_top){
        MPI_Isend(above, env->N, MPI_Entity, env->id - 1, tag_above,
                  MPI_COMM_WORLD, requests + request_index_above);
    }
    if(env->is_not_bottom){
        MPI_Isend(below, env->N, MPI_Entity, env->id + 1, tag_below,
                  MPI_COMM_WORLD, requests + request_index_below);
    }
}

void empty_ghost_lines(Environment* env){
    if(env->is_not_top)
        for(int j=0; j<env->N;j++)
            env->temp_board[0][j] = (Entity) {.type = EMPTY, .age=0, .starve=0, .moved=0};

    if(env->is_not_bottom)
        for(int j=0; j<env->N;j++)
            env->temp_board[env->block_size_ghost - 1][j] = (Entity) {.type =
                EMPTY, .age=0, .starve=0, .moved=0};
}

void update_generation(Environment* env, int color, Entity* row_below, Entity* row_above){
    // extern MPI_Datatype MPI_Entity;
    MPI_Request requests1[4] = {[0 ... 3] = MPI_REQUEST_NULL};
    MPI_Request requests2[4] = {[0 ... 3] = MPI_REQUEST_NULL};
    int index;

    empty_ghost_lines(env);
    
    // MPI recieve ghost line from process above and below
    recieve_above_below(env, row_above, 0, row_below, 1, requests1, 0, 1);

    // Update block (without updating ghost)
    for(int i=env->is_not_top; i<env->block_size + env->is_not_top; i++){
        for(int j=(i+env->row_low+color-env->is_not_top)%2; j<env->N; j+= 2){
            move_entity(env, i, j);
        }
    }

    send_above_below(env, env->temp_board[0], 1,
                     env->temp_board[env->block_size_ghost - 1], 0, requests1, 2, 3);

    recieve_above_below(env, env->temp_board[0], 3,
                        env->temp_board[env->block_size_ghost-1], 2,
                        requests2, 2, 3);

    for(int k=0; k<env->is_not_bottom+env->is_not_top;k++){
        //Wait for any of the recieves
        MPI_Waitany(2, requests1, &index, MPI_STATUS_IGNORE);
        //If index is 0, it recieved a row from above
        if (!index){
            merge_rows(env, 1, row_above);
            MPI_Isend(env->temp_board[1], env->N, MPI_Entity, env->id - 1, 2,
                      MPI_COMM_WORLD, requests2);
        }
        else{ //Else, it recieved a row from below
            merge_rows(env, env->block_size_ghost - 2, row_below);
            MPI_Isend(env->temp_board[env->block_size_ghost - 2], env->N,
                      MPI_Entity, env->id + 1, 3, MPI_COMM_WORLD, requests2 + 1);
        }
    }

    MPI_Waitall(2, requests2+2, MPI_STATUS_IGNORE);
}

void run_simulation(struct Environment* env){
    // Board -> Static
    // Temp - board -> Edit

    Entity* row_below = malloc(env->N * sizeof(Entity));
    Entity* row_above = malloc(env->N * sizeof(Entity));

    // Save current board and increase ages
    for(int i=0; i<env->block_size_ghost; i++){
        for(int j=0; j<env->N; j++){
            env->board[i][j].age++;
            env->temp_board[i][j].age++;
            env->board[i][j].starve++;
            env->temp_board[i][j].starve++;
        }
    }

    for(int gen = 0; gen<env->generations; ++gen){
        // Update red
        update_generation(env, RED, row_below, row_above);

        for(int i=0; i<env->block_size_ghost; i++){
            for(int j=0; j<env->N; j++)
                env->board[i][j] = env->temp_board[i][j];
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // Update black
        update_generation(env, BLACK, row_below, row_above);

        MPI_Barrier(MPI_COMM_WORLD);

        reset_generation(env);
    }
    free(row_below);
    free(row_above);
}


int main (int argc, char *argv[])
{
    int id, p;
    // MPI_Datatype MPI_Entity;

    MPI_Init (&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (argc != 11) {
        printf("Incorrect number of parameters\n");
        return -1;
    }

    const int nfields = 1;
    const int block_lens[] = {4};
    MPI_Aint displacements = offsetof(Entity, type);
    MPI_Datatype types[] = {MPI_INT}; 
    MPI_Type_create_struct(nfields, block_lens, &displacements, types, &MPI_Entity);
    MPI_Type_commit(&MPI_Entity);

    Environment env;
    generate_world(&env, argv, id, p);
    // printf("%d generated world!\n",id);
    // print_board(&env);fflush(stdout);


    MPI_Barrier(MPI_COMM_WORLD);
    double exec_time = -MPI_Wtime();

    run_simulation(&env);

    MPI_Barrier(MPI_COMM_WORLD);
    exec_time  += MPI_Wtime();

    // printf("%d final world!\n",id);
    // print_board(&env);fflush(stdout);
    int rocks, rabbits, foxes;
    int t_rocks=0, t_rabbits=0, t_foxes=0;
    
    get_results(&env, &rocks, &rabbits, &foxes);

    MPI_Reduce(&rocks, &t_rocks, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&rabbits, &t_rabbits, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&foxes, &t_foxes, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(!env.id){
        printf("%d %d %d\n", t_rocks, t_rabbits, t_foxes);
        fprintf(stderr, "%.8fs\n", exec_time);
    }

    //Free memory
    for (int i = 0; i < env.block_size_ghost; i++) {
        free(env.board[i]);
        free(env.temp_board[i]);
    }
    free(env.board);
    free(env.temp_board);

    MPI_Type_free(&MPI_Entity);
    MPI_Finalize();
    return 0;
}
