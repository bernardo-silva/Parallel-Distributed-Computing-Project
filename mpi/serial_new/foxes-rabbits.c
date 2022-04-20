#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
#include "environment.h"


#define EMPTY  0
#define ROCK   1
#define RABBIT 2
#define FOX    3

#define RED 0
#define BLACK 1


int choose_next_pos(struct Environment* env, int (*criteria)(), int i, int j, int* k, int* l){
    int left = j>0?-1:0;
    int right = j<env->N-1?1:0;

    int up = i>0?-1:0;
    int down = i<env->M-1?1:0;

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

    int choice = (i * env->N + j)%possible;

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
                // new.starve++;
                if(!(new.moved = choose_next_pos(env, position_empty, i, j, &k, &l))){
                    // env->temp_board[i][j].starve++;
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
    // env->temp_board[k][l] = new;
}

void merge_rows(Environment* env, int row_index, Entity* row){
    for(int i=0; i<env->N; i++){
        check_conflict(env, &row[i], &env->temp_board[row_index][i]);
    }
}

void reset_generation(Environment* env){
    for(int i=0; i<env->M; i++){
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

void update_generation(Environment* env, int color){
    for(int i=0; i<env->M; i++){
        for(int j=(i+color)%2; j<env->N; j+=2){
            move_entity(env, i, j);
        }
    }

}

void run_simulation(struct Environment* env){
    // Board -> Static
    // Temp - board -> Edit

    // Save current board and increase ages
    for(int i=0; i<env->M; i++){
        for(int j=0; j<env->N; j++){
            env->board[i][j].age++;
            env->temp_board[i][j].age++;
            env->board[i][j].starve++;
            env->temp_board[i][j].starve++;
        }
    }

    for(int gen = 0; gen<env->generations; ++gen){
        // Update red
        update_generation(env, RED);

        for(int i=0; i<env->M; i++){
            for(int j=0; j<env->N; j++)
                env->board[i][j] = env->temp_board[i][j];

        }
        print_board(env);

        // Update black
        update_generation(env, BLACK);

        reset_generation(env);
    }
}


int main (int argc, char *argv[])
{
    if (argc != 11) {
        printf("Incorrect number of parameters\n");
        return -1;
    }

    Environment env;
    generate_world(&env, argv);

    // printf("%d generated world!\n",id); 
    // print_board(&env);fflush(stdout);


    double exec_time = -omp_get_wtime();

    run_simulation(&env);

    exec_time  += omp_get_wtime();

    // print_board(&env);fflush(stdout);
    print_results(&env);
    fprintf(stderr, "%.2fs\n", exec_time);

    //Free memory
    for (int i = 0; i < env.M; i++) {
        free(env.board[i]);
        free(env.temp_board[i]);
    }
    free(env.board);
    free(env.temp_board);

    return 0;
}
