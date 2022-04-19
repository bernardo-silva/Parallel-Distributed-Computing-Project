#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
// #include <string.h>
#include <omp.h>
#include "environment.h"

#define EMPTY  0//' '
#define ROCK   1//'*'
#define RABBIT 2//'R'
#define FOX    3//'F'

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

void check_conflict(struct Environment* env, int k, int l, Entity* ent){
    Entity *previous = &env->temp_board[k][l];

    //Do nothing if previous was rock or empty (rock should never happen)
    if(previous->type == EMPTY || previous->type == ROCK) return;

    // Check if they are the same type
    if(previous->type == ent->type){
        //Both Foxes
        if(ent->type == FOX){
            if(!ent->starve || !previous->starve) ent->starve = 0;
            if(ent->age == previous->age){
                ent->starve = (ent->starve < previous->starve)?ent->starve:previous->starve;
                return;
            }
        }
        ent->age = (ent->age > previous->age)?ent->age:previous->age;
        return;
    }
    //If not, pick the fox
    if(previous->type == FOX){
        *ent = *previous;
    }
    ent->starve = 0;

    return;
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
    check_conflict(env, k, l, &new);
    env->temp_board[k][l] = new;
}

void reset_generation(struct Environment* env){
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
    // printf("Generation %d, red\n",0);
    // print_board(env);

    for(int gen = 0; gen<env->generations; ++gen){
        // Update red
        for(int i=0; i<env->M; i++){
            for(int j=i%2; j<env->N; j+=2){
                move_entity(env, i, j);
            }
        }
        // printf("Red updated\n");

        for(int i=0; i<env->M; i++){
            for(int j=0; j<env->N; j++)
                env->board[i][j] = env->temp_board[i][j];
        }

        // printf("Generation %d, red\n",gen+1);
        // print_board(env);
        // Update black
        for(int i=0; i<env->M; i++){
            for(int j=(i+1)%2; j<env->N; j+=2){
                move_entity(env, i, j);
            }
        }
        // printf("Updated black\n");

        reset_generation(env);
        // printf("Generation %d, black\n",gen+1);
        // print_board(env);
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

    // printf("Generated world\n");
    // print_board(&env);

    double exec_time = -omp_get_wtime();

    run_simulation(&env);

    exec_time  += omp_get_wtime();

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
