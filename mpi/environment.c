#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
#include <mpi.h>
#include "environment.h"

#define EMPTY  0//' '
#define ROCK   1//'*'
#define RABBIT 2//'R'
#define FOX    3//'F'

#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n) - 1)
#define BLOCK_SIZE(id,p,n) \
  (BLOCK_HIGH(id,p,n) - BLOCK_LOW(id,p,n) + 1)
#define BLOCK_OWNER(index,p,n) \
(((p)*((index)+1)-1)/(n)) 
    
#define N_GHOST_LINES(id,p) (!id || id==p-1)? 1:2
#define BLOCK_LOW_GHOST(id, p, n) (BLOCK_LOW(id,p,n) - ((!id)?0:1))
#define BLOCK_HIGH_GHOST(id, p, n) (BLOCK_HIGH(id,p,n) + ((id != p-1)?1:0))


int position_empty(Environment* env, int i, int j) {
    return env->board[i][j].type == EMPTY;
}
int position_rabbit(Environment* env, int i, int j) {
    return env->board[i][j].type == RABBIT;
}
void insert_animal(Environment* env, int i, int j, char atype) {
    env->board[i][j].type = atype;
    env->temp_board[i][j].type = atype;
}

float r4_uni(uint32_t *seed)
{
    int seed_input, sseed;
    float r;

    seed_input = *seed;
    *seed ^= (*seed << 13);
    *seed ^= (*seed >> 17);
    *seed ^= (*seed << 5);
    sseed = *seed;
    r = 0.5 + 0.2328306e-09 * (seed_input + sseed);

    return r;
}

void generate_element(Environment* env, int n, char atype, uint32_t *seed)
{
    int i, j, k;

    for(k = 0; k < n; k++){
	i = env->M * r4_uni(seed);
	j = env->N * r4_uni(seed);

        if (i >= env->row_low_ghost && i <= env->row_high_ghost)
	if(position_empty(env, i-env->row_low_ghost, j))
	    insert_animal(env, i-env->row_low_ghost, j, atype);
    }
}

void generate_world(Environment* env, char *argv[], int id, int p){
    env->generations      = atoi(argv[1]);
    env->M                = atoi(argv[2]);
    env->N                = atoi(argv[3]);
    env->n_rocks          = atoi(argv[4]);
    env->n_rabbits        = atoi(argv[5]);
    env->rabbits_breeding = atoi(argv[6]);
    env->n_foxes          = atoi(argv[7]);
    env->foxes_breeding   = atoi(argv[8]);
    env->foxes_starvation = atoi(argv[9]);
    env->seed             = atoi(argv[10]);

    env->id = id;
    env->p = p;

    env->row_low = BLOCK_LOW(id, p, env->M);
    env->row_high = BLOCK_HIGH(id, p, env->M);
    env->block_size = BLOCK_SIZE(id, p, env->M);

    env->n_ghost_lines = N_GHOST_LINES(id, p);
    env->row_low_ghost = BLOCK_LOW_GHOST(id, p, env->M);
    env->row_high_ghost = BLOCK_HIGH_GHOST(id, p, env->M);
    env->block_size_ghost = env->block_size + env->n_ghost_lines;
    
    env->is_not_top = (!id)?0:1;
    env->is_not_bottom = (id == p-1)?0:1;

    // printf("%d with row low %d and high %d\n",
    //        id,env->row_low,env->row_high);
    //
    // printf("%d with block ghost size %d and n_ghost_lines %d\n",
    //        id,env->block_size_ghost,env->n_ghost_lines);
    //
    // printf("%d with row ghots low %d and high %d\n",
    //        id,env->row_low_ghost,env->row_high_ghost);

    env->board = malloc(sizeof *env->board * env->block_size_ghost);
    env->temp_board = malloc(sizeof *env->board * env->block_size_ghost);

    for (int i=0;i<env->block_size_ghost;i++){
        env->board[i] = malloc(env->N * sizeof(Entity));
        env->temp_board[i] = malloc(env->N * sizeof(Entity));
    }

    for (int i=0;i<env->block_size_ghost;i++){
        for (int j=0;j<env->N;j++){
            env->board[i][j].type = EMPTY;
            env->board[i][j].age = 0;
            env->board[i][j].starve = 0;
            env->board[i][j].moved = 0;
            env->temp_board[i][j].type = EMPTY;
            env->temp_board[i][j].age = 0;
            env->temp_board[i][j].starve = 0;
            env->temp_board[i][j].moved = 0;
        }
    }

    generate_element(env, env->n_rocks, ROCK, &env->seed);
    generate_element(env, env->n_rabbits, RABBIT, &env->seed);
    generate_element(env, env->n_foxes, FOX, &env->seed);

}

int kill_fox(Environment* env, int i, int j){
    if(env->temp_board[i][j].starve >= env->foxes_starvation) {
        env->temp_board[i][j] = (Entity) {.type = EMPTY, .age=0, .starve=0, .moved=0};
        env->board[i][j] = (Entity) {.type = EMPTY, .age=0, .starve=0, .moved=0};
        return 1;
    }
    return 0;
}
void print_board(Environment* env){
    char * animal = " *RF";
    printf("---------------------\n");fflush(stdout);
    printf("   ");
    for(int j=0; j<env->N; j++)
        printf("%02d|",j);

    printf("\n");

    for(int i=0; i<env->block_size_ghost; i++){
        printf("%02d:",i);
        for(int j=0; j<env->N; j++){
            printf(" %c|", animal[env->temp_board[i][j].type]);
        }
        printf("\n");fflush(stdout);
    }
    printf("---------------------\n");fflush(stdout);
}

void print_temp_board(Environment* env){
    printf("-----------------\n");
    for(int i=0; i<env->M; i++){
        printf("|");
        for(int j=0; j<env->N; j++){
            printf(" %c %d %d %d|",
                   env->temp_board[i][j].type,
                   env->temp_board[i][j].age,
                   env->temp_board[i][j].starve,
                   env->temp_board[i][j].moved);
        }
        printf("\n");
    }
    printf("-----------------\n");
}

void print_results(Environment* env){
    int rocks=0, rabbits=0, foxes=0;
    for(int i=env->is_not_top; i<env->block_size + env->is_not_top; i++){
        for(int j=0; j<env->N; j++){
            if(env->board[i][j].type == ROCK){ rocks++; continue; }
            if(env->board[i][j].type == RABBIT){ rabbits++; continue; }
            if(env->board[i][j].type == FOX){ foxes++; continue; }
        }
    }
    printf("%d %d %d\n", rocks, rabbits, foxes);
}
