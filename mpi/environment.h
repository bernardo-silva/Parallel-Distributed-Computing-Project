#ifndef ENVIRONMENT
#define ENVIRONMENT

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <mpi.h>

#define EMPTY  0//' '
#define ROCK   1//'*'
#define RABBIT 2//'R'
#define FOX    3//'F'

typedef struct Entity{
    int type;
    int age;
    int starve;
    int moved;
} Entity;

typedef struct Environment {
    int generations;
    int M;
    int N;
    int n_rocks;
    int n_rabbits;
    int rabbits_breeding;
    int n_foxes;
    int foxes_breeding;
    int foxes_starvation;
    uint32_t seed;

    int id;
    int p;

    int row_low;
    int row_high;
    int column_low;
    int column_high;

    int row_block_size;
    int column_block_size;

    int n_ghost_rows;
    int n_ghost_columns;

    int row_low_ghost;
    int row_high_ghost;
    int column_low_ghost;
    int column_high_ghost;

    int row_block_size_ghost;
    int column_block_size_ghost;

    int is_not_bottom;
    int is_not_top;
    int is_not_right;
    int is_not_left;

    Entity** board;
    Entity** temp_board;
} Environment;

int position_empty(Environment* env, int i, int j);

int position_rabbit(Environment* env, int i, int j);

void insert_animal(Environment* env, int i, int j, char atype);

float r4_uni(uint32_t *seed);

void generate_element(Environment* env, int n, char atype, uint32_t *seed) ;

void generate_world(Environment* env, char *argv[], int id, int p);

int kill_fox(Environment* env, int i, int j);

void print_board( Environment* env);

void print_temp_board( Environment* env);

void get_results(Environment* env, int* rocks, int* rabbits, int* foxes);

#endif // !ENVIRONMENT
