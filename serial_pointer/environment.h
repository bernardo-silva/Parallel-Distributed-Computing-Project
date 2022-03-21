#ifndef ENVIRONMENT
#define ENVIRONMENT

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

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

    Entity** board;
    Entity** temp_board;
} Environment;

int position_empty(Environment* env, int i, int j);

int position_rabbit(Environment* env, int i, int j);

void insert_animal(Environment* env, int i, int j, char atype);

float r4_uni(uint32_t *seed);

void generate_element(Environment* env, int n, char atype, uint32_t *seed) ;

void generate_world(Environment* env, char *argv[]);

int kill_fox(Environment* env, int i, int j);

void print_results( Environment* env);

#endif // !ENVIRONMENT
