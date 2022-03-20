#ifndef ENVIRONMENT
#define ENVIRONMENT

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
// #include <worldGen.h>

#define EMPTY  0//' '
#define ROCK   1//'*'
#define RABBIT 2//'R'
#define FOX    3//'F'


typedef struct Entity{
    unsigned int type;
    unsigned int age;
    unsigned int starve;
    unsigned int moved;
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

int position_empty(Environment* const env, int i, int j);

int position_rabbit(Environment* env, int i, int j);

void insert_animal(Environment* env, int i, int j, char atype);

float r4_uni(uint32_t *seed);

void generate_element(Environment* env, int n, char atype, uint32_t *seed) ;

void generate_world(Environment* env, char *argv[]);

int kill_fox(Environment* env, int i, int j);
// void print_board(){
//     printf("-----------------\n");
//     for(int i=0; i<env.M; i++){
//         printf("|");
//         for(int j=0; j<env.N; j++){
//             printf(" %c |", env.board[i][j].type);
//         }
//         printf("\n");
//     }
//     printf("-----------------\n");
// }
// void print_temp_board(){
//     printf("-----------------\n");
//     for(int i=0; i<env.M; i++){
//         printf("|");
//         for(int j=0; j<env.N; j++){
//             printf(" %c %d %d %d|",
//                    env.temp_board[i][j].type,
//                    env.temp_board[i][j].age,
//                    env.temp_board[i][j].starve,
//                    env.temp_board[i][j].moved);
//         }
//         printf("\n");
//     }
//     printf("-----------------\n");
// }
void print_results( Environment* env);

#endif // !ENVIRONMENT