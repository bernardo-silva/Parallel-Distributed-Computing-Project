#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
// #include <worldGen.h>
#include <omp.h>

#define EMPTY  0//' '
#define ROCK   1//'*'
#define RABBIT 2//'R'
#define FOX    3//'F'


struct Entity{
    int type;
    unsigned int age;
    unsigned int starve;
    unsigned int moved;
};

struct Environment {
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

    struct Entity** board;
    struct Entity** temp_board;
} env;

int position_empty(int i, int j) {
    return !env.board[i][j].type;
}
int position_rabbit(int i, int j) {
    return env.board[i][j].type == RABBIT;
}
void insert_animal(int i, int j, char atype) {
    env.board[i][j].type = atype;
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

void generate_element(int n, char atype, uint32_t *seed)
{
    int i, j, k;

    for(k = 0; k < n; k++){
	i = env.M * r4_uni(seed);
	j = env.N * r4_uni(seed);
	if(position_empty(i, j))
	    insert_animal(i, j, atype);
    }
}

void generate_world(char *argv[]){
    env.generations      = atoi(argv[1]);
    env.M                = atoi(argv[2]);
    env.N                = atoi(argv[3]);
    env.n_rocks          = atoi(argv[4]);
    env.n_rabbits        = atoi(argv[5]);
    env.rabbits_breeding = atoi(argv[6]);
    env.n_foxes          = atoi(argv[7]);
    env.foxes_breeding   = atoi(argv[8]);
    env.foxes_starvation = atoi(argv[9]);
    env.seed             = atoi(argv[10]);

    env.board = malloc(sizeof *env.board * env.M);
    env.temp_board = malloc(sizeof *env.temp_board * env.M);

    for (int i=0;i<env.M;i++){
        env.board[i] = malloc(env.N * sizeof(struct Entity));
        env.temp_board[i] = malloc(env.N * sizeof(struct Entity));
    }
    for (int i=0;i<env.M;i++){
        for (int j=0;j<env.N;j++){
            env.board[i][j].type = EMPTY;
            env.board[i][j].age = 0;
            env.board[i][j].starve = 0;
            env.board[i][j].moved = 0;
            env.temp_board[i][j].type = EMPTY;
            env.temp_board[i][j].age = 0;
            env.temp_board[i][j].starve = 0;
            env.temp_board[i][j].moved = 0;
        }
    }

    generate_element(env.n_rocks, ROCK, &env.seed);
    generate_element(env.n_rabbits, RABBIT, &env.seed);
    generate_element(env.n_foxes, FOX, &env.seed);
}

void print_board(){
    char * labels = " *RF";
    printf("-----------------\n");
    for(int i=0; i<env.M; i++){
        printf("|");
        for(int j=0; j<env.N; j++){
            printf(" %c |", labels[env.board[i][j].type]);
        }
        printf("\n");
    }
    printf("-----------------\n");
}
void print_temp_board(){
    char * labels = " *RF";
    printf("-----------------\n");
    for(int i=0; i<env.M; i++){
        printf("|");
        for(int j=0; j<env.N; j++){
            printf(" %c %d %d %d|",
                   labels[env.temp_board[i][j].type],
                   env.temp_board[i][j].age,
                   env.temp_board[i][j].starve,
                   env.temp_board[i][j].moved);
        }
        printf("\n");
    }
    printf("-----------------\n");
}
void print_results(){
    int rocks=0, rabbits=0, foxes=0;
    for(int i=0; i<env.M; i++){
        for(int j=0; j<env.N; j++){
            if(env.board[i][j].type == ROCK){ rocks++; continue; }
            if(env.board[i][j].type == RABBIT){ rabbits++; continue; }
            if(env.board[i][j].type == FOX){ foxes++; continue; }
        }
    }
    printf("%d %d %d\n", rocks, rabbits, foxes);
}

int choose_next_pos(int (*criteria)(), int i, int j, int* k, int* l){
    int left = j>0?-1:0;
    int right = j<env.N-1?1:0;

    int up = i>0?-1:0;
    int down = i<env.M-1?1:0;

    int possible = 0;
    int moves[4][2] = {{0,0},{0,0},{0,0},{0,0}};


    if(up) if(criteria(i+up, j)){
        moves[possible][0] = i+up;
        moves[possible][1] = j;
        possible++;
    }
    if(right) if(criteria(i, j+right)){
        moves[possible][0] = i;
        moves[possible][1] = j+right;
        possible++;
    }
    if(down) if(criteria(i+down, j)){
        moves[possible][0] = i+down;
        moves[possible][1] = j;
        possible++;
    }
    if(left) if(criteria(i, j+left)){
        moves[possible][0] = i;
        moves[possible][1] = j+left;
        possible++;
    }
    if(!possible) return 0;//{printf("Not possible\n"); return 0;}

    int choice = (i * env.N + j)%possible;
    // choice %= possible;
    *k = moves[choice][0];
    *l = moves[choice][1];
    // printf("Moving %c on %d %d with possible %d and choice %d\n",
           // env.board[i][j].type, i,j,possible, choice);

    return 1;
}

void check_conflict(int k, int l, struct Entity* ent){
    struct Entity *previous = &env.temp_board[k][l];

    //Do nothing if previous was rock or empty (rock should never happen)
    if(previous->type == EMPTY || previous->type == ROCK) return;

    // Check if they are the same time and pick the one with the highest age
    if(previous->type == ent->type){
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

void move_entity(int i, int j){
    int k = 0 , l = 0;
    struct Entity *ent = &env.board[i][j];
    struct Entity new = env.board[i][j];
    int type = ent->type;

    // Do nothing if rock or empty or if moved previously
    if (type == ROCK || type == EMPTY || ent->moved) return;

    // Pick move criteria based on entity type
    int (*criteria)() = (type == RABBIT? position_empty : position_rabbit);

    //Check if moved and save new pos on (k,l)
    new.moved = choose_next_pos(criteria, i, j, &k, &l);

    //For foxes, try to move again to empty spot
    if(type==FOX && !new.moved){
        new.moved = choose_next_pos(position_empty, i, j, &k, &l);
        new.starve++;
    }
    //If first move was successful, reset starvation
    else new.starve = 0;

    //If no move, do nothing
    if (!new.moved) return;

    //If it reached the breading age, reset age and keep a copy there
    if(ent->age >= (type==RABBIT?env.rabbits_breeding:env.foxes_breeding)){
        new.age = 0;
        env.temp_board[i][j] = (struct Entity) {.type = type, .age=0, .starve=0, .moved=0};
    }
    else{
        //Otherwise, leave empty entity behind
        env.temp_board[i][j] = (struct Entity) {.type = EMPTY, .age=0, .starve=0, .moved=0};
    }
    check_conflict(k, l, &new);
    env.temp_board[k][l] = new;
}

int kill_fox(int i, int j){
    if(env.temp_board[i][j].starve >= env.foxes_starvation) {
        env.temp_board[i][j] = (struct Entity) {.type = EMPTY, .age=0, .starve=0, .moved=0};
        env.board[i][j] = (struct Entity) {.type = EMPTY, .age=0, .starve=0, .moved=0};
        return 1;
    }
    return 0;
}

void reset_generation(){
    for(int i=0; i<env.M; i++){
        for(int j=0; j<env.N; j++){
            switch (env.temp_board[i][j].type) {
                case ROCK:
                    continue;
                case FOX:
                    if(kill_fox(i, j)) continue;
                    break;
                case EMPTY:
                    env.board[i][j] = env.temp_board[i][j];
                    continue;
                default:
                    break;
            }
            env.temp_board[i][j].moved = 0;
            ++env.temp_board[i][j].age;
            env.board[i][j] = env.temp_board[i][j];
            // env.board[i][j].type = env.temp_board[i][j].type;
            // env.board[i][j].age = env.temp_board[i][j].age;
            // env.board[i][j].starve = env.temp_board[i][j].starve;
            // env.board[i][j].moved = 0;
            // if (env.temp_board[i][j].type == ROCK || env.temp_board[i][j].type == EMPTY) continue;
            // if(env.temp_board[i][j].type == FOX) 
            //     if(env.temp_board[i][j].starve >= env.foxes_starvation) 
            //         env.temp_board[i][j] = (struct Entity) {.type = EMPTY, .age=0, .starve=0, .moved=0};

        }
    }
}

void increase_ages(){
    for(int i=0; i<env.M; i++){
        for(int j=0; j<env.N; j++){
            if (env.board[i][j].type == EMPTY || env.board[i][j].type == ROCK) continue;
            ++env.board[i][j].age;
            ++env.temp_board[i][j].age;
        }
    }
}


void run_simulation(){
    // Board -> Estatico
    // Temp - board -> Editavel

    // Save current board and increase ages
    for(int i=0; i<env.M; i++){
        for(int j=0; j<env.N; j++){
            env.temp_board[i][j] = env.board[i][j];
            if (env.temp_board[i][j].type == ROCK || env.temp_board[i][j].type == EMPTY) continue;
            env.board[i][j].age++;
            env.temp_board[i][j].age++;
        }
    }

    for(int gen = 0; gen<env.generations; ++gen){
        // printf("--------------\n");
        // printf("Generation %d      \n", gen);
        // printf("--------------\n");
        //
        //
        // printf("\033[0;34m");
        // print_temp_board();
        // printf("\033[0m");

        // Update red
        for(int i=0; i<env.M; i++){
            for(int j=i%2; j<env.N; j+=2){
                move_entity(i, j);
            }
        }
        // printf("Updated red\n");
        // printf("\033[0;31m");
        // print_temp_board();
        // printf("\033[0m");

        for(int i=0; i<env.M; i++){
            for(int j=0; j<env.N; j++)
                env.board[i][j] = env.temp_board[i][j];
        }
        // Update black
        for(int i=0; i<env.M; i++){
            for(int j=(i+1)%2; j<env.N; j+=2){
                move_entity(i, j);
            }
        }
        // printf("Updated black\n");
        // print_temp_board();

        reset_generation();
    }
    print_results();
}


int main (int argc, char *argv[])
{
    if (argc != 11) {
        printf("Incorrect number of parameters\n");
        return -1;
    }
    generate_world(argv);

    double exec_time = -omp_get_wtime();

    run_simulation();

    exec_time  += omp_get_wtime();
    fprintf(stderr, "%.8fs\n", exec_time);

    //Free memory
    for (int i = 0; i < env.M; i++) {
        free(env.board[i]);
        free(env.temp_board[i]);
    }
    free(env.board);
    free(env.temp_board);

    return 0;
}
