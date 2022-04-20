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
    
// #define N_GHOST_LINES(coord,dim) (!coord || id==p-1)? 1:2
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
        if (j >= env->column_low_ghost && j <= env->column_high_ghost)
	if(position_empty(env, i-env->row_low_ghost, j-env->column_low_ghost))
	    insert_animal(env, i-env->row_low_ghost, j-env->column_low_ghost, atype);
    }
}

void get_block_dimentions(Environment* env){
    env->row_low  = BLOCK_LOW(env->coords[0], env->grid_dims[0], env->M);
    env->row_high = BLOCK_HIGH(env->coords[0], env->grid_dims[0], env->M);

    env->column_low  = BLOCK_LOW(env->coords[1], env->grid_dims[1], env->N);
    env->column_high = BLOCK_HIGH(env->coords[1], env->grid_dims[1], env->N);

    env->row_block_size    = BLOCK_SIZE(env->coords[0], env->grid_dims[0], env->M);
    env->column_block_size = BLOCK_SIZE(env->coords[1], env->grid_dims[1], env->N);

    env->n_ghost_rows    = env->is_not_top + env->is_not_bottom;
    env->n_ghost_columns = env->is_not_left + env->is_not_right;

    env->row_low_ghost  = BLOCK_LOW_GHOST(env->coords[0], env->grid_dims[0], env->M);
    env->row_high_ghost = BLOCK_HIGH_GHOST(env->coords[0], env->grid_dims[0], env->M);
    env->column_low_ghost  = BLOCK_LOW_GHOST(env->coords[1], env->grid_dims[1], env->N);
    env->column_high_ghost = BLOCK_HIGH_GHOST(env->coords[1], env->grid_dims[1], env->N);

    env->row_block_size_ghost    = env->row_block_size    + env->n_ghost_rows;
    env->column_block_size_ghost = env->column_block_size + env->n_ghost_columns;
}

void get_neigs_ids(Environment* env){
    int neig_coords[2];
    neig_coords[0] = env->coords[0];
    neig_coords[1] = env->coords[1];

    if(env->is_not_top){
        neig_coords[0] = env->coords[0]-1;
        neig_coords[1] = env->coords[1];
        MPI_Cart_rank(env->cart_comm, neig_coords, &env->neigs_ids[0]);
    }
    if(env->is_not_bottom){
        neig_coords[0] = env->coords[0]+1;
        neig_coords[1] = env->coords[1];
        MPI_Cart_rank(env->cart_comm, neig_coords, &env->neigs_ids[1]);
    }
    if(env->is_not_left){
        neig_coords[0] = env->coords[0];
        neig_coords[1] = env->coords[1]-1;
        MPI_Cart_rank(env->cart_comm, neig_coords, &env->neigs_ids[2]);
    }
    if(env->is_not_right){
        neig_coords[0] = env->coords[0];
        neig_coords[1] = env->coords[1]+1;
        MPI_Cart_rank(env->cart_comm, neig_coords, &env->neigs_ids[3]);
    }
}

void create_datatypes(Environment* env){
    extern MPI_Datatype MPI_Entity;

    const int nfields = 1;
    const int block_lens[] = {4};
    MPI_Aint displacements = offsetof(Entity, type);
    MPI_Datatype types[] = {MPI_INT}; 
    MPI_Type_create_struct(nfields, block_lens, &displacements, types, &MPI_Entity);
    MPI_Type_commit(&MPI_Entity);

    // extern MPI_Datatype MPI_Column;
    // MPI_Type_vector(env->row_block_size_ghost, 1, env->column_block_size_ghost,
    //                 MPI_Entity, &MPI_Column);
    // MPI_Type_commit(&MPI_Column);
}



void generate_world(Environment* env, char *argv[], int id, int p, int*
                    grid_dims, MPI_Comm cart_comm){
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

    env->cart_comm = cart_comm;
    env->id = id;
    env->p = p;
    env->grid_dims[0] = grid_dims[0];
    env->grid_dims[1] = grid_dims[1];
    MPI_Cart_coords(cart_comm, id, 2, env->coords);

    
    env->is_not_top = (!env->coords[0])?0:1;
    env->is_not_bottom = (env->coords[0] == env->grid_dims[0]-1)?0:1;
    env->is_not_left = (!env->coords[1])?0:1;
    env->is_not_right = (env->coords[1] == env->grid_dims[1]-1)?0:1;

    env->n_neigs = env->is_not_top + env->is_not_bottom + env->is_not_left +
        env->is_not_right;

    get_block_dimentions(env);

    get_neigs_ids(env);
    create_datatypes(env);

    env->board = (Entity**) malloc(sizeof(Entity*) * env->row_block_size_ghost);
    env->temp_board = (Entity**) malloc(sizeof(Entity*) * env->row_block_size_ghost);

    for (int i=0;i<env->row_block_size_ghost;i++){
        env->board[i] = (Entity*) malloc(env->column_block_size_ghost * sizeof(Entity));
        env->temp_board[i] = (Entity*) malloc(env->column_block_size_ghost * sizeof(Entity));
    }


    for (int i=0;i<env->row_block_size_ghost;i++){
        for (int j=0;j<env->column_block_size_ghost;j++){
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
    for(int j=0; j<env->column_block_size_ghost; j++)
        printf("%02d|",j);

    printf("\n");

    for(int i=0; i<env->row_block_size_ghost; i++){
        printf("%02d:",i);
        for(int j=0; j<env->column_block_size_ghost; j++){
            printf(" %c|", animal[env->board[i][j].type]);
        }
        printf("\n");fflush(stdout);
    }
    printf("---------------------\n");fflush(stdout);
}

void print_temp_board(Environment* env){
    char * animal = " *RF";
    printf("---------------------\n");fflush(stdout);
    printf("   ");
    for(int j=0; j<env->column_block_size_ghost; j++)
        printf("%02d|",j);

    printf("\n");

    for(int i=0; i<env->row_block_size_ghost; i++){
        printf("%02d:",i);
        for(int j=0; j<env->column_block_size_ghost; j++){
            printf(" %c|", animal[env->temp_board[i][j].type]);
        }
        printf("\n");fflush(stdout);
    }
    printf("---------------------\n");fflush(stdout);
}

void get_results(Environment* env, int* rocks, int* rabbits, int* foxes){
    *rocks=0;
    *rabbits=0;
    *foxes=0;
    for(int i=env->is_not_top; i<env->row_block_size + env->is_not_top; i++){
        for(int j=env->is_not_left; j<env->column_block_size + env->is_not_left; j++){
            if(env->board[i][j].type == ROCK){ (*rocks)++; continue; }
            if(env->board[i][j].type == RABBIT){ (*rabbits)++; continue; }
            if(env->board[i][j].type == FOX){ (*foxes)++; continue; }
        }
    }
    // printf("%d %d %d\n", rocks, rabbits, foxes);
}
