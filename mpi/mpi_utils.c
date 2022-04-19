#include "mpi_utils.h"
#include "environment.h"
#include <mpi.h>

#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

extern MPI_Datatype MPI_Entity;

void recieve_from_neigs(Environment* env, 
                        Entity* above, int tag_above, 
                        Entity* below, int tag_below, 
                        Entity* left, int tag_left, 
                        Entity* right, int tag_right, 
                        MPI_Request* requests,
                        int idx_above, int idx_below, int idx_left, int idx_right
                        ){
    // MPI request ghost line from process above
    if(env->is_not_top)
        MPI_Irecv(above, env->column_block_size_ghost, MPI_Entity,
                  env->neigs_ids[UP], tag_above,
                  env->cart_comm, requests + idx_above);

    // MPI request ghost line from process below
    if(env->is_not_bottom)
        MPI_Irecv(below, env->column_block_size_ghost, MPI_Entity,
                  env->neigs_ids[DOWN], tag_below,
                  env->cart_comm, requests + idx_below);

    // MPI request ghost right from process to the left
    if(env->is_not_left)
        MPI_Irecv(left, env->row_block_size_ghost, MPI_Entity,
                  env->neigs_ids[LEFT], tag_left,
                  env->cart_comm, requests + idx_left);

    // MPI request ghost column from process to the right
    if(env->is_not_right)
        MPI_Irecv(right, env->row_block_size_ghost, MPI_Entity,
                  env->neigs_ids[RIGHT], tag_right,
                  env->cart_comm, requests + idx_right);
}

void send_to_neigs(Environment* env, 
                        Entity* above, int tag_above, 
                        Entity* below, int tag_below, 
                        Entity* left, int tag_left, 
                        Entity* right, int tag_right, 
                        MPI_Request* requests,
                        int idx_above, int idx_below, int idx_left, int idx_right
                        ){
    if(env->is_not_top)
        MPI_Send(above, env->column_block_size_ghost, MPI_Entity,
                  env->neigs_ids[UP], tag_above, env->cart_comm) ;
                  // requests + idx_above);
    
    if(env->is_not_bottom)
        MPI_Send(below, env->column_block_size_ghost, MPI_Entity,
                  env->neigs_ids[DOWN], tag_below, env->cart_comm);
                  // requests + idx_below);

    if(env->is_not_left)
        MPI_Send(left, env->row_block_size_ghost, MPI_Entity,
                  env->neigs_ids[LEFT], tag_left, env->cart_comm);
                  // requests + idx_left);
    
    if(env->is_not_right)
        MPI_Send(right, env->row_block_size_ghost, MPI_Entity,
                  env->neigs_ids[RIGHT], tag_right, env->cart_comm);
                  // requests + idx_right);
    
}

void empty_ghost_lines(Environment* env){
    if(env->is_not_top)
        for(int j=0; j<env->column_block_size_ghost;j++)
            env->temp_board[0][j] = 
                (Entity) {.type = EMPTY, .age=0, .starve=0, .moved=0};

    if(env->is_not_bottom)
        for(int j=0; j<env->column_block_size_ghost;j++)
            env->temp_board[env->row_block_size_ghost - 1][j] = 
                (Entity) {.type = EMPTY, .age=0, .starve=0, .moved=0};
                
    if(env->is_not_left)
        for(int j=0; j<env->row_block_size_ghost;j++)
            env->temp_board[j][0] = 
                (Entity) {.type = EMPTY, .age=0, .starve=0, .moved=0};

    if(env->is_not_right)
        for(int j=0; j<env->row_block_size_ghost;j++)
            env->temp_board[j][env->column_block_size_ghost - 1] =
                (Entity) {.type = EMPTY, .age=0, .starve=0, .moved=0};
}

void column_to_row(Environment* env, int col_idx, Entity* row){
        for(int i=0; i<env->row_block_size_ghost; i++)
            row[i] = env->temp_board[i][col_idx];
}

void columns_to_rows(Environment* env, Entity* col_left, Entity* col_right){
    if(env->is_not_left)
        column_to_row(env, 0, col_left);

    if(env->is_not_right)
        column_to_row(env, env->column_block_size_ghost-1, col_right);
}

void row_to_column(Environment* env, int col_idx, Entity* row){
        for(int i=0; i<env->row_block_size_ghost; i++)
            env->temp_board[i][col_idx] = row[i];
}

void rows_to_columns(Environment* env, Entity* col_left, Entity* col_right){
    if(env->is_not_left)
        row_to_column(env, 0, col_left);

    if(env->is_not_right)
        row_to_column(env, env->column_block_size_ghost-1, col_right);
}
