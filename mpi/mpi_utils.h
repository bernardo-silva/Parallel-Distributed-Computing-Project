#ifndef MPIUTILS
#define MPIUTILS
#include <mpi.h>
// #include "environment.c"

typedef struct Entity Entity;
typedef struct Environment Environment;

void recieve_from_neigs(Environment* env, 
                        Entity* above, int tag_above, 
                        Entity* below, int tag_below, 
                        Entity* left, int tag_left, 
                        Entity* right, int tag_right, 
                        MPI_Request* requests,
                        int idx_above, int idx_below, int idx_left, int idx_right
                        );

void send_to_neigs(Environment* env, 
                        Entity* above, int tag_above, 
                        Entity* below, int tag_below, 
                        Entity* left, int tag_left, 
                        Entity* right, int tag_right, 
                        MPI_Request* requests,
                        int idx_above, int idx_below, int idx_left, int idx_right
                        );

void empty_ghost_lines(Environment* env);

void column_to_row(Environment* env, int col_idx, Entity* row);

void columns_to_rows(Environment* env, Entity* col_left, Entity* col_right);

void row_to_column(Environment* env, int col_idx, Entity* row);

void rows_to_columns(Environment* env, Entity* col_left, Entity* col_right);
#endif // !MPIUTILS
