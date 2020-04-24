/**
 * @file    mpi_jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements MPI functions for distributing vectors and matrixes,
 *          parallel distributed matrix-vector multiplication and Jacobi's
 *          method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "mpi_jacobi.h"
#include "jacobi.h"
#include "utils.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>

/*
 * TODO: Implement your solutions here
 */


void distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm)
{
    // TODO
    //Set Variable before processing
    int cordas[2], dimens[2], timeslots[2];
    int rank, tag, localrank;
    int receivecnt;
    int restdimens[2] = {0, 0};
    int *sendcnt = NULL, *displays = NULL;
    MPI_Cart_get(comm, 2, dimens, timeslots, cordas);
    MPI_Comm_rank(comm, &rank);


    MPI_Comm firstcolm;
    //int temp = dimens[0];
    if(cordas[1] != 0){
        tag = 1;
    }else{
        tag = 0;
    }
    int pval, qval; //get value of p and q
    MPI_Comm_size(comm, &pval);
    qval = (int) sqrt(pval);

    MPI_Comm_split(comm, tag, 1, &firstcolm);
    MPI_Cart_rank(comm, restdimens, &localrank);

    if(rank == localrank){
        sendcnt = new int[dimens[0]]; 
        displays = new int[dimens[0]];
        for(int i = 0; i < dimens[0]; ++i){
            sendcnt[i] = block_decompose(n, dimens[0], i);
            displays[i] = (i == 0) ? 0 : displays[i - 1] + sendcnt[i - 1];   //set value of display[i] 0 or new
        }
    }

    if(tag == 0){
        receivecnt = block_decompose(n, dimens[0], cordas[0]);
        double *receivebuffer = new double[receivecnt];
        MPI_Scatterv(&input_vector[0], sendcnt, displays, MPI_DOUBLE, receivebuffer, receivecnt, MPI_DOUBLE, localrank, firstcolm); //scatterv
        *local_vector = receivebuffer;
        //free(receivebuffer);
    }
    free(sendcnt);
    free(displays);

    //MPI_Comm_free(&firstcolm); //Not workable
    return;
}


// gather the local vector distributed among (i,0) to the processor (0,0)
void gather_vector(const int n, double* local_vector, double* output_vector, MPI_Comm comm)
{
    // TODO
    int cordas[2], dimens[2], timeslots[2];
    int sendcnt;
    int *receivecnt = NULL, *displays = NULL;
    //int restdimens[2] = {true, false};
    int restdimens[2] = {0, 0};
    MPI_Cart_get(comm, 2, dimens, timeslots, cordas);

    MPI_Comm column_comm;
    //int temp = dimens[0];
    sendcnt = block_decompose(n, dimens[0], cordas[0]);

    MPI_Comm_split(comm, cordas[1], cordas[0], &column_comm);

    int pval2, qval2; //get value of p and q
    MPI_Comm_size(comm, &pval2);
    qval2 = (int) sqrt(pval2);

    if(cordas[0] == 0 && cordas[1] == 0){
        receivecnt = new int[dimens[0]];
        displays = new int[dimens[0]];
        for(int i = 0; i < dimens[0]; ++i){
            receivecnt[i] = block_decompose(n, dimens[0], i);
            displays[i] = (i == 0) ? 0 : displays[i - 1] + receivecnt[i - 1];
        }
    }

    if(cordas[1] == 0){
        MPI_Gatherv(local_vector, sendcnt, MPI_DOUBLE, output_vector, receivecnt, displays, MPI_DOUBLE, 0, column_comm);
    }
    free(receivecnt);
    free(displays);

    MPI_Comm_free(&column_comm);
    return;
}

void distribute_matrix(const int n, double* input_matrix, double** local_matrix, MPI_Comm comm)
{
    // TODO
    int cordas[2], dimens[2], timeslots[2];
    int rank, localrank, tag;
    int rowreceivecnt, colreceivecnt;
    int *rowsendcnt = NULL, *rowdisplays = NULL;
    int restdimens[2] = {0, 0};
    MPI_Cart_get(comm, 2, dimens, timeslots, cordas);
    MPI_Comm_rank(comm, &rank);

    rowreceivecnt = block_decompose(n, dimens[0], cordas[0]);
    colreceivecnt = block_decompose(n, dimens[0], cordas[1]);

    MPI_Cart_rank(comm, restdimens, &localrank);

    int pval3, qval3; //get value of p and q
    MPI_Comm_size(comm, &pval3);
    qval3 = (int) sqrt(pval3);

    if(localrank == rank){
        rowsendcnt = new int[dimens[0] * dimens[0]];
        rowdisplays = new int[dimens[0] * dimens[0]];
        for(int i = 0; i < dimens[0]; ++i){
            for(int j = 0; j < dimens[0]; ++j){
                rowsendcnt[i*(dimens[0]) + j] = block_decompose(n, dimens[0], j);
                rowdisplays[i*(dimens[0]) + j] = ((j == 0) ? ((i == 0) ? 0 : rowdisplays[(i-1)*(dimens[0]) + j] + n*block_decompose(n, dimens[0], i-1)) : rowdisplays[i*dimens[0] + j - 1] + rowsendcnt[i*dimens[0] + j - 1]);
            }
        }
    }

    double *receivebuffer = new double[rowreceivecnt * colreceivecnt];
    int temp1 = n / dimens[0], temp2 = n % dimens[0];
    for(int i = 0; i < temp1; ++i){
        MPI_Scatterv(&input_matrix[i*n], rowsendcnt, rowdisplays, MPI_DOUBLE, &receivebuffer[i*colreceivecnt], colreceivecnt, MPI_DOUBLE, localrank, comm);
    }

    if(temp2){
        MPI_Comm row_comm;
        // if(cordas[0] >= temp2){
        //     tag = 1;
        // }else{
        //     tag = 0;
        // }
        tag = (cordas[0] >= temp2) ? 1 : 0;
        MPI_Comm_split(comm, tag, 1, &row_comm);

        if(tag == 0){
            MPI_Scatterv(&input_matrix[n*temp1], rowsendcnt, rowdisplays, MPI_DOUBLE, &receivebuffer[colreceivecnt*temp1], colreceivecnt, MPI_DOUBLE, localrank, row_comm);
        }
        MPI_Comm_free(&row_comm);
    }


    *local_matrix = receivebuffer;

    free(rowsendcnt);
    free(rowdisplays);
    //free(receivebuffer);
    return;
}


void transpose_bcast_vector(const int n, double* col_vector, double* row_vector, MPI_Comm comm)
{
    // TODO
    int d[2], loop[2], mesh_location[2], m;
    int rowcount, num_cols;
    MPI_Cart_get(comm, 2, d, loop, mesh_location);

    m = d[0];
    rowcount = block_decompose(n, m, mesh_location[0]);
    num_cols = block_decompose(n, m, mesh_location[1]); 
    
    if (mesh_location[0] == 0 && mesh_location[1] == 0)
    {
        for (int i = 0; i < rowcount; ++i) {
            row_vector[i] = col_vector[i];
        }
    }else if (mesh_location[1] == 0){
        int locationDIAG[] = {mesh_location[0], mesh_location[0]};
        int rankDcopy, rankD; 
        MPI_Cart_rank(comm, locationDIAG, &rankD);
        MPI_Send(&col_vector[0], rowcount, MPI_DOUBLE, rankD, 1, comm);
        rankDcopy = rankD;
    }else if (mesh_location[0] == mesh_location[1]){
        int original_coordinates[] = {mesh_location[0], 0};
        int rankBegin; 
        MPI_Cart_rank(comm, original_coordinates, &rankBegin);
        MPI_Recv(&row_vector[0], rowcount, MPI_DOUBLE, rankBegin, 1, comm, MPI_STATUS_IGNORE);
    }

    MPI_Comm column_comm;
    MPI_Comm_split(comm, mesh_location[1], mesh_location[0], &column_comm);
    MPI_Bcast(row_vector, num_cols, MPI_DOUBLE, mesh_location[1], column_comm);
    
    // FREE COMM
    MPI_Comm_free(&column_comm);
    return;

}


void distributed_matrix_vector_mult(const int n, double* local_A, double* local_x, double* local_y, MPI_Comm comm)
{
    // TODO
    MPI_Comm COM_ROW;
    int d[2], loop[2], mesh_location[2], m;
    int rowcount, num_cols;
    double *transposed_x = NULL, *result = NULL;
    MPI_Cart_get(comm, 2, d, loop, mesh_location);
    MPI_Comm_split(comm, mesh_location[0], mesh_location[1], &COM_ROW);

    m = d[0]; 
    rowcount = block_decompose(n, m, mesh_location[0]);
    num_cols = block_decompose(n, m, mesh_location[1]);
    transposed_x = new double[num_cols]; 
    result = new double[rowcount];
    transpose_bcast_vector(n, local_x, transposed_x, comm);

    for (int i = 0; i < rowcount; ++i)
    {
        result[i] = 0;
        for (int j = 0; j < num_cols; ++j){
            result[i] =result[i] + local_A[i*num_cols + j] * transposed_x[j];
        }
    }

    MPI_Reduce(result, local_y, rowcount, MPI_DOUBLE, MPI_SUM, 0, COM_ROW);
    return;
}

// Solves Ax = b using the iterative jacobi method
void distributed_jacobi(const int n, double* local_A, double* local_b, double* local_x,
                MPI_Comm comm, int max_iter, double l2_termination)
{
    // TODO
    int cordas[2], dimens[2], timeslots[2];
    int localrank, rowcnt, colcnt;
    double errsum, errlocal;
    int restdimens[2] = {0, 0};
    bool status = false;
    MPI_Cart_get(comm, 2, dimens, timeslots, cordas);
    MPI_Cart_rank(comm, restdimens, &localrank);

    MPI_Comm row_comm, column_comm;

    rowcnt = block_decompose(n, dimens[0], cordas[0]);
    colcnt = block_decompose(n, dimens[0], cordas[1]);

    double* R = new double[rowcnt*colcnt];

    for(int i = 0; i < rowcnt; ++i){
        for(int j = 0; j < colcnt; ++j){
            R[i*colcnt + j] = (cordas[0] == cordas[1] && i == j) ? 0 : local_A[i*colcnt + j];
        }
    }

    MPI_Comm_split(comm, cordas[0], cordas[1], &row_comm);
    MPI_Comm_split(comm, cordas[1], cordas[0], &column_comm);

    double *temp = new double[rowcnt];
    double *Diag = NULL;
    double *Rsum = NULL;
    double *Asum = NULL;

    if(cordas[1] == 0){
        Diag = new double[rowcnt];
    }

    for(int i = 0; i < rowcnt; ++i){
        temp[i] = (cordas[0] != cordas[1]) ? 0.0 : local_A[i*colcnt + i];
    }

    MPI_Reduce(temp, Diag, rowcnt, MPI_DOUBLE, MPI_SUM, 0, row_comm);

    for(int i = 0; i < rowcnt; ++i){
        local_x[i] = 0.0;
    }

    if(cordas[1] == 0){
        Rsum = new double[rowcnt];
        Asum = new double[rowcnt];
    }

    for(int i = 0; i < max_iter; ++i){
        distributed_matrix_vector_mult(n, R, local_x, Rsum, comm);
        distributed_matrix_vector_mult(n, local_A, local_x, Asum, comm);

        if(cordas[1] == 0){
            errsum = 0.0;
            errlocal = 0.0;
            for(int j = 0; j < rowcnt; ++j){
                errlocal += (Asum[j] - local_b[j]) * (Asum[j] - local_b[j]);
            }
            MPI_Reduce(&errlocal, &errsum, 1, MPI_DOUBLE, MPI_SUM, 0, column_comm);
            if(cordas[0] == 0 && errsum < l2_termination){
                status = true;
            }
        }

        MPI_Bcast(&status, 1, MPI::BOOL, localrank, comm);

        if(status){
            break;
        }else if(cordas[1] == 0){
            for(int p = 0; p < rowcnt; ++p){
                local_x[p] = (local_b[p] - Rsum[p])/Diag[p];
            }
        }
    }

    free(R);
    free(temp);
    free(Diag);
    free(Rsum);
    free(Asum);

    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&column_comm);

    return;
}


// wraps the distributed matrix vector multiplication
void mpi_matrix_vector_mult(const int n, double* A,
                            double* x, double* y, MPI_Comm comm)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_x = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &x[0], &local_x, comm);

    // allocate local result space
    double* local_y = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_matrix_vector_mult(n, local_A, local_x, local_y, comm);

    // gather results back to rank 0
    gather_vector(n, local_y, y, comm);
}

// wraps the distributed jacobi function
void mpi_jacobi(const int n, double* A, double* b, double* x, MPI_Comm comm,
                int max_iter, double l2_termination)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_b = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &b[0], &local_b, comm);

    // allocate local result space
    double* local_x = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_jacobi(n, local_A, local_b, local_x, comm, max_iter, l2_termination);

    // gather results back to rank 0
    gather_vector(n, local_x, x, comm);
}
