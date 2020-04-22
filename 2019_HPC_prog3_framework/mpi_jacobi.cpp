// /**
//  * @file    mpi_jacobi.cpp
//  * @author  Patrick Flick <patrick.flick@gmail.com>
//  * @brief   Implements MPI functions for distributing vectors and matrixes,
//  *          parallel distributed matrix-vector multiplication and Jacobi's
//  *          method.
//  *
//  * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
//  */

// #include "mpi_jacobi.h"
// #include "jacobi.h"
// #include "utils.h"

// #include <stdlib.h>
// #include <assert.h>
// #include <math.h>
// #include <vector>

// /*
//  * TODO: Implement your solutions here
//  */

// void distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm)
// {
//     // TODO
//     // retrieve Cartesian topology information
//     int d[2];
//     int loop[2];
//     int mesh_location[2];
//     MPI_Cart_get(comm, 2, d, loop, mesh_location);
//     int m = d[0];
//     int rank; 
//     MPI_Comm_rank(comm, &rank);

//     //create comm for the first column
//     int color;
//     if (mesh_location[1]==0){
//         color = 0;
//     }
//     else{
//         color = 1;
//     }
//     MPI_Comm first_column;
//     MPI_Comm_split(comm, color, 1, &first_column);

//     //calculate the parameters for scatter
//     int rankzero; 
//     int firstgridcoords[2] = {0,0};
//     MPI_Cart_rank(comm, firstgridcoords, &rankzero);
//     int *sendcounts = NULL;
//     int *displs = NULL;
//     if (rank == rankzero)
//     {
//         sendcounts = new int[m];
//         displs = new int[m];
//         for (int i = 0; i < m; i++)
//         {
//             sendcounts[i] = block_decompose(n, m, i);
//             if (i == 0){
//                 displs[i] = 0;
//             }
//             else {
//                 displs[i] = displs[i-1] + sendcounts[i-1];
//             }
//         }
//     }

//     //send data from processor(0,0) onto processors(i,0)
//     if (color == 0)
//     {
//         int recvcount = block_decompose(n, m, mesh_location[0]);
//         double *recvbuf = new double[recvcount];
//         MPI_Scatterv(&input_vector[0], sendcounts, displs, MPI_DOUBLE, 
//         recvbuf, recvcount, MPI_DOUBLE, rankzero, first_column);
//         *local_vector = recvbuf; 
//     }
// }

// // gather the local vector distributed among (i,0) to the processor (0,0)
// void gather_vector(const int n, double* local_vector, double* output_vector, MPI_Comm comm)
// {
//     // TODO
//     // retrieve Cartesian topology information
//     int d[2];
//     int loop[2];
//     int mesh_location[2];
//     MPI_Cart_get(comm, 2, d, loop, mesh_location);

//     // get sendcount for gathering
//     int m = d[0];
//     int sendcount = block_decompose(n, m, mesh_location[0]);

//     //create comm for the column
//     MPI_Comm column_comm;
//     MPI_Comm_split(comm, mesh_location[1], mesh_location[0], &column_comm);

//     //calculate the parameters for gathering
//     int *recvcounts = NULL;
//     int *displs = NULL;
//     if (mesh_location[0]==0 && mesh_location[1]==0)
//     {
//         recvcounts = new int[m];
//         displs = new int[m];
//         for (int i = 0; i < m; i++)
//         {
//             recvcounts[i] = block_decompose(n, m, i);
//             if (i == 0){
//                 displs[i] = 0;
//             }
//             else {
//                 displs[i] = displs[i-1] + recvcounts[i-1];
//             }
//         }
//     }
    
//     //Gather data from processors (i,0) onto the processor (0,0)
//     if(mesh_location[1]==0){
//         MPI_Gatherv(local_vector, sendcount, MPI_DOUBLE, output_vector, 
//         recvcounts, displs, MPI_DOUBLE, 0, column_comm);
//     }
//     // release the column_comm
//     MPI_Comm_free(&column_comm);
// }

// void distribute_matrix(const int n, double* input_matrix, double** local_matrix, MPI_Comm comm)
// {
//     // TODO
//     // retrieve Cartesian topology information
//     int d[2];
//     int loop[2];
//     int mesh_location[2];
//     MPI_Cart_get(comm, 2, d, loop, mesh_location);
//     int rank; 
//     MPI_Comm_rank(comm, &rank);
    
//     //calculate the parameters for distribution
//     //number of rows or cols sent to other processors
//     int m = d[0]; 
//     int recvcountrow = block_decompose(n, m, mesh_location[0]);
//     int recvcountcol = block_decompose(n, m, mesh_location[1]);

//     //setup parameters for distribution
//     int firstgridcoords[2] = {0,0};
//     int rankzero; 
//     MPI_Cart_rank(comm, firstgridcoords, &rankzero);
//     int *rowsendcount = NULL;
//     int *rowdispls = NULL;

//     if (rank == rankzero)
//     {
//         //number of elements sent to each processor
//         rowsendcount = new int[m*m]; 
//         //displacement relative to sendbuf
//         rowdispls = new int[m*m]; 
//         for (int i = 0; i < m; i++)
//         {
//             for (int j = 0; j < m; j++)
//             {
//                 rowsendcount[i*m+j] = block_decompose(n, m, j); 
//                 if (j == 0) {
//                     if(i==0){
//                         rowdispls[i*m+j] = 0;
//                     }
//                     else {
//                         rowdispls[i*m+j] = rowdispls[(i-1)*m+j] + n*block_decompose(n, m, i-1);
//                     }
//                 }
//                 else {
//                     rowdispls[i*m+j] = rowdispls[i*m+j-1] + rowsendcount[i*m+j-1];
//                 }
//             }
//         }
        
//     }
//     //create recvbuf
//     double *recvbuf = new double[recvcountrow*recvcountcol];

//     //scatter data to other processors
//     for (int i = 0; i < n / m; i++) {
//         MPI_Scatterv(&input_matrix[n*i], rowsendcount, rowdispls, MPI_DOUBLE, 
//         &recvbuf[recvcountcol*i], recvcountcol, MPI_DOUBLE, rankzero, comm);
//     }
    
//     //scatter the remainder data
//     if (n%m != 0) 
//     {
//         //create new communicators
//         MPI_Comm rows_comm;
//         int color;
//         if (mesh_location[0] < n%m){
//             color = 0;
//         }
//         else {
//             color = 1;
//         }
//         MPI_Comm_split(comm, color, 1, &rows_comm);

//         int temp = n / m;
//         if(color==0) {
//             MPI_Scatterv(&input_matrix[n*temp], rowsendcount, rowdispls, MPI_DOUBLE, 
//             &recvbuf[recvcountcol*temp], recvcountcol, MPI_DOUBLE, rankzero, rows_comm);
//         }
//         //release the comm
//         MPI_Comm_free(&rows_comm);
//     }

//     *local_matrix = recvbuf;
// }

// void transpose_bcast_vector(const int n, double* col_vector, double* row_vector, MPI_Comm comm)
// {
//     int d[2];
//     int loop[2];
//     int mesh_location[2];
//     MPI_Cart_get(comm, 2, d, loop, mesh_location);

//     // get # of col and row
//     int m = d[0];
//     int rowcount = block_decompose(n, m, mesh_location[0]);
//     int num_cols = block_decompose(n, m, mesh_location[1]); 
    
//     if (mesh_location[0] == 0 && mesh_location[1] == 0)
//     {
//         for (int i = 0; i < rowcount; i++) {
//             row_vector[i] = col_vector[i];
//         }
//     }else if (mesh_location[1] == 0){
    
//         int locationDIAG[] = {mesh_location[0], mesh_location[0]};
//         int rankDcopy; 
//         int rankD; 
//         MPI_Cart_rank(comm, locationDIAG, &rankD);
//         MPI_Send(&col_vector[0], rowcount, MPI_DOUBLE, rankD, 1, comm);
//         rankDcopy = rankD;
//     }else if (mesh_location[0] == mesh_location[1]){
//         int original_coordinates[] = {mesh_location[0], 0};
//         int rankBegin; 
//         MPI_Cart_rank(comm, original_coordinates, &rankBegin);
//         MPI_Recv(&row_vector[0], rowcount, MPI_DOUBLE, rankBegin, 1, comm, MPI_STATUS_IGNORE);
//     }
//     MPI_Comm column_comm;
//     MPI_Comm_split(comm, mesh_location[1], mesh_location[0], &column_comm);
//     MPI_Bcast(row_vector, num_cols, MPI_DOUBLE, mesh_location[1], column_comm);
    
//     // FREE COMM
//     MPI_Comm_free(&column_comm);
// }

// void distributed_matrix_vector_mult(const int n, double* local_A, double* local_x, double* local_y, MPI_Comm comm)
// {
//     MPI_Comm COM_ROW;
//     int d[2];
//     int loop[2];
//     int mesh_location[2];
//     MPI_Cart_get(comm, 2, d, loop, mesh_location);
//     MPI_Comm_split(comm, mesh_location[0], mesh_location[1], &COM_ROW);


//     int m = d[0]; 
//     int rowcount = block_decompose(n, m, mesh_location[0]);


//     int num_cols = block_decompose(n, m, mesh_location[1]);

//     double *transposed_x = new double[num_cols]; 
//     transpose_bcast_vector(n, local_x, transposed_x, comm);


//     double *result = new double[rowcount];

//     for (int i = 0; i < rowcount; i++)
//     {
//         result[i] = 0;
//         for (int j = 0; j < num_cols; j++){
//             result[i] =result[i] + local_A[i*num_cols + j] * transposed_x[j];
//         }
//     }

//     MPI_Reduce(result, local_y, rowcount, MPI_DOUBLE, MPI_SUM, 0, COM_ROW);
// }

// // Solves Ax = b using the iterative jacobi method
// void distributed_jacobi(const int n, double* local_A, double* local_b, double* local_x, MPI_Comm comm, int max_iter, double l2_termination)
// {
//     // TODO
//     // retrieve Cartesian topology information
//     int d[2];
//     int loop[2];
//     int mesh_location[2];
//     MPI_Cart_get(comm, 2, d, loop, mesh_location);

//     // rank of coordinate_zero
//     int coordinate_zero[] = {0,0};
//     int rankzero; 
//     MPI_Cart_rank(comm, coordinate_zero, &rankzero);

//     // number of rows or columns
//     int m = d[0];
//     int rowcount = block_decompose(n, m, mesh_location[0]); 
//     int num_cols = block_decompose(n, m, mesh_location[1]); 

//     // initialize R, diagonal elements as 0, others as A(i,j)
//     double* R = new double[rowcount*num_cols];
//     for (int i = 0; i < rowcount; i++){
//         for (int j = 0; j < num_cols; j++){
//             if (mesh_location[0]==mesh_location[1] && i==j){
//                 R[i*num_cols + j] = 0;
//             }
//             else {
//                 R[i*num_cols + j] = local_A[i*num_cols + j];
//             }
//         }
//     }
        
//     // calculate matrix D = diag(A)
//     MPI_Comm COM_ROW;
//     MPI_Comm_split(comm, mesh_location[0], mesh_location[1], &COM_ROW);
//     MPI_Comm column_comm;
//     MPI_Comm_split(comm, mesh_location[1], mesh_location[0], &column_comm);
//     double *temp = new double[rowcount];
//     for (int i = 0; i < rowcount; i++) {
//         if (mesh_location[0]==mesh_location[1]){
//             temp[i] = local_A[i*num_cols + i];
//         }
//         else {
//             temp[i] = 0.0;
//         }
//     }
//     double *D_diag = NULL;
//     if (mesh_location[1] == 0) {
//         D_diag = new double[rowcount];
//     }
//     MPI_Reduce(temp, D_diag, rowcount, MPI_DOUBLE, MPI_SUM, 0, COM_ROW);

//     // initialize local_x
//     for (int i = 0; i < rowcount; i++){
//         local_x[i] = 0.0;
//     }

//     // product of R*x at first column
//     double *sum_Rx = NULL; 
//     if (mesh_location[1] == 0){
//         sum_Rx = new double[rowcount];
//     }

//     // product of A*x at first column
//     double *sum_Ax = NULL;
//     if (mesh_location[1] == 0){
//         sum_Ax = new double[rowcount];
//     }

//     // iteration check, continue or stop
//     bool iteration_status = false;

//     // iterative calcualtion of x
//     for (int iter = 0; iter < max_iter; iter++)
//     {
//         // calculate product of R*x at first column
//         distributed_matrix_vector_mult(n, R, local_x, sum_Rx, comm);

//         // calculate product of A*x at first column
//         distributed_matrix_vector_mult(n, local_A, local_x, sum_Ax, comm);
        
//         // compare the error with l2_termination
//         if (mesh_location[1] == 0)
//         {
//             double sum_error = 0;
//             double local_error = 0;
//             for (int i = 0; i < rowcount; i++) {
//                 local_error += (sum_Ax[i]-local_b[i])*(sum_Ax[i]-local_b[i]);
//             }
//             // reduce/sum local_error to processor 0
//             MPI_Reduce(&local_error, &sum_error, 1, MPI_DOUBLE, MPI_SUM, 0, column_comm);

//             if (mesh_location[0] == 0 && sum_error < l2_termination)
//             {
//                 iteration_status = true;
//             }
//         }

//         // brocast result of iteration check to other processors
//         MPI_Bcast(&iteration_status, 1, MPI::BOOL, rankzero, comm);

//         // iteration check, stop or continue
//         if (iteration_status)
//         {
//             break;
//         }
//         else if (mesh_location[1] == 0) 
//         {
//             // update local_x with D^(-1)*(b-R*x)
//             for (int i = 0; i < rowcount; i++) {
//                 local_x[i] = (local_b[i]-sum_Rx[i]) / D_diag[i];
//             }
//         }
//     }
// }

// // wraps the distributed matrix vector multiplication
// void mpi_matrix_vector_mult(const int n, double* A,
//                             double* x, double* y, MPI_Comm comm)
// {
//     // distribute the array onto local processors!
//     double* local_A = NULL;
//     double* local_x = NULL;
//     distribute_matrix(n, &A[0], &local_A, comm);
//     distribute_vector(n, &x[0], &local_x, comm);

//     // allocate local result space
//     double* local_y = new double[block_decompose_by_dim(n, comm, 0)];
//     distributed_matrix_vector_mult(n, local_A, local_x, local_y, comm);

//     // gather results back to rank 0
//     gather_vector(n, local_y, y, comm);
// }

// // wraps the distributed jacobi function
// void mpi_jacobi(const int n, double* A, double* b, double* x, MPI_Comm comm,
//                 int max_iter, double l2_termination)
// {
//     // distribute the array onto local processors!
//     double* local_A = NULL;
//     double* local_b = NULL;
//     distribute_matrix(n, &A[0], &local_A, comm);
//     distribute_vector(n, &b[0], &local_b, comm);

//     // allocate local result space
//     double* local_x = new double[block_decompose_by_dim(n, comm, 0)];
//     distributed_jacobi(n, local_A, local_b, local_x, comm, max_iter, l2_termination);

//     // gather results back to rank 0
//     gather_vector(n, local_x, x, comm);
// }

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
    //retrieve Cartesian topology information
    // int dimension[2];
    // int periods[2];
    // int coordinates[2];
    // MPI_Cart_get(comm, 2, dimension, periods, coordinates);
    // int m = dimension[0];
    // int rank; 
    // MPI_Comm_rank(comm, &rank);

    // //create comm for the first column
    // int color;
    // if (coordinates[1]==0){
    //     color = 0;
    // }
    // else{
    //     color = 1;
    // }
    // MPI_Comm first_column;
    // MPI_Comm_split(comm, color, 1, &first_column);

    // //calculate the parameters for scatter
    // int rankzero; 
    // int firstgridcoords[2] = {0,0};
    // MPI_Cart_rank(comm, firstgridcoords, &rankzero);
    // int *sendcounts = NULL;
    // int *displs = NULL;
    // if (rank == rankzero)
    // {
    //     sendcounts = new int[m];
    //     displs = new int[m];
    //     for (int i = 0; i < m; i++)
    //     {
    //         sendcounts[i] = block_decompose(n, m, i);
    //         if (i == 0){
    //             displs[i] = 0;
    //         }
    //         else {
    //             displs[i] = displs[i-1] + sendcounts[i-1];
    //         }
    //     }
    // }

    // //send data from processor(0,0) onto processors(i,0)
    // if (color == 0)
    // {
    //     int recvcount = block_decompose(n, m, coordinates[0]);
    //     double *recvbuf = new double[recvcount];
    //     MPI_Scatterv(&input_vector[0], sendcounts, displs, MPI_DOUBLE, 
    //     recvbuf, recvcount, MPI_DOUBLE, rankzero, first_column);
    //     *local_vector = recvbuf; 
    // }

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
    MPI_Comm_split(comm, tag, 1, &firstcolm);
    MPI_Cart_rank(comm, restdimens, &localrank);

    if(rank == localrank){
        sendcnt = new int[dimens[0]];
        displays = new int[dimens[0]];
        for(int i = 0; i < dimens[0]; ++i){
            sendcnt[i] = block_decompose(n, dimens[0], i);
            displays[i] = (i == 0) ? 0 : displays[i - 1] + sendcnt[i - 1];   
        }
    }

    if(tag == 0){
        receivecnt = block_decompose(n, dimens[0], cordas[0]);
        double *receivebuffer = new double[receivecnt];
        MPI_Scatterv(&input_vector[0], sendcnt, displays, MPI_DOUBLE, receivebuffer, receivecnt, MPI_DOUBLE, localrank, firstcolm);
        *local_vector = receivebuffer;
    }
}

// gather the local vector distributed among (i,0) to the processor (0,0)
void gather_vector(const int n, double* local_vector, double* output_vector, MPI_Comm comm)
{
    // TODO
    // retrieve Cartesian topology information
    // int dimension[2];
    // int periods[2];
    // int coordinates[2];
    // MPI_Cart_get(comm, 2, dimension, periods, coordinates);

    // // get sendcount for gathering
    // int m = dimension[0];
    // int sendcount = block_decompose(n, m, coordinates[0]);

    // //create comm for the column
    // MPI_Comm column_comm;
    // MPI_Comm_split(comm, coordinates[1], coordinates[0], &column_comm);

    // //calculate the parameters for gathering
    // int *recvcounts = NULL;
    // int *displs = NULL;
    // if (coordinates[0]==0 && coordinates[1]==0)
    // {
    //     recvcounts = new int[m];
    //     displs = new int[m];
    //     for (int i = 0; i < m; i++)
    //     {
    //         recvcounts[i] = block_decompose(n, m, i);
    //         if (i == 0){
    //             displs[i] = 0;
    //         }
    //         else {
    //             displs[i] = displs[i-1] + recvcounts[i-1];
    //         }
    //     }
    // }
    
    // //Gather data from processors (i,0) onto the processor (0,0)
    // if(coordinates[1]==0){
    //     MPI_Gatherv(local_vector, sendcount, MPI_DOUBLE, output_vector, 
    //     recvcounts, displs, MPI_DOUBLE, 0, column_comm);
    // }
    // // release the column_comm
    // MPI_Comm_free(&column_comm);

    int cordas[2], dimens[2], timeslots[2];
    int sendcnt;
    int *receivecnt = NULL, *displays = NULL;
    MPI_Cart_get(comm, 2, dimens, timeslots, cordas);

    MPI_Comm column_comm;
    //int temp = dimens[0];
    sendcnt = block_decompose(n, dimens[0], cordas[0]);

    MPI_Comm_split(comm, cordas[1], cordas[0], &column_comm);

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
    MPI_Comm_free(&column_comm);
}

void distribute_matrix(const int n, double* input_matrix, double** local_matrix, MPI_Comm comm)
{
    // TODO
    // retrieve Cartesian topology information
    // int dimension[2];
    // int periods[2];
    // int coordinates[2];
    // MPI_Cart_get(comm, 2, dimension, periods, coordinates);
    // int rank; 
    // MPI_Comm_rank(comm, &rank);
    
    // //calculate the parameters for distribution
    // //number of rows or cols sent to other processors
    // int m = dimension[0]; 
    // int recvcountrow = block_decompose(n, m, coordinates[0]);
    // int recvcountcol = block_decompose(n, m, coordinates[1]);

    // //setup parameters for distribution
    // int firstgridcoords[2] = {0,0};
    // int rankzero; 
    // MPI_Cart_rank(comm, firstgridcoords, &rankzero);
    // int *rowsendcount = NULL;
    // int *rowdispls = NULL;

    // if (rank == rankzero)
    // {
    //     //number of elements sent to each processor
    //     rowsendcount = new int[m*m]; 
    //     //displacement relative to sendbuf
    //     rowdispls = new int[m*m]; 
    //     for (int i = 0; i < m; i++)
    //     {
    //         for (int j = 0; j < m; j++)
    //         {
    //             rowsendcount[i*m+j] = block_decompose(n, m, j); 
    //             if (j == 0) {
    //                 if(i==0){
    //                     rowdispls[i*m+j] = 0;
    //                 }
    //                 else {
    //                     rowdispls[i*m+j] = rowdispls[(i-1)*m+j] + n*block_decompose(n, m, i-1);
    //                 }
    //             }
    //             else {
    //                 rowdispls[i*m+j] = rowdispls[i*m+j-1] + rowsendcount[i*m+j-1];
    //             }
    //         }
    //     }
        
    // }
    // //create recvbuf
    // double *recvbuf = new double[recvcountrow*recvcountcol];

    // //scatter data to other processors
    // for (int i = 0; i < n / m; i++) {
    //     MPI_Scatterv(&input_matrix[n*i], rowsendcount, rowdispls, MPI_DOUBLE, 
    //     &recvbuf[recvcountcol*i], recvcountcol, MPI_DOUBLE, rankzero, comm);
    // }
    
    // //scatter the remainder data
    // if (n%m != 0) 
    // {
    //     //create new communicators
    //     MPI_Comm rows_comm;
    //     int color;
    //     if (coordinates[0] < n%m){
    //         color = 0;
    //     }
    //     else {
    //         color = 1;
    //     }
    //     MPI_Comm_split(comm, color, 1, &rows_comm);

    //     int temp = n / m;
    //     if(color==0) {
    //         MPI_Scatterv(&input_matrix[n*temp], rowsendcount, rowdispls, MPI_DOUBLE, 
    //         &recvbuf[recvcountcol*temp], recvcountcol, MPI_DOUBLE, rankzero, rows_comm);
    //     }
    //     //release the comm
    //     MPI_Comm_free(&rows_comm);
    // }

    // *local_matrix = recvbuf;

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

    if(localrank == rank){
        rowsendcnt = new int[dimens[0] * dimens[0]];
        rowdisplays = new int[dimens[0] * dimens[0]];
        for(int i = 0; i < dimens[0]; ++i){
            for(int j = 0; j < dimens[0]; ++j){
                rowsendcnt[i*(dimens[0]) + j] = block_decompose(n, dimens[0], j);
                if(j == 0){
                    rowdisplays[i*(dimens[0]) + j] = (i == 0) ? 0 : rowdisplays[(i-1)*(dimens[0]) + j] + n*block_decompose(n, dimens[0], i-1);
                }else{
                    rowdisplays[i*(dimens[0]) + j] = rowdisplays[i*dimens[0] + j - 1] + rowsendcnt[i*dimens[0] + j - 1];
                }
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
        if(cordas[0] >= temp2){
            tag = 1;
        }else{
            tag = 0;
        }
        MPI_Comm_split(comm, tag, 1, &row_comm);

        if(tag == 0){
            MPI_Scatterv(&input_matrix[n*temp1], rowsendcnt, rowdisplays, MPI_DOUBLE, &receivebuffer[colreceivecnt*temp1], colreceivecnt, MPI_DOUBLE, localrank, row_comm);
        }
        MPI_Comm_free(&row_comm);
    }
    *local_matrix = receivebuffer;

}

void transpose_bcast_vector(const int n, double* col_vector, double* row_vector, MPI_Comm comm)
{
    int d[2];
    int loop[2];
    int mesh_location[2];
    MPI_Cart_get(comm, 2, d, loop, mesh_location);

    // get # of col and row
    int m = d[0];
    int rowcount = block_decompose(n, m, mesh_location[0]);
    int num_cols = block_decompose(n, m, mesh_location[1]); 
    
    if (mesh_location[0] == 0 && mesh_location[1] == 0)
    {
        for (int i = 0; i < rowcount; i++) {
            row_vector[i] = col_vector[i];
        }
    }else if (mesh_location[1] == 0){
    
        int locationDIAG[] = {mesh_location[0], mesh_location[0]};
        int rankDcopy; 
        int rankD; 
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
}

void distributed_matrix_vector_mult(const int n, double* local_A, double* local_x, double* local_y, MPI_Comm comm)
{
    MPI_Comm COM_ROW;
    int d[2];
    int loop[2];
    int mesh_location[2];
    MPI_Cart_get(comm, 2, d, loop, mesh_location);
    MPI_Comm_split(comm, mesh_location[0], mesh_location[1], &COM_ROW);


    int m = d[0]; 
    int rowcount = block_decompose(n, m, mesh_location[0]);


    int num_cols = block_decompose(n, m, mesh_location[1]);

    double *transposed_x = new double[num_cols]; 
    transpose_bcast_vector(n, local_x, transposed_x, comm);


    double *result = new double[rowcount];

    for (int i = 0; i < rowcount; i++)
    {
        result[i] = 0;
        for (int j = 0; j < num_cols; j++){
            result[i] =result[i] + local_A[i*num_cols + j] * transposed_x[j];
        }
    }

    MPI_Reduce(result, local_y, rowcount, MPI_DOUBLE, MPI_SUM, 0, COM_ROW);
}

// Solves Ax = b using the iterative jacobi method
void distributed_jacobi(const int n, double* local_A, double* local_b, double* local_x, MPI_Comm comm, int max_iter, double l2_termination)
{
    // TODO
    // retrieve Cartesian topology information
    // int dimension[2];
    // int periods[2];
    // int coordinates[2];
    // MPI_Cart_get(comm, 2, dimension, periods, coordinates);

    // // rank of coordinate_zero
    // int coordinate_zero[] = {0,0};
    // int rankzero; 
    // MPI_Cart_rank(comm, coordinate_zero, &rankzero);

    // // number of rows or columns
    // int m = dimension[0];
    // int num_rows = block_decompose(n, m, coordinates[0]); 
    // int num_cols = block_decompose(n, m, coordinates[1]); 

    // // initialize R, diagonal elements as 0, others as A(i,j)
    // double* R = new double[num_rows*num_cols];
    // for (int i = 0; i < num_rows; i++){
    //     for (int j = 0; j < num_cols; j++){
    //         if (coordinates[0]==coordinates[1] && i==j){
    //             R[i*num_cols + j] = 0;
    //         }
    //         else {
    //             R[i*num_cols + j] = local_A[i*num_cols + j];
    //         }
    //     }
    // }
        
    // // calculate matrix D = diag(A)
    // MPI_Comm row_comm;
    // MPI_Comm_split(comm, coordinates[0], coordinates[1], &row_comm);
    // MPI_Comm column_comm;
    // MPI_Comm_split(comm, coordinates[1], coordinates[0], &column_comm);
    // double *temp = new double[num_rows];
    // for (int i = 0; i < num_rows; i++) {
    //     if (coordinates[0]==coordinates[1]){
    //         temp[i] = local_A[i*num_cols + i];
    //     }
    //     else {
    //         temp[i] = 0.0;
    //     }
    // }
    // double *D_diag = NULL;
    // if (coordinates[1] == 0) {
    //     D_diag = new double[num_rows];
    // }
    // MPI_Reduce(temp, D_diag, num_rows, MPI_DOUBLE, MPI_SUM, 0, row_comm);

    // // initialize local_x
    // for (int i = 0; i < num_rows; i++){
    //     local_x[i] = 0.0;
    // }

    // // product of R*x at first column
    // double *sum_Rx = NULL; 
    // if (coordinates[1] == 0){
    //     sum_Rx = new double[num_rows];
    // }

    // // product of A*x at first column
    // double *sum_Ax = NULL;
    // if (coordinates[1] == 0){
    //     sum_Ax = new double[num_rows];
    // }

    // // iteration check, continue or stop
    // bool iteration_status = false;

    // // iterative calcualtion of x
    // for (int iter = 0; iter < max_iter; iter++)
    // {
    //     // calculate product of R*x at first column
    //     distributed_matrix_vector_mult(n, R, local_x, sum_Rx, comm);

    //     // calculate product of A*x at first column
    //     distributed_matrix_vector_mult(n, local_A, local_x, sum_Ax, comm);
        
    //     // compare the error with l2_termination
    //     if (coordinates[1] == 0)
    //     {
    //         double sum_error = 0;
    //         double local_error = 0;
    //         for (int i = 0; i < num_rows; i++) {
    //             local_error += (sum_Ax[i]-local_b[i])*(sum_Ax[i]-local_b[i]);
    //         }
    //         // reduce/sum local_error to processor 0
    //         MPI_Reduce(&local_error, &sum_error, 1, MPI_DOUBLE, MPI_SUM, 0, column_comm);

    //         if (coordinates[0] == 0 && sum_error < l2_termination)
    //         {
    //             iteration_status = true;
    //         }
    //     }

    //     // brocast result of iteration check to other processors
    //     MPI_Bcast(&iteration_status, 1, MPI::BOOL, rankzero, comm);

    //     // iteration check, stop or continue
    //     if (iteration_status)
    //     {
    //         break;
    //     }
    //     else if (coordinates[1] == 0) 
    //     {
    //         // update local_x with D^(-1)*(b-R*x)
    //         for (int i = 0; i < num_rows; i++) {
    //             local_x[i] = (local_b[i]-sum_Rx[i]) / D_diag[i];
    //         }
    //     }
    // }
    int cordas[2], dimens[2], timeslots[2];
    int localrank;
    int restdimens[2] = {0, 0};
    bool status = false;
    MPI_Cart_get(comm, 2, dimens, timeslots, cordas);
    MPI_Cart_rank(comm, restdimens, &localrank);

    MPI_Comm row_comm;
    MPI_Comm column_comm;

    int rowcnt = block_decompose(n, dimens[0], cordas[0]);
    int colcnt = block_decompose(n, dimens[0], cordas[1]);

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
    for(int i = 0; i < rowcnt; ++i){
        temp[i] = (cordas[0] != cordas[1]) ? 0.0 : local_A[i*colcnt + i];
    }

    if(cordas[1] == 0){
        Diag = new double[rowcnt];
    }

    MPI_Reduce(temp, Diag, rowcnt, MPI_DOUBLE, MPI_SUM, 0, row_comm);

    for(int i = 0; i < rowcnt; ++i){
        local_x[i] = 0.0;
    }

    if(cordas[1] == 0){
        Rsum = new double[rowcnt];
        Asum = new double[rowcnt];
    }

    // if(cordas[1] == 0){
    //     Asum = new double[rowcnt];
    // }

    for(int i = 0; i < max_iter; ++i){
        distributed_matrix_vector_mult(n, R, local_x, Rsum, comm);
        distributed_matrix_vector_mult(n, local_A, local_x, Asum, comm);

        if(cordas[1] == 0){
            double errsum = 0.0;
            double errlocal = 0.0;
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

