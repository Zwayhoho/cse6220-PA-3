/**
 * @file    jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements matrix vector multiplication and Jacobi's method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
#include "jacobi.h"

/*
 * TODO: Implement your solutions here
 */

// my implementation:
#include <iostream>
#include <math.h>
using namespace std;

// Calculates y = A*x for a square n-by-n matrix A, and n-dimensional vectors x
// and y
void matrix_vector_mult(const int n, const double* A, const double* x, double* y)
{
    // TODO
    for (int i = 0; i < n; ++i) {
        y[i] = 0;
        for (int j = 0; j < n; ++j) {
            y[i] += x[j] * A[i * n + j];
        }
    }
}

// Calculates y = A*x for a n-by-m matrix A, a m-dimensional vector x
// and a n-dimensional vector y
void matrix_vector_mult(const int n, const int m, const double* A, const double* x, double* y)
{
    // TODO
    for (int i = 0; i < n; ++i) {
        y[i] = 0;
        for (int j = 0; j < m; ++j) {
            y[i] += x[j] * A[i * m + j];
        }
    }
}

// implements the sequential jacobi method
void jacobi(const int n, double* A, double* b, double* x, int max_iter, double l2_termination)
{
    // TODO
    double Diag[n], R[n * n], tempx[n];
    double templ2 = 0.0;
    int cnt = 0, stoppoint;

    //Set the matrix
    for(int i = 0; i < n; ++i){
        x[i] = 0;
        for(int j =0; j < n; ++j){
            if(i != j){
                R[i * n + j] = A[i * n + j];
            }else{
                R[i * n + j] = 0;
                Diag[i] = A[i * n + j];
            }
        }
    }

    while(cnt++ < max_iter){
        matrix_vector_mult(n, R, x, tempx);
        for(int i = 0; i < n; ++i){
            x[i] = ((b[i]-tempx[i]) / Diag[i]);
        }
        matrix_vector_mult(n, A, x, tempx);
        for(int i = 0; i < n; ++i){
            templ2 += ((tempx[i] - b[i]) * (tempx[i] - b[i]));
        }
        stoppoint = sqrt(templ2);
        if(l2_termination > stoppoint){
            break;
        }
    }

}
