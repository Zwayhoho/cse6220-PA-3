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
#include <vector>

// Calculates y = A*x for a square n-by-n matrix A, and n-dimensional vectors x
// and y
void matrix_vector_mult(const int n, const double* A, const double* x, double* y)
{
    // TODO
	for (int i=0; i<n; i++){
		y[i]=0;
		for (int j=0; j<n; j++){
			y[i] = y[i] + A[i*n+j]*x[j];
		}
	}
}

// Calculates y = A*x for a n-by-m matrix A, a m-dimensional vector x
// and a n-dimensional vector y
void matrix_vector_mult(const int n, const int m, const double* A, const double* x, double* y)
{
    // TODO
	for (int i=0; i<n; i++){
		y[i]=0;
		for (int j=0; j<m; j++){
			y[i] = y[i] + A[i*n+j]*x[j];
		}
	}
}

// implements the sequential jacobi method
void jacobi(const int n, double* A, double* b, double* x, int max_iter, double l2_termination)
{
    // TODO
    	// Create array for Diag, R and y
	std::vector<double> Diag(n);
	std::vector<double> R(n*n);
	std::vector<double> y(n);

	// initialization for iteration, x, R and Diag
	int iteration = 0;
	for (int i=0; i<n; i++){
		x[i] = 0.0;
	}

	for (int i = 0; i < n; ++i) {
		Diag[i] = A[i*n + i];
		for (int j = 0; j < n; ++j) {
			if (i==j){
				R[i*n + j] = 0.0;
			} 
			else{
				R[i*n + j] = A[i*n + j];
			}
		}
	}

	//Do the interation
	while (iteration < max_iter) {
	
		//Calculate value of y
		matrix_vector_mult(n, &A[0], &x[0], &y[0]);

		//initialization of l2
		double l2 = 0.0;

		for (int i = 0; i < n; i++){
			l2 += (y[i] - b[i])*(y[i] - b[i]);
		}
		l2 = sqrt(l2);
		
		//termination criterion
		if (l2 > l2_termination) {
			matrix_vector_mult(n, &R[0], &x[0], &y[0]);
			for (int i = 0; i < n; i++){
				x[i] = (b[i] - y[i]) / Diag[i];
			}
		} 
		else{
				break;
			}
		iteration = iteration + 1;
	}
}
