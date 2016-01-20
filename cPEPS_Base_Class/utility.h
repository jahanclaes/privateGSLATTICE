#ifndef My_UTILITY_FUNCTIONS_H
#define My_UTILITY_FUNCTIONS_H

#include <cmath>
#include <cassert>
#include <algorithm>
#include <Eigen/Dense>

#include "mpo.h"
#include "iterCompress.h"

typedef Eigen::MatrixXd Mxd;


// B = A1 cross A2
void KroneckerProd(Mxd& A1, Mxd& A2, Mxd& B);

// Apply A to B. Results are saved to B
// A is thought of as the unitary gates
// A is applied through sites Site to Site+L(A)-1 of B
// from the top or the bottom direction
// direc =T (top), B (bottom), op = 'I' (kept the same), 'T' (tranposed) 
void applyMPO(MPO& A, MPO& B, int site, char direc, char op);

#endif