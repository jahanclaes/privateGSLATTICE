#ifndef My_Iterative_Compression_H
#define My_Iterative_Compression_H

#include <cstdlib>
#include <Eigen/Dense>
#include "mpo.h"

typedef Eigen::MatrixXd Mxd;

void buildR(MPO& H, MPO& Hc, Mxd* CR);

void updateSite(MPO& H, MPO& Hc, Mxd* CL, Mxd* CR, int& site, char& direc, double& dt);

void updateEnvr(MPO& H, MPO& Hc, Mxd* CL, Mxd* CR, int& site, char& direc, double& dt);

void iterCompress(bool is_Random, int maxIter, double tol, int nbD, MPO& H, bool print_info);

#endif