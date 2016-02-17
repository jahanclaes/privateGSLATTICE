// This is a class for cylindrical boundary cPEPS
// cPEPS of size Lx * Ly
// cPEPS is periodic along x-driection, open in y-direction

// Bond dimension of the tensor along x-direction stays unchanged
// Bond dimension of the tensor along y-direction changes due to MPO-MPO multiplication

#ifndef My_cPEPS_CLASS_H
#define My_cPEPS_CLASS_H

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <Eigen/Dense>

#include "utility.h"

typedef Eigen::MatrixXd Mxd;


// Storage for cPEPS's tensors
class TensorList
{
public:
	TensorList();
	TensorList(int& xs, int& ylen, int& phy, int& xd, int& yd);
	~TensorList();
	
	// set methods
	void setTensorList(int& xs, int& ylen, int& phy, int& xd, int& yd);
	void setRandom();
	void setZero();
	
	void boostBD(int dxBD, int dyBD);
	
	// Data
	Mxd*** T; // (1) y-direction site label; (2) true physical index; (3) x-direction index label (pseudo physical);
	int* Dim;
	int xSite;
	int yL;
	int pD;
	int xBD;
	int yBD;
	int numParams;
	bool initted;
};

// Storage for cPEPS's derivative tensors
class dTensorList
{
public:
	dTensorList();
	dTensorList(int& xs, int& ylen, int& xd, int& yd);
	~dTensorList();
	
	// set methods
	void setdTensorList(int& xs, int& ylen, int& xd, int& yd);
	void setZero();
	
	void boostBD(int dxBD, int dyBD);
	
	// Data
	Mxd** T; // (1) y-direction site label; (2) x-direction index label (pseudo physical);
	int* Dim;
	int xSite;
	int yL;
	int xBD;
	int yBD;
	int numParams;
	bool initted;
};

class cPEPS
{
public:
	// Con/De-structor
	cPEPS();
	cPEPS(int xlen, int ylen, int phy, int xd, int yd, int maxBD, double tl);
	~cPEPS();
	
	// set methods
	void setcPEPS(int xlen, int ylen, int phy, int xd, int yd, int maxBD, double tl);
	void setUniform();
	void setNearUniform();
	void setProductState();
	void setNearProductState();
	
	void boostBD(int dxBD, int dyBD);
	
	bool initted;
	
	// size
	int xL;
	int yL;
	int pD; // phyical dimension
	int numParams;
	
	// bond dimensions
	int xBD;
	int yBD;
	
	int max_yBD; // for compression
	double tol; // for compression
	
	// Tensor networks contained in cPEPS
	TensorList*   TN;   // ordered by column
	dTensorList* dTN;  // ordered by column
	
	// given a configuration, pick out the tensors and wrap them as MPOs
	MPO* tMPO;
	void buildMPO(int** phyC);
	
	// contraction method
	double contractPEPS(int** phyC);  // contract cPEPS given a physical configuration
	
	// derivative method 
	void diffPEPS(int** phyC);         // All derivative

	// print
	void printTN(int** phyC);
	void printDiffTN();
	
};


#endif