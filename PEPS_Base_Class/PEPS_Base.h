#ifndef PEPS_BASE
#define PEPS_BASE

// #include <Eigen/Dense>
// using namespace Eigen;
//#include "lapack_wrapper_v2.h"
//#define Mxd Eigen::MatrixXd
using Eigen::MatrixXd;
typedef MatrixXd Mxd;
class PEPS_Base{
	public:
	
	int L;
	int W;
	int phyD;
	int D;
	int VD;
	
	bool initted_;
	double SVD_Tolerance;
	double Var_Tolerance;
	
	Mxd **** PEPS_;
	Mxd *** DiffT_;
	Mxd ** UM;// upper MPS
	Mxd ** LM;// lower MPS
	////////////////////////////////////////////////////////////
	//                     Work Space                         //
	////////////////////////////////////////////////////////////
	Mxd ** VM;
	Mxd ** TVM;
	////////////////////////////////////////////////////////////
	
	//public:
	
	PEPS_Base(int length, int width, int physicalDim, int BondDim, int MaxBondDim);
	~PEPS_Base();
	
	////////////////////////////////////////////////////////////
	//                     Main Routines                      //
	////////////////////////////////////////////////////////////
	double build_UpLow_MPS(int ** PhyC, int Urow, int Lrow); // returns norm
	double contract(int ** PhyC, int ptm);
	double diff(int ** PhyC);
	void setValue(){}; // need to implement data structure later
	////////////////////////////////////////////////////////////
	//                    Sub Routines                        //
	////////////////////////////////////////////////////////////
	void setRandom();
	void setProductState(); // staggered 
	void setNearProductState(); // staggered 
	void setNearUniform();
	void setUniform(); 
	void M_setZero(int phyDim, int Dim, Mxd ** MPS);
	void M_setRandom(int phyDim, int Dim, Mxd ** MPS);
	double psiphi (int pD, Mxd ** Psi, Mxd ** Phi);
	double RC (int Dim, int pD, Mxd ** MPS);
	void CompressL (int Dim, int DT, int pD, Mxd ** VMPS, Mxd ** MPS);
	double IterCompress (int IfRandom, int MaxIter, double tol, int VDim, int Dim, int pD, Mxd ** VMPS, Mxd ** MPS);
	////////////////////////////////////////////////////////////
};


#endif
