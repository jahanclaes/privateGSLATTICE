//////////////////////////////////////////////////////
// Version 1.0 -- Xiongjie Yu, 07/13/2014
// This is a wrapper file for useful Lapack functions
// Using Eigen3's matrix classes
// Containes QR factorization, SVD, density matrix decomposition
//////////////////////////////////////////////////////
// IMPORTANT NOTICE:
// The SVD and density matrix decomposition functions only support compression starting from the left, and only return thin U and (U^dagger * M) as V^dagger
//////////////////////////////////////////////////////
// debug mode (controls whether to output svd errors or not)
#define DEBUG_ENABLED_LapackWrapper false
#define Mxd Eigen::MatrixXd
//////////////////////////////////////////////////////
#include <algorithm>
#include <math.h>
//////////////////////////////////////////////////////
// Lapack functions
// Part A: real matrices
extern "C" void dgeqrf_(int*, int*, double*, int*, double*, double*, int*, int*);
extern "C" void dorgqr_(int*, int*, int*, double*, int*, double*, double*, int*, int*);
extern "C" void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
extern "C" void dsyev_(char*,char*,int*,double*,int*,double*,double*,int*,int*);
//extern "C" void dgemm_(char*,char*,int*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
//////////////////////////////////////////////////////
// A. Real matrices
// A.0 Matrices multiplication
//void dgemm(Mxd& M1, Mxd& M2, Mxd& R)
//{
//	char TRANSA = 'N';
//	char TRANSB = 'N';
//	int M = M1.rows();
//	int N = M2.cols();
//	if(M1.cols()!=M2.rows())
//	{
//		std::cout<<"Please make sure the matrix dimensions match!\n";
//		exit(1);
//	}
//	int K = M1.cols();
//	double ALPHA = 1.0;
//
//	int LDA = M;
//
//	int LDB = K;
//	double BETA = 0.0;
//	
//	R.setZero(M,N);
//	
//	int LDC = M;
//	///////////////////////////////
//	// Call Lapack
//	dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,M1.data(),&LDA,M2.data(),&LDB,&BETA,R.data(),&LDC);
//	///////////////////////////////
//}


// A.1 QR factorization
void rQR(const Mxd *MT, Mxd *Q){
	///////////////////////////////////
	// INITIALIZATION
	int M = MT[0].rows();
	int N = MT[0].cols();
	
	Mxd tQ = MT[0];

	int LDA, K;
	LDA = max(1,M);
	K = min(M,N);

	double * TAU = new double [K];
	int lwork = max(M,N)*max(M,N);//The dimension of the array WORK.
	double * work = new double [lwork];
	int INFO;
	///////////////////////////////////
	dgeqrf_(&M, &N, tQ.data(), &LDA, TAU, work, &lwork, &INFO);
	if (INFO!=0) {
		std::cout<<"Illegal value at "<<INFO<<std::endl;
	}
	if (N>M) {
		N = M;
		dorgqr_(&M, &N, &K, tQ.data(), &LDA, TAU, work, &lwork, &INFO);
		if (INFO!=0) {
			std::cout<<"QR Failed!\n";
		}
		Q[0]=tQ.block(0,0,M,N);
	} else {
		dorgqr_(&M, &N, &K, tQ.data(), &LDA, TAU, work, &lwork, &INFO);
		if (INFO!=0) {
			std::cout<<"QR Failed!\n";
		}
		Q[0]=tQ;
	}
	///////////////////////////////////
	//Free Space
	delete [] TAU;
	delete [] work;
}


// A.2 Density matrix decomposition
void dsyev(Mxd *M, double * evals)
{
	int N = M[0].cols();
    char jobz = 'V';
    char uplo = 'U';
    int lwork = max(1,3*N-1);//max(1, 1+6*N+2*N*N);
    double *work = new double [lwork];
    int info;
	dsyev_(&jobz,&uplo,&N,M[0].data(),&N,evals,work,&lwork,&info);
	if(info!=0)
	{
		std::cout<<"dsyev error info is "<<info<<std::endl;
	}
	delete [] work;
}

int rDenMatDecomp(const Mxd *M, int ds, double * sv, Mxd *UM, Mxd *VTM, int truncD, double ElasticLimit)
{
	int lowerTD = truncD;
	Mxd tp;
	if(ds<min(M[0].rows(),M[0].cols()))
	{
		std::cout<<"Please assign more space to the singular value container!\n";
		exit(1);
	}
	for(int i = 0; i < ds; i++)
	{
		sv[i]=0.0;
	}
	double svderror = 0.0;
	////////////////////////////////////////
	if(M[0].cols()>=M[0].rows()) // 1. row <= col
	{
		UM[0].noalias()=-1.0*M[0]*M[0].transpose();
		dsyev(UM,sv);
		if(truncD>=M[0].rows()) // presumably do not need to do truncation
		{
			for(int i = 0; i < M[0].rows(); i++)
			{
				if(sv[i]<0)
				{
					sv[i]=sqrt(-sv[i]);
				}else
				{
					sv[i]=0.0;
				}
			}
			VTM[0].noalias() = UM[0].transpose()*M[0];
			lowerTD = M[0].rows();
		}else if(truncD>0)  // do truncation
		{
			for(int i = 0; i < M[0].rows(); i++)
			{
				if(sv[i]<0)
				{
					sv[i]=sqrt(-sv[i]);
				}else
				{
					sv[i]=0.0;
				}
			}
			// int lowerTD = truncD;
			if(ElasticLimit>sv[truncD]) // maybe can do a more ambitious truncation
			{
				for(int i = truncD; i >= 0; i--)
				{
					if(ElasticLimit<=sv[i])
					{
						lowerTD = i+1;
						i = -1; //break
					}
				}
			}
			for(int i = lowerTD; i < M[0].rows(); i++)
			{
				svderror += sv[i];
			}
			tp=UM[0].block(0,0,M[0].rows(),truncD);
			UM[0]=tp;
			if(lowerTD<truncD) UM[0].block(0,lowerTD,M[0].rows(),truncD-lowerTD).setZero();
			VTM[0].noalias() = UM[0].transpose() * M[0];
		}
	}else // 2. row > col (this is the actual but more troublesome case for svd starting from the right hand side)
	{
		VTM[0].noalias()=-1.0*M[0].transpose()*M[0];
		dsyev(VTM,sv);
		if(truncD>=M[0].cols()) // do nothing
		{
			for(int i = 0; i < M[0].cols(); i++)
			{
				if(sv[i]<0)
				{
					sv[i]=sqrt(-sv[i]);
				}else
				{
					sv[i]=0.0;
				}
			}
			UM[0].noalias() = M[0]*VTM[0];
			for(int i = 0; i < M[0].cols(); i++)
			{
				UM[0].block(0,i,M[0].rows(),1) /= sv[i];
				VTM[0].block(i,0,1,M[0].cols()) *= sv[i];
			}
			lowerTD = M[0].cols();
		}else if(truncD>0) // do truncation
		{
			for(int i = 0; i < M[0].cols(); i++)
			{
				if(sv[i]<0)
				{
					sv[i]=sqrt(-sv[i]);
				}else
				{
					sv[i]=0.0;
				}
			}
			// int lowerTD = truncD;
			if(ElasticLimit>sv[truncD])
			{
				for(int i = truncD; i >= 0; i--)
				{
					if(ElasticLimit<=sv[i])
					{
						lowerTD = i+1;
						i = -1; //break
					}
				}
			}
			for(int i = lowerTD; i < M[0].cols(); i++)
			{
				svderror += sv[i];
			}
			UM[0].setZero(M[0].rows(),truncD);
			UM[0].block(0,0,M[0].rows(),lowerTD).noalias() = M[0]*VTM[0].block(0,0,VTM[0].rows(),lowerTD);
			tp=VTM[0].block(0,0,truncD,M[0].cols());
			VTM[0]=tp;
			if(truncD>lowerTD) VTM[0].block(lowerTD,0,truncD-lowerTD,M[0].cols()).setZero();
			for(int i = 0; i < lowerTD; i++)
			{
				UM[0].block(0,i,M[0].rows(),1) /= sv[i];
				VTM[0].block(i,0,1,M[0].cols()) *= sv[i];
			}
		}
	}
	if(DEBUG_ENABLED_LapackWrapper) std::cout<<"SVD error: "<<svderror<<"; actual truncation dimension: "<<lowerTD<<std::endl;
	return lowerTD;
}



// A.3 SVD
void rSVD(const Mxd &M, int ds, double * sv, Mxd &UM, Mxd &VTM, int truncD, double ElasticLimit){
	///////////////////////////////////
	// INITIALIZATION
	char jobU = 'S'; //all M columns of U are returned in array U
	char jobV = 'N'; //all N rows of V^T are returned in the array VT
	
	int rowN = M.rows();
	int colN = M.cols();
	int LDA = rowN;//The leading dimension of the array A
	
	int DofS;
	if (rowN >= colN) {
		DofS = colN;
	}else {
		DofS = rowN;
	}
	
	if (ds<DofS) {
		std::cout<<"Please assign more space to the singular value container!\n";
		exit(1);
	}
	for (int i = 0; i < ds; i++) {
		sv[i]=0.0;
	}
	
	UM.setZero(rowN,DofS);
	VTM.setZero(DofS,colN);
	
	int LDU = rowN;//The leading dimension of the array U
	
	int LDV = DofS;//The leading dimension of the array VT
	
	int lwork = max(3*min(rowN,colN)+max(rowN,colN),5*min(rowN,colN));
	
	double * work = new double [lwork];
	
	int INFO;
	
	Mxd tp = M;
	Mxd temp;
	///////////////////////////////////
	// CALL LAPACK
	dgesvd_(&jobU, &jobV, &rowN, &colN, tp.data(), &LDA, sv, UM.data(), &LDU, VTM.data(), &LDV, work, &lwork, &INFO);
	///////////////////////////////////
	if (INFO!=0) {
		std::cout<<"SVD Failed!\n";
	}
	///////////////////////////////////
	// Truncate
	if(truncD>=DofS) // do nothing
	{
		VTM.noalias() = UM.transpose() * M;
	}else if(truncD>0)
	{
		double svderror = 0.0;
		int lowerTD = truncD;
		if(ElasticLimit>sv[truncD])
		{
			for(int i = truncD; i >= 0; i--)
			{
				if(ElasticLimit<=sv[i])
				{
					lowerTD = i+1;
					i = -1; //break
				}
			}
		}
		for(int i = lowerTD; i < DofS; i++)
		{
			svderror += sv[i];
		}
		if(DEBUG_ENABLED_LapackWrapper) std::cout<<"SVD error: "<<svderror<<"; actual truncation dimension: "<<lowerTD<<std::endl;
		temp=UM.block(0,0,M.rows(),truncD);
		UM=temp;
		if(lowerTD<truncD) UM.block(0,lowerTD,M.rows(),truncD-lowerTD).setZero();
		VTM.noalias() = UM.transpose() * M;
	}
	///////////////////////////////////
	//Free Space
	delete [] work;
}
