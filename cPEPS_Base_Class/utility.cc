#ifndef My_UTILITY_FUNCTIONS
#define My_UTILITY_FUNCTIONS

#include "utility.h"

typedef Eigen::MatrixXd Mxd;

// B = A1 cross A2
void KroneckerProd(Mxd& A1, Mxd& A2, Mxd& B)
{
	assert(A1.size()>0 && A2.size()>0);
	
	B.setZero(A1.rows()*A2.rows(), A1.cols()*A2.cols());
	
	int r1 = A1.rows();
	int c1 = A1.cols();

	int r2 = A2.rows();
	int c2 = A2.cols();
	
	for(int i = 0; i < r1; ++i)
	{
		for(int j = 0; j < c1; ++j)
		{
			B.block(i*r2, j*c2, r2, c2) = A1(i,j)*A2;
		}
	}
}

// Apply A to B. Results are saved to B
// A is thought of as the unitary gates
// A is applied through sites Site to Site+L(A)-1 of B
// from the top or the bottom direction
// direc =T (top), B (bottom), op = 'I' (kept the same), 'T' (tranposed) 
void applyMPO(MPO& A, MPO& B, int site, char direc, char op)
{
	assert(A.if_init && B.if_init);
	assert(site+A.Len-1<B.Len && A.pD==B.pD);
	assert(direc=='T'||direc=='B');
	assert(op=='I'||op=='T');
	
	int pD = A.pD;
	
	Mxd* tp = new Mxd [pD*pD];
	Mxd temp;
	
	for(int i = 0; i < A.Len; ++i)
	{
		for(int j = 0; j < pD*pD; ++j)
		{
			int topIdx = j/pD; int botIdx = j%pD;
			for(int k = 0; k < pD; ++k)
			{
				if(direc=='T' && op=='I')
					KroneckerProd(A.M[i][topIdx*pD+k], B.M[site+i][k*pD+botIdx], temp);
				else if(direc=='B' && op=='I')
					KroneckerProd(B.M[site+i][topIdx*pD+k], A.M[i][k*pD+botIdx], temp);
				else if(direc=='T' && op=='T')
					KroneckerProd(A.M[i][k*pD+topIdx], B.M[site+i][k*pD+botIdx], temp);
				else if(direc=='B' && op=='T')
					KroneckerProd(B.M[site+i][topIdx*pD+k], A.M[i][botIdx*pD+k], temp);
				
				if(k==0)
					tp[j] = temp;
				else
					tp[j] += temp;
			}
		}
		for(int j = 0; j < pD*pD; ++j)
		{
			B.M[site+i][j] = tp[j];
		}
	}
	
	for(int i = 0; i < A.Len+1; ++i)
	{
		B.Dim[site+i] *= A.Dim[i];
	}
	
	B.bD = *std::max_element(B.Dim,B.Dim+B.Len+1);

	delete [] tp;
}


#endif