#include <iostream>
#include <cmath>

#include "../lapack_wrapper.h"
#include "../lapack_wrapper.cpp"
#include "../mpo.h"
#include "../mpo.cpp"
#include "../iterCompress.h"
#include "../iterCompress.cpp"
#include "../utility.h"
#include "../utility.cpp"
#include "../cpeps.h"
#include "../cpeps.cpp"

using namespace std;

int main (int argc, char const *argv[])
{
	int xL = 2;
	int yL = 2;
	int pD = 2;
	int xD = 4;
	int yD = 4;
	int maxD = 12;
	double tol = 1E-10;
	
	int** phyC = new int* [xL];
	for(int i = 0; i < xL; ++i)
	{
		phyC[i] = new int [yL]();
	}
	
	cPEPS cp(xL,yL,pD,xD,yD,maxD,tol);
	// for(int i = 1; i < 16; ++i)
	// {
	// 	cp.TN[0].T[0][0][i].setZero();
	// 	cp.TN[0].T[1][0][i].setZero();
	// 	cp.TN[1].T[0][0][i].setZero();
	// 	cp.TN[1].T[1][0][i].setZero();
	// }
	// cout<<cp.TN[0].T[0][0][0]<<endl;
	// cout<<cp.TN[0].T[1][0][0]<<endl;
	cout<<cp.TN[0].T[0][0][0] * cp.TN[0].T[1][0][0]<<endl;
	
	// cout<<cp.TN[1].T[0][0][0]<<endl;
	// cout<<cp.TN[1].T[1][0][0]<<endl;
	cout<<cp.TN[1].T[0][0][0] * cp.TN[1].T[1][0][0]<<endl;
	cout<<(cp.TN[0].T[0][0][0] * cp.TN[0].T[1][0][0]) * (cp.TN[1].T[0][0][0] * cp.TN[1].T[1][0][0])<<endl;
	// cp.setUniform();
	// cp.setNearUniform();
	// cp.setProductState();
	// cp.setNearProductState();
	cout<<"Number of parameters per config = "<<cp.numParams<<endl;
	cp.buildMPO(phyC);
	// cp.printTN(phyC);
	double ct = cp.contractPEPS(phyC);
	cout<<"Contraction of the PEPS = "<<ct<<endl;
	cp.diffPEPS(phyC);
	cout<<"diff done!"<<endl;
	cp.printDiffTN();
	return 0;
}