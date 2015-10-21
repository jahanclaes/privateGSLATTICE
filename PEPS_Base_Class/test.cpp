#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <stdint.h>
#include <ctime>
#include <Eigen/Dense>
#include <limits>
// using namespace Eigen;
using namespace std;
#include "lapack_wrapper_v2.h"
#include "PEPS_Base.cpp"


int ** PhyConfig;
int W = 4;
int L = 4;
int phyD = 4;
int D = 5;
int VD = 20;

int main (int argc, char const *argv[])
{
	cout.precision(std::numeric_limits<double>::digits10+2);
	////////////////////////////////////////////////
	PhyConfig = new int * [W];
	for(int i = 0; i < W; i++)
	{
		PhyConfig[i] = new int [L]();
	}
	////////////////////////////////////////////////
	// PEPS_Base peps(W, L, phyD, D, VD);
	PEPS_Base * peps;
	peps = new PEPS_Base(W, L, phyD, D, VD);
	cout<<peps->contract(PhyConfig, 2)<<endl<<endl;
	cout<<peps->contract(PhyConfig, 1)<<endl<<endl;
	cout<<peps->contract(PhyConfig, 0)<<endl<<endl;
	peps->diff(PhyConfig);
	////////////////////////////////////////////////
	delete peps;
	return 0;
}

