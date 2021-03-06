#ifndef My_MPO_CLASS_H
#define My_MPO_CLASS_H

#include <string>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "lapack_wrapper.h"

typedef Eigen::MatrixXd Mxd;

struct SpM
{
	std::vector<int> r;
	std::vector<int> c;
	std::vector<double> v;
};

class MPO{
public:
	int Len;
	int pD;
	int bD;
	Mxd ** M;
	SpM ** H;
	int *  Dim;
	
	double norm;
	
	bool if_init;
	bool if_shell;
	
	MPO ();
	MPO (int l, int pd, int bd);
	~ MPO ();

	void setMPO(int l, int pd, int bd);
	void setShellMPO(int l, int pd, int bd);
	void clearMPO();	
	void copyMPO(const MPO& other);
	void addMPO(double coeff, const MPO& other);
	
	void buildSpMPO();
	
	void setZero();
	void setZero(std::string s);
	void setZero(int nbD);
	void setRand();
	void setRand(int nbD);
	void setIdentity();
	void setIdentity(int nbD);
	void square();
	
	void LC();               // Left  Canonicalization
	void RC();               // right Canonicalization
	void moveRight(int site);
	void moveLeft(int site);
	void compressL(int nbD); // Compress from left
	void EE();
	
	double trace();
};

#endif