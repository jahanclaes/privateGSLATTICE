// This is a class for cylindrical boundary cPEPS
// cPEPS of size Lx * Ly
// cPEPS is periodic along x-driection, open in y-direction

// Bond dimension of the tensor along x-direction stays unchanged
// Bond dimension of the tensor along y-direction changes due to MPO-MPO multiplication

#ifndef My_cPEPS_CLASS
#define My_cPEPS_CLASS

typedef Eigen::MatrixXd Mxd;

TensorList::TensorList()
{	
	initted = false;
}

TensorList::TensorList(int& xs, int& ylen, int& phy, int& xd, int& yd): xSite(xs), yL(ylen), pD(phy), xBD(xd), yBD(yd)
{
	Dim = new int [yL+1]();
	Dim[0] = 1;
	for(int i = 1; i < yL; ++i)
	{
		int pw = std::min(i,yL-i);
		Dim[i] = std::pow(xBD*xBD,pw);
		if(Dim[i]<0||log2(Dim[i])<pw)
		{
			Dim[i] = yBD;
		}else
		{
			Dim[i] = std::min(Dim[i],yBD);
		}
	}
	Dim[yL] = 1;
	
	T = new Mxd** [yL];
	for(int i = 0; i < yL; ++i)
	{
		T[i] = new Mxd* [pD];
		for(int j = 0; j < pD; ++j)
		{
			T[i][j] = new Mxd [xBD*xBD];
		}
	}
	
	// compute number of parameters in the TensorList
	numParams = 0;
	for(int i = 0; i < yL; ++i)
	{
		numParams += Dim[i] * Dim[i+1];
	}
	numParams *= (pD*xBD*xBD);
	
	initted = true;
	
	setRandom();
}


TensorList::~TensorList()
{
	// std::cout<<"TensorList destructor"<<std::endl;
	if(initted)
	{
		for(int i = 0; i < yL; ++i)
		{
			for(int j = 0; j < pD; ++j)
			{
				delete [] T[i][j];
			}
			delete [] T[i];
		}
		delete [] T;
		delete [] Dim;
		initted = false;
	}

}

void TensorList::setTensorList(int& xs, int& ylen, int& phy, int& xd, int& yd)
{
	if(!initted)
	{
		xSite = xs;
		yL    = ylen;
		pD    = phy;
		xBD   = xd;
		yBD   = yd;
		
		Dim = new int [yL+1]();
		Dim[0] = 1;
		for(int i = 1; i < yL; ++i)
		{
			int pw = std::min(i,yL-i);
			Dim[i] = std::pow(xBD*xBD,pw);
			if(Dim[i]<0||log2(Dim[i])<pw)
			{
				Dim[i] = yBD;
			}else
			{
				Dim[i] = std::min(Dim[i],yBD);
			}
		}
		Dim[yL] = 1;
	
		T = new Mxd** [yL];
		for(int i = 0; i < yL; ++i)
		{
			T[i] = new Mxd* [pD];
			for(int j = 0; j < pD; ++j)
			{
				T[i][j] = new Mxd [xBD*xBD];
			}
		}
	
		// compute number of parameters in the TensorList
		numParams = 0;
		for(int i = 0; i < yL; ++i)
		{
			numParams += Dim[i] * Dim[i+1];
		}
		numParams *= (pD*xBD*xBD);
	
		initted = true;
	}
}

void TensorList::setRandom()
{
	if(!initted)
	{
		std::cout<<"Tensors' dimensions are not set yet!"<<std::endl;
		abort();
	}
	
	for(int i = 0; i < yL; ++i)
	{
		for(int j = 0; j < pD; ++j)
		{
			for(int k = 0; k < xBD*xBD; ++k)
			{
				T[i][j][k].setRandom(Dim[i],Dim[i+1]);
			}
		}
	}
}

void TensorList::setZero()
{
	if(!initted)
	{
		std::cout<<"Tensors' dimensions are not set yet!"<<std::endl;
		abort();
	}
	
	for(int i = 0; i < yL; ++i)
	{
		for(int j = 0; j < pD; ++j)
		{
			for(int k = 0; k < xBD*xBD; ++k)
			{
				T[i][j][k].setZero(Dim[i],Dim[i+1]);
			}
		}
	}
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////


dTensorList::dTensorList()
{	
	initted = false;
}

dTensorList::dTensorList(int& xs, int& ylen, int& xd, int& yd): xSite(xs), yL(ylen), xBD(xd), yBD(yd)
{
	Dim = new int [yL+1]();
	Dim[0] = 1;
	for(int i = 1; i < yL; ++i)
	{
		int pw = std::min(i,yL-i);
		Dim[i] = std::pow(xBD*xBD,pw);
		if(Dim[i]<0||log2(Dim[i])<pw)
		{
			Dim[i] = yBD;
		}else
		{
			Dim[i] = std::min(Dim[i],yBD);
		}
	}
	Dim[yL] = 1;
	
	T = new Mxd* [yL];
	for(int i = 0; i < yL; ++i)
	{
		T[i] = new Mxd [xBD*xBD];
	}
	
	// compute number of parameters in the TensorList
	numParams = 0;
	for(int i = 0; i < yL; ++i)
	{
		numParams += Dim[i] * Dim[i+1];
	}
	numParams *= (xBD*xBD);
	
	initted = true;
	
	setZero();
}


dTensorList::~dTensorList()
{
	// std::cout<<"dTensorList destructor"<<std::endl;
	if(initted)
	{
		for(int i = 0; i < yL; ++i)
		{
			delete [] T[i];
		}
		delete [] T;
		delete [] Dim;
		initted = false;
	}

}

void dTensorList::setdTensorList(int& xs, int& ylen, int& xd, int& yd)
{
	if(!initted)
	{
		xSite = xs;
		yL    = ylen;
		xBD   = xd;
		yBD   = yd;
		
		Dim = new int [yL+1]();
		Dim[0] = 1;
		for(int i = 1; i < yL; ++i)
		{
			int pw = std::min(i,yL-i);
			Dim[i] = std::pow(xBD*xBD,pw);
			if(Dim[i]<0||log2(Dim[i])<pw)
			{
				Dim[i] = yBD;
			}else
			{
				Dim[i] = std::min(Dim[i],yBD);
			}
		}
		Dim[yL] = 1;
	
		T = new Mxd* [yL];
		for(int i = 0; i < yL; ++i)
		{
			T[i] = new Mxd [xBD*xBD];
		}
	
		// compute number of parameters in the TensorList
		numParams = 0;
		for(int i = 0; i < yL; ++i)
		{
			numParams += Dim[i] * Dim[i+1];
		}
		numParams *= (xBD*xBD);
	
		initted = true;
	}
}

void dTensorList::setZero()
{
	if(!initted)
	{
		std::cout<<"Tensors' dimensions are not set yet!"<<std::endl;
		abort();
	}
	
	for(int i = 0; i < yL; ++i)
	{
		for(int k = 0; k < xBD*xBD; ++k)
		{
			T[i][k].setZero(Dim[i],Dim[i+1]);
		}
	}
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////


cPEPS::cPEPS()
{
	initted = false;
}

cPEPS::cPEPS(int xlen, int ylen, int phy, int xd, int yd, int maxBD, double tl): xL(xlen), yL(ylen), pD(phy), xBD(xd), yBD(yd), max_yBD(maxBD), tol(tl)
{
	TN  = new TensorList [xL];
	dTN = new dTensorList [xL];
	
	numParams = 0;
	for(int i = 0; i < xL; ++i)
	{
		TN[i].setTensorList(i,yL,pD,xBD,yBD);
		TN[i].setRandom();
		numParams += TN[i].numParams;
		
		dTN[i].setdTensorList(i,yL,xBD,yBD);
		dTN[i].setZero();
	}
	
	tMPO = new MPO [xL];
	for(int i = 0; i < xL; ++i)
	{
		tMPO[i].setMPO(yL,xBD,yBD,0);
		tMPO[i].setZero(yBD);
	}
	
	initted = true;
}

cPEPS::~cPEPS()
{
	// std::cout<<"cPEPS destructor"<<std::endl;
	if(initted)
	{
		delete [] TN;
		delete [] dTN;
		delete [] tMPO;
	}
}

void cPEPS::setcPEPS(int xlen, int ylen, int phy, int xd, int yd, int maxBD, double tl)
{
	if(!initted)
	{
		xL      = xlen;
		yL      = ylen;
		pD      = phy;
		xBD     = xd;
		yBD     = yd;
		max_yBD = maxBD;
		tol     = tl;
		
		TN  = new TensorList [xL];
		dTN = new dTensorList [xL];
		
		numParams = 0;
		for(int i = 0; i < xL; ++i)
		{
			TN[i].setTensorList(i,yL,pD,xBD,yBD);
			TN[i].setRandom();
			numParams += TN[i].numParams;
		
			dTN[i].setdTensorList(i,yL,xBD,yBD);
			dTN[i].setZero();
		}
		
		tMPO = new MPO [xL];
		for(int i = 0; i < xL; ++i)
		{
			tMPO[i].setMPO(yL,xBD,yBD,0);
			tMPO[i].setZero(yBD);
		}
	
		initted = true;
	}
}

void cPEPS::buildMPO(int** phyC)
{
	if(!initted)
	{
		std::cout<<"cPEPS is not initialized yet!"<<std::endl;
		abort();
	}
	
	for(int x = 0; x < xL; ++x)
	{
		for(int y = 0; y < yL; ++y)
		{
			for(int i = 0; i < xBD*xBD; ++i)
			{
				tMPO[x].M[y][i] = TN[x].T[y][phyC[x][y]][i];
			}
		}
		tMPO[x].RC();
		// std::cout<<"Norm of MPO "<<x<<" = "<<tMPO[x].norm<<std::endl;
	}
	// for(int i = 0; i < tMPO[0].Len; ++i)
	// {
	// 	for(int j = 0; j < tMPO[0].pD*tMPO[0].pD; ++j)
	// 	{
	// 		std::cout<<tMPO[0].M[i][j]<<std::endl;
	// 	}
	// }
}

double cPEPS::contractPEPS(int** phyC)
{
	buildMPO(phyC);
	MPO H;
	H.copyMPO(tMPO[xL-1]);
	for(int x = xL-2; x >=0 ; --x)
	{
		applyMPO(tMPO[x],H,0,'T','I');
		// std::cout<<"BD of MPO at step "<<x<<" = "<<H.bD<<std::endl;
		if(H.bD>max_yBD) iterCompress(true, 6, tol, max_yBD, H, true);
	}
	// for(int i = 0; i < H.Len; ++i)
	// {
	// 	for(int j = 0; j < H.pD*H.pD; ++j)
	// 	{
	// 		std::cout<<tMPO[0].M[i][j].norm()<<std::endl;
	// 	}
	// }
	// std::cout<<H.norm<<std::endl;
	return H.trace();
}

double cPEPS::diffPEPS(int** phyC)
{
	buildMPO(phyC);
	MPO* LEnv = new MPO [xL];
	MPO* REnv = new MPO [xL];
	
	// Build envioroment MPOs
	MPO HR;
	HR.copyMPO(tMPO[xL-1]);
	REnv[xL-1].copyMPO(HR);
	for(int x = xL-2; x >=0 ; --x)
	{
		applyMPO(tMPO[x],HR,0,'T','I');
		// std::cout<<"BD of MPO at step "<<x<<" = "<<HR.bD<<std::endl;
		if(HR.bD>max_yBD) iterCompress(true, 6, tol, max_yBD, HR, true);
		REnv[x].copyMPO(HR);
	}
	// Build envioroment MPOs
	MPO HL;
	HL.copyMPO(tMPO[0]);
	LEnv[0].copyMPO(HL);
	for(int x = 1; x < xL ; ++x)
	{
		applyMPO(tMPO[x],HL,0,'B','I');
		// std::cout<<"BD of MPO at step "<<x<<" = "<<HL.bD<<std::endl;
		if(HL.bD>max_yBD) iterCompress(true, 6, tol, max_yBD, HL, true);
		LEnv[x].copyMPO(HL);
	}

	// Build the derivative tensors
	for(int x = 0; x < xL; ++x)
	{
		// Build the surrounding MPO
		MPO H;
		if(x==0)
			H.copyMPO(REnv[1]);
		else if(x==xL-1)
			H.copyMPO(LEnv[xL-2]);
		else
		{
			H.copyMPO(REnv[x+1]);
			applyMPO(LEnv[x-1],H,0,'B','I');
			// std::cout<<"BD of MPO at step "<<x<<" = "<<H.bD<<std::endl;
			if(H.bD>max_yBD) iterCompress(true, 6, tol, max_yBD, H, true);
		}
		
		// Build top and botttom environments
		Mxd* TEnv = new Mxd [yL];
		Mxd* BEnv = new Mxd [yL];
		for(int i = 0; i < yL-1; ++i)
		{
			Mxd tp;
			BEnv[i].setZero(tMPO[x].M[i][0].cols(), H.M[i][0].cols());
			TEnv[yL-1-i].setZero(tMPO[x].M[yL-1-i][0].rows(), H.M[yL-1-i][0].rows());
			if(i==0)
			{
				// Bottom
				for(int j = 0; j < xBD; ++j)
				{
					for(int k = 0; k < xBD; ++k)
					{
						BEnv[i].noalias() += tMPO[x].M[i][j*xBD+k].transpose() * H.M[i][k*xBD+j];
					}
				}
				
				// Top
				for(int j = 0; j < xBD; ++j)
				{
					for(int k = 0; k < xBD; ++k)
					{
						TEnv[yL-1-i].noalias() += tMPO[x].M[yL-1-i][j*xBD+k] * H.M[yL-1-i][k*xBD+j].transpose();
					}
				}
			}else
			{
				// Bottom
				for(int j = 0; j < xBD; ++j)
				{
					for(int k = 0; k < xBD; ++k)
					{
						BEnv[i] += tMPO[x].M[i][j*xBD+k].transpose() * BEnv[i-1] * H.M[i][k*xBD+j];
					}
				}
				
				// Top
				for(int j = 0; j < xBD; ++j)
				{
					for(int k = 0; k < xBD; ++k)
					{
						TEnv[yL-1-i] += tMPO[x].M[yL-1-i][j*xBD+k] * TEnv[yL-i] * H.M[yL-1-i][k*xBD+j].transpose();
					}
				}
			}
		}
		
		for(int y = 0; y < yL; ++y)
		{
			if(y==0)
			{
				for(int j = 0; j < xBD; ++j)
				{
					for(int k = 0; k < xBD; ++k)
					{
						dTN[x].T[y][k*xBD+j] = H.M[y][j*xBD+k] * TEnv[y+1].transpose();
					}
				}
			}else if(y==yL-1)
			{
				for(int j = 0; j < xBD; ++j)
				{
					for(int k = 0; k < xBD; ++k)
					{
						dTN[x].T[y][k*xBD+j] = BEnv[y-1] * H.M[y][j*xBD+k];
					}
				}
			}else
			{
				for(int j = 0; j < xBD; ++j)
				{
					for(int k = 0; k < xBD; ++k)
					{
						dTN[x].T[y][k*xBD+j] = BEnv[y-1] * H.M[y][j*xBD+k] * TEnv[y+1].transpose();
					}
				}
			}
		}
		
		delete [] TEnv;
		delete [] BEnv;
	}
	
	delete [] LEnv;
	delete [] REnv;
}

void cPEPS::setUniform()
{
	if(!initted)
	{
		std::cout<<"cPEPS is not initialized yet!"<<std::endl;
		abort();
	}
	
	for(int x = 0; x < xL; ++x)
	{
		for(int y = 0; y < yL; ++y)
		{
			for(int i = 0; i < pD; ++i)
			{
				for(int j = 0; j < xBD; ++j)
				{
					for(int k = 0; k < xBD; ++k)
					{
						TN[x].T[y][i][j*xBD+k].setIdentity();
					}
				}
			}
		}
	}
}

void cPEPS::setNearUniform()
{
	if(!initted)
	{
		std::cout<<"cPEPS is not initialized yet!"<<std::endl;
		abort();
	}
	
	for(int x = 0; x < xL; ++x)
	{
		TN[x].setRandom();
	}
	
	Mxd tp;
	
	for(int x = 0; x < xL; ++x)
	{
		for(int y = 0; y < yL; ++y)
		{
			for(int i = 0; i < pD; ++i)
			{
				for(int j = 0; j < xBD; ++j)
				{
					for(int k = 0; k < xBD; ++k)
					{
						tp = 0.1 * TN[x].T[y][i][j*xBD+k];
						TN[x].T[y][i][j*xBD+k].setIdentity();
						TN[x].T[y][i][j*xBD+k] += tp;
					}
				}
			}
		}
	}
}

void cPEPS::setProductState()
{
	if(!initted)
	{
		std::cout<<"cPEPS is not initialized yet!"<<std::endl;
		abort();
	}
	
	for(int x = 0; x < xL; ++x)
	{
		TN[x].setZero();
	}
	
	for(int x = 0; x < xL; ++x)
	{
		for(int y = 0; y < yL; ++y)
		{
			for(int i = 0; i < pD; ++i)
			{
				for(int j = 0; j < xBD; ++j)
				{
					for(int k = 0; k < xBD; ++k)
					{
						if( ((x+y)%2==0 && i==1 && j==k) || ((x+y)%2==1 && i==2  && j==k) )
						{
							TN[x].T[y][i][j*xBD+k].setIdentity();
						}else
						{
							TN[x].T[y][i][j*xBD+k].setZero();
						}
					}
				}
			}
		}
	}
}

void cPEPS::setNearProductState()
{
	if(!initted)
	{
		std::cout<<"cPEPS is not initialized yet!"<<std::endl;
		abort();
	}
	
	for(int x = 0; x < xL; ++x)
	{
		TN[x].setRandom();
	}
	
	Mxd tp;
	for(int x = 0; x < xL; ++x)
	{
		for(int y = 0; y < yL; ++y)
		{
			for(int i = 0; i < pD; ++i)
			{
				for(int j = 0; j < xBD; ++j)
				{
					for(int k = 0; k < xBD; ++k)
					{
						if( ((x+y)%2==0 && i==1  && j==k ) || ((x+y)%2==1 && i==2  && j==k ) )
						{
							tp = 0.1 * TN[x].T[y][i][j*xBD+k];
							TN[x].T[y][i][j*xBD+k].setIdentity();
							TN[x].T[y][i][j*xBD+k] += tp;
						}else
						{
							TN[x].T[y][i][j*xBD+k] *= 0.1;
						}
					}
				}
			}
		}
	}
}

void cPEPS::printTN(int** phyC)
{
	for(int x = 0; x < xL; ++x)
	{
		for(int y = 0; y < yL; ++y)
		{
			for(int i = 0; i < xBD; ++i)
			{
				for(int j = 0; j < xBD; ++j)
				{
					std::cout<<TN[x].T[y][phyC[x][y]][i*xBD+j]<<std::endl;
				}
			}
		}
		std::cout<<std::endl;
	}
}

void cPEPS::printDiffTN()
{
	for(int x = 0; x < xL; ++x)
	{
		for(int y = 0; y < yL; ++y)
		{
			for(int i = 0; i < xBD; ++i)
			{
				for(int j = 0; j < xBD; ++j)
				{
					std::cout<<dTN[x].T[y][i*xBD+j]<<std::endl;
				}
			}
		}
		std::cout<<std::endl;
	}
}





#endif