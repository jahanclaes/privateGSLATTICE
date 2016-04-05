#ifndef PEPS_BASE_C
#define PEPS_BASE_C

//#include "PEPS_Base.h"
// #include <Eigen/Dense>
#define Mxd Eigen::MatrixXd

PEPS_Base::PEPS_Base(int length, int width, int physicalDim, int BondDim, int MaxBondDim)
{
	initted_ = false;
	SVD_Tolerance = 1e-10;
	Var_Tolerance = 1e-10;
	cout<<"Initializing PEPS_Base object\n";
	if(initted_)
	{
		std::cout<<"Please call Clear() to clear data member before changing the system size"<<std::endl;
	}else
	{
		L = length;
		W = width;
		phyD = physicalDim;
		D = BondDim;
		VD = MaxBondDim;
		////////////////////////////////////////////////////////////
		PEPS_ = new Mxd *** [W];
		for(int i = 0; i < W; i++)
		{
			PEPS_[i] = new Mxd ** [L];
			for(int j = 0; j < L; j++)
			{
				PEPS_[i][j] = new Mxd * [phyD];
				for(int k = 0; k < phyD; k++)
				{
					if(i==0||i==W-1)
					{
						PEPS_[i][j][k] = new Mxd [D];
						for(int l = 0; l < D; l++)
						{
							if(j==0)
							{
								PEPS_[i][j][k][l].setRandom(1,D);
							}else if(j==L-1)
							{
								PEPS_[i][j][k][l].setRandom(D,1);
							}else
							{
								PEPS_[i][j][k][l].setRandom(D,D);
							}
						}
					}else
					{
						PEPS_[i][j][k] = new Mxd [D*D];
						for(int l = 0; l < D*D; l++)
						{
							if(j==0)
							{
								PEPS_[i][j][k][l].setRandom(1,D);
							}else if(j==L-1)
							{
								PEPS_[i][j][k][l].setRandom(D,1);
							}else
							{
								PEPS_[i][j][k][l].setRandom(D,D);
							}
						}
					}
				}
			}
		}
		////////////////////////////////////////////////////////////
		DiffT_ = new Mxd ** [W];
		for(int i = 0; i < W; i++)
		{
			DiffT_[i] = new Mxd * [L];
			for(int j = 0; j < L; j++)
			{
				if(i==0||i==W-1)
				{
					DiffT_[i][j] = new Mxd [D];
					for(int l = 0; l < D; l++)
					{
						if(j==0)
						{
							DiffT_[i][j][l].setZero(1,D);
						}else if(j==L-1)
						{
							DiffT_[i][j][l].setZero(D,1);
						}else
						{
							DiffT_[i][j][l].setZero(D,D);
						}
					}
				}else
				{
					DiffT_[i][j] = new Mxd [D*D];
					for(int l = 0; l < D*D; l++)
					{
						if(j==0)
						{
							DiffT_[i][j][l].setZero(1,D);
						}else if(j==L-1)
						{
							DiffT_[i][j][l].setZero(D,1);
						}else
						{
							DiffT_[i][j][l].setZero(D,D);
						}
					}
				}
			}
		}
		////////////////////////////////////////////////////////////
		UM = new Mxd * [L];
		for (int i = 0; i < L; i++) {
			UM[i] = new Mxd [D]; // Physical degrees of freedom: 00,10,01,11
		}

		LM = new Mxd * [L];
		for (int i = 0; i < L; i++) {
			LM[i] = new Mxd [D]; // Physical degrees of freedom: 00,10,01,11
		}
		////////////////////////////////////////////////
		VM = new Mxd * [L];
		for (int i = 0; i < L; i++) {
			VM[i] = new Mxd [D]; // Physical degrees of freedom: 00,10,01,11
		}

		TVM = new Mxd * [L];
		for (int i = 0; i < L; i++) {
			TVM[i] = new Mxd [D]; // Physical degrees of freedom: 00,10,01,11
		}
		////////////////////////////////////////////////
		initted_ = true;
	}
}

PEPS_Base::~PEPS_Base()
{
	if(initted_)
	{
		for(int i = 0; i < W; i++)
		{
			for(int j = 0; j < L; j++)
			{
				for(int k = 0; k < phyD; k++)
				{
					delete [] PEPS_[i][j][k];
				}
				delete [] PEPS_[i][j];
			}
			delete [] PEPS_[i];
		}
		delete [] PEPS_;
		for(int i = 0; i < L; i++)
		{
			delete [] UM[i];
			delete [] LM[i];
			delete [] VM[i];
			delete [] TVM[i];
		}
		delete [] UM;
		delete [] LM;
		delete [] VM;
		delete [] TVM;
		for(int i = 0; i < W; i++)
		{
			for(int j = 0; j < L; j++)
			{
				delete [] DiffT_[i][j];
			}
			delete [] DiffT_[i];
		}
		delete [] DiffT_;
	}
}

double PEPS_Base::build_UpLow_MPS(int ** PhyC, int Urow, int Lrow)
{
	// initialization
	double unorm, lnorm;
	unorm=lnorm=1.0;
	int TD;
	for(int i = 0; i < D; i++)
	{
		// top (upper)
		for(int j = 0; j < L; j++)
		{
			UM[j][i] = PEPS_[0][j][PhyC[0][j]][i];
		}
		// bottom (lower)
		for(int j = 0; j < L; j++)
		{
			LM[j][i] = PEPS_[W-1][j][PhyC[W-1][j]][i];
		}
	}
	unorm *= RC(UM[0][0].cols(), D, UM);
	lnorm *= RC(LM[0][0].cols(), D, LM);
	// contraction
	if(Urow<0||Lrow>W-1||Urow>=Lrow)
	{
		cout<<"Input error! Make sure the row parameters are reasonable."<<endl;
	}else
	{
		/////////////////////////////
		// from 0 to ptm
		for(int l = 0; l < Urow; l++)
		{
			int tempD = VD;
			TD = UM[0][0].cols()*D;
			for(int k = 0; k < D; k++)
			{
				VM[0][k].setZero(1,TD);
				VM[L-1][k].setZero(TD,1);
				for (int i = 1; i < L-1; i++) {
					VM[i][k].setZero(TD,TD);
				}

				TVM[0][k].setZero(1,TD);
				TVM[L-1][k].setZero(TD,1);
				for (int i = 1; i < L-1; i++) {
					TVM[i][k].setZero(TD,TD);
				}
			}
			for (int i = 0; i < TD; i++) {
				int p1,p2;
				p1=i%D;
				p2=i/D;
				for(int j = 0; j < D*D; j++)
				{
					int p3, p4;
					p3=j%D;
					p4=j/D;
					VM[0][p3](0,i) += PEPS_[l+1][0][PhyC[l+1][0]][j](0,p1) * UM[0][p4](0,p2);
					VM[L-1][p3](i,0) += PEPS_[l+1][L-1][PhyC[l+1][L-1]][j](p1,0) * UM[L-1][p4](p2,0);
				}
			}
			// 2nd to last-but-1 site //
			for (int i = 1; i < L-1; i++) {
				//cout<<i<<endl;
				for (int j = 0; j < TD*TD; j++) {
					int p1,p2,p3,p4,p5,p6;
					p1=j/TD;
					p3=p1%D;
					p4=p1/D;
					p2=j%TD;
					p5=p2%D;
					p6=p2/D;
					for(int k = 0; k < D*D; k++)
					{
						int p7, p8;
						p7=k%D;
						p8=k/D;
						VM[i][p7](p1,p2) += PEPS_[l+1][i][PhyC[l+1][i]][k](p3,p5) * UM[i][p8](p4,p6);
					}
				}
			}
			// Compress and right canonicalize the enlarged MPS
			unorm = unorm * RC(TD, D, VM);
			if(TD<=tempD)
			{
				for (int i = 0; i < D; i++) {
					UM[0][i]=VM[0][i];
					UM[L-1][i]=VM[L-1][i];
				}
				for (int i = 1; i < L-1; i++) {
					for (int j = 0; j < D; j++) {
						UM[i][j]=VM[i][j];
					}
				}
			}else
			{
				double MPSerror;
				bool Done = false;
				while(!Done)
				{
					if(TD<=tempD)
					{
						for (int i = 0; i < D; i++) {
							UM[0][i]=VM[0][i];
							UM[L-1][i]=VM[L-1][i];
						}
						for (int i = 1; i < L-1; i++) {
							for (int j = 0; j < D; j++) {
								UM[i][j]=VM[i][j];
							}
						}
						Done = true;
					}else
					{
						for (int i = 0; i < D; i++) { //Save a copy of the enlarged MPS
							TVM[0][i]=VM[0][i];
							TVM[L-1][i]=VM[L-1][i];
						}
						for (int i = 1; i < L-1; i++) {
							for (int j = 0; j < D; j++) {
								TVM[i][j]=VM[i][j];
							}
						}
						for (int i = 0; i < D; i++) { //Save a copy of the enlarged MPS
							UM[0][i].resize(1,tempD);
							UM[L-1][i].resize(tempD,1);
						}
						for (int i = 1; i < L-1; i++) {
							for (int j = 0; j < D; j++) {
								UM[i][j].resize(tempD,tempD);
							}
						}
						CompressL(TD, tempD, D, TVM, UM);
						RC(tempD, D, UM);
						MPSerror = IterCompress(0, 1, Var_Tolerance, TD, tempD, D, VM, UM);
						if(abs(MPSerror)<Var_Tolerance)
						{
							Done = true;
						}
						else
						{
							tempD ++;
						}
					}
				}
				unorm = unorm * RC(tempD, D, UM);
			}
		}
		/////////////////////////////
		// from ptm to W-1
		for(int l = W-2; l >= Lrow; l--)
		{
			int tempD = VD;
			TD = LM[0][0].cols()*D;
			for(int k = 0; k < D; k++)
			{
				VM[0][k].setZero(1,TD);
				VM[L-1][k].setZero(TD,1);
				for (int i = 1; i < L-1; i++) {
					VM[i][k].setZero(TD,TD);
				}

				TVM[0][k].setZero(1,TD);
				TVM[L-1][k].setZero(TD,1);
				for (int i = 1; i < L-1; i++) {
					TVM[i][k].setZero(TD,TD);
				}
			}
			for (int i = 0; i < TD; i++) {
				int p1,p2;
				p1=i%D;
				p2=i/D;
				for(int j = 0; j < D*D; j++)
				{
					int p3, p4;
					p3=j/D;
					p4=j%D;
					VM[0][p3](0,i) += PEPS_[l][0][PhyC[l][0]][j](0,p1) * LM[0][p4](0,p2);
					VM[L-1][p3](i,0) += PEPS_[l][L-1][PhyC[l][L-1]][j](p1,0) * LM[L-1][p4](p2,0);
				}
			}
			// 2nd to last-but-1 site //
			for (int i = 1; i < L-1; i++) {
				for (int j = 0; j < TD*TD; j++) {
					int p1,p2,p3,p4,p5,p6;
					p1=j/TD;
					p3=p1%D;
					p4=p1/D;
					p2=j%TD;
					p5=p2%D;
					p6=p2/D;
					for(int k = 0; k < D*D; k++)
					{
						int p7, p8;
						p7=k/D;
						p8=k%D;
						VM[i][p7](p1,p2) += PEPS_[l][i][PhyC[l][i]][k](p3,p5) * LM[i][p8](p4,p6);
					}
				}
			}
			// Compress and right canonicalize the enlarged MPS
			lnorm = lnorm * RC(TD, D, VM);
			if(TD<=tempD)
			{
				for (int i = 0; i < D; i++) {
					LM[0][i]=VM[0][i];
					LM[L-1][i]=VM[L-1][i];
				}
				for (int i = 1; i < L-1; i++) {
					for (int j = 0; j < D; j++) {
						LM[i][j]=VM[i][j];
					}
				}
			}else
			{
				double MPSerror;
				bool Done = false;
				while(!Done)
				{
					if(TD<=tempD)
					{
						for (int i = 0; i < D; i++) {
							LM[0][i]=VM[0][i];
							LM[L-1][i]=VM[L-1][i];
						}
						for (int i = 1; i < L-1; i++) {
							for (int j = 0; j < D; j++) {
								LM[i][j]=VM[i][j];
							}
						}
						Done = true;
					}else
					{
						for (int i = 0; i < D; i++) { //Save a copy of the enlarged MPS
							TVM[0][i]=VM[0][i];
							TVM[L-1][i]=VM[L-1][i];
						}
						for (int i = 1; i < L-1; i++) {
							for (int j = 0; j < D; j++) {
								TVM[i][j]=VM[i][j];
							}
						}
						for (int i = 0; i < D; i++) { //Save a copy of the enlarged MPS
							LM[0][i].resize(1,tempD);
							LM[L-1][i].resize(tempD,1);
						}
						for (int i = 1; i < L-1; i++) {
							for (int j = 0; j < D; j++) {
								LM[i][j].resize(tempD,tempD);
							}
						}
						CompressL(TD, tempD, D, TVM, LM);
						RC(tempD, D, LM);
						MPSerror = IterCompress(0, 1, Var_Tolerance, TD, tempD, D, VM, LM);
						if(abs(MPSerror)<Var_Tolerance)
						{
							Done = true;
						}
						else
						{
							tempD ++;
						}
					}
				}
				lnorm = lnorm * RC(tempD, D, LM);
			}
		}
		/////////////////////////////
	}
	return unorm*lnorm;
}

double PEPS_Base::contract(int ** PhyC, int ptm)
{
	double norm = build_UpLow_MPS(PhyC, ptm, ptm+1);
	return norm*psiphi(D, UM, LM);
}

double PEPS_Base::diff(int ** PhyC)
{
	int r1, c1, r2, c2;
	int r;
	double norm;
	for(int r = 0; r < W; r++)
	{
		////////////////////////////////////////////
		if(r==0)
		{
			norm = build_UpLow_MPS(PhyC, r, r+1);
			////////////////////////////////////////////////
			for(int i = 0; i < D; i++)
			{
				DiffT_[r][0][i].setZero();
				DiffT_[r][L-1][i].setZero();
				for(int j = 1; j < L-1; j++)
				{
					DiffT_[r][j][i].setZero();
				}
			}
			Mxd * CRM = new Mxd [L];
			Mxd * CLM = new Mxd [L];
			Mxd * TM1 = new Mxd [D];
			////////////////////////////////////////////////
			// Use UM[0] and LM[1] to build CRM
			// Last site //
			r1=UM[L-1][0].rows();
			c1=UM[L-1][0].cols();
			r2=LM[L-1][0].rows();
			c2=LM[L-1][0].cols();
			CRM[L-1].setZero(r1,r2);
			for(int i = 0; i < D; i++)
			{
				CRM[L-1] += UM[L-1][i] * LM[L-1][i].transpose();
			}
			// Last-but-1 to the first site //
			for (int i = L-2; i >= 1; i--) {
				r1=UM[i][0].rows();
				c1=UM[i][0].cols();
				r2=LM[i][0].rows();
				c2=LM[i][0].cols();
				CRM[i].setZero(r1,r2);
				for(int j = 0; j < D; j++)
				{
					TM1[j].setZero(c1,r2);
				}
				for(int tid = 0; tid < D; tid++)
				{
					TM1[tid].noalias() = CRM[i+1] * LM[i][tid].transpose();
				}
				for(int j = 0; j < D; j++)
				{
					CRM[i].noalias() += UM[i][j] * TM1[j];
				}
			}
			////////////////////////////////////////////////
			// Use UM[0] and LM[1] to build CLM
			// 1st site //
			r1=UM[0][0].rows();
			c1=UM[0][0].cols();
			r2=LM[0][0].rows();
			c2=LM[0][0].cols();
			CLM[0].setZero(c1,c2);
			for(int i = 0; i < D; i++)
			{
				CLM[0].noalias() += UM[0][i].transpose() * LM[0][i];
			}
			// 2nd to the last site //
			for (int i = 1; i < L-1; i++) {
				r1=UM[i][0].rows();
				c1=UM[i][0].cols();
				r2=LM[i][0].rows();
				c2=LM[i][0].cols();
				CLM[i].setZero(c1,c2);
				for(int j = 0; j < D; j++)
				{
					TM1[j].setZero(r1,c2);
				}
				for(int tid = 0; tid < D; tid++)
				{
					TM1[tid].noalias() = CLM[i-1] * LM[i][tid];
				}
				for(int j = 0; j < D; j++)
				{
					CLM[i].noalias() += UM[i][j].transpose() * TM1[j];
				}
			}
			////////////////////////////////////////////////
			// Build the derivative matrix
			for(int i = 0; i < L; i++)
			{
				if(i==0)
				{
					for(int j = 0; j < D; j++)
					{
						DiffT_[r][i][j].noalias() = LM[i][j] * CRM[i+1].transpose();
					}
				}else if(i==L-1)
				{
					for(int j = 0; j < D; j++)
					{
						DiffT_[r][i][j].noalias() = CLM[i-1] * LM[i][j];
					}
				}else
				{
					for(int j = 0; j < D; j++)
					{
						DiffT_[r][i][j].noalias() = CLM[i-1] * LM[i][j] * CRM[i+1].transpose();
					}
				}
			}
			////////////////////////////////////////////////
			delete [] CRM;
			delete [] CLM;
			delete [] TM1;
			////////////////////////////////////////////////
		}else if(r==W-1)
		{
			norm = build_UpLow_MPS(PhyC, r-1, r);
			////////////////////////////////////////////////
			for(int i = 0; i < D; i++)
			{
				DiffT_[r][0][i].setZero();
				DiffT_[r][L-1][i].setZero();
				for(int j = 1; j < L-1; j++)
				{
					DiffT_[r][j][i].setZero();
				}
			}
			Mxd * CRM = new Mxd [L];
			Mxd * CLM = new Mxd [L];
			Mxd * TM1 = new Mxd [D];
			// DimW is not referenced
			////////////////////////////////////////////////
			// Use UM[W-2] and LM[W-1] to build CRM
			// Last site //
			r1=UM[L-1][0].rows();
			c1=UM[L-1][0].cols();
			r2=LM[L-1][0].rows();
			c2=LM[L-1][0].cols();
			CRM[L-1].setZero(r1,r2);
			for(int i = 0; i < D; i++)
			{
				CRM[L-1].noalias() += UM[L-1][i] * LM[L-1][i].transpose();
			}
			// Last-but-1 to the first site //
			for (int i = L-2; i >= 1; i--) {
				r1=UM[i][0].rows();
				c1=UM[i][0].cols();
				r2=LM[i][0].rows();
				c2=LM[i][0].cols();
				CRM[i].setZero(r1,r2);
				for(int j = 0; j < D; j++)
				{
					TM1[j].setZero(c1,r2);
				}
				for(int tid = 0; tid < D; tid++)
				{
					TM1[tid].noalias() = CRM[i+1] * LM[i][tid].transpose();
				}
				for(int j = 0; j < D; j++)
				{
					CRM[i].noalias() += UM[i][j] * TM1[j];
				}
			}
			////////////////////////////////////////////////
			// Use UM[W-2] and LM[W-1] to build CLM
			// 1st site //
			r1=UM[0][0].rows();
			c1=UM[0][0].cols();
			r2=LM[0][0].rows();
			c2=LM[0][0].cols();
			CLM[0].setZero(c1,c2);
			for(int i = 0; i < D; i++)
			{
				CLM[0].noalias() += UM[0][i].transpose() * LM[0][i];
			}
			// 2nd to the last site //
			for (int i = 1; i < L-1; i++) {
				r1=UM[i][0].rows();
				c1=UM[i][0].cols();
				r2=LM[i][0].rows();
				c2=LM[i][0].cols();
				CLM[i].setZero(c1,c2);
				for(int j = 0; j < D; j++)
				{
					TM1[j].setZero(r1,c2);
				}
				for(int tid = 0; tid < D; tid++)
				{
					TM1[tid].noalias() = CLM[i-1] * LM[i][tid];
				}
				for(int j = 0; j < D; j++)
				{
					CLM[i].noalias() += UM[i][j].transpose() * TM1[j];
				}
			}
			////////////////////////////////////////////////
			// Build the derivative matrix
			for(int i = 0; i < L; i++)
			{
				if(i==0)
				{
					for(int j = 0; j < D; j++)
					{
						DiffT_[r][i][j].noalias() = UM[i][j] * CRM[i+1];
					}
				}else if(i==L-1)
				{
					for(int j = 0; j < D; j++)
					{
						DiffT_[r][i][j].noalias() = CLM[i-1].transpose() * UM[i][j];
					}
				}else
				{
					for(int j = 0; j < D; j++)
					{
						DiffT_[r][i][j].noalias() = CLM[i-1].transpose() * UM[i][j] * CRM[i+1];
					}
				}
			}
			////////////////////////////////////////////////
			delete [] CRM;
			delete [] CLM;
			delete [] TM1;
			////////////////////////////////////////////////
		}else
		{
			norm = build_UpLow_MPS(PhyC, r-1, r+1);
			////////////////////////////////////////////////
			for(int i = 0; i < D*D; i++)
			{
				DiffT_[r][0][i].setZero();
				DiffT_[r][L-1][i].setZero();
				for(int j = 1; j < L-1; j++)
				{
					DiffT_[r][j][i].setZero();
				}
			}
			Mxd ** CRM = new Mxd * [L];
			Mxd ** CLM = new Mxd * [L];
			for(int i = 0; i < L; i++)  
			{
				CRM[i] = new Mxd [D];
				CLM[i] = new Mxd [D];
			}
			Mxd ** TM1 = new Mxd * [D];
			Mxd ** TM2 = new Mxd * [D];
			for(int i = 0; i < D; i++)
			{
				TM1[i] = new Mxd [D];
				TM2[i] = new Mxd [D];
			}
			// DimW is not referenced
			////////////////////////////////////////////////
			// Use UM[r-1], PEPS[r], LM[r+1] to build CRM
			// Last site //
			r1=UM[L-1][0].rows();
			c1=UM[L-1][0].cols();
			r2=LM[L-1][0].rows();
			c2=LM[L-1][0].cols();
			for(int i = 0; i < D; i++)
			{
				CRM[L-1][i].setZero(r1,r2);
			}
			for(int j = 0; j < D; j++)
			{
				for(int k = 0; k < D; k++)
				{
					Mxd tp;
					tp.setZero(c2,r2);
					for(int l = 0; l < D; l++)
					{
						tp += PEPS_[r][L-1][PhyC[r][L-1]][k*D+l](j,0) * LM[L-1][l].transpose();
					}
					CRM[L-1][j].noalias() +=  UM[L-1][k] * tp;
				}
			}
			// Last-but-1 to the first site //
			for (int i = L-2; i >= 1; i--) {
				r1=UM[i][0].rows();
				c1=UM[i][0].cols();
				r2=LM[i][0].rows();
				c2=LM[i][0].cols();
				for(int j = 0; j < D; j++)
				{
					CRM[i][j].setZero(r1,r2);
					for(int k = 0; k < D; k++)
					{
						TM1[k][j].setZero(c1,r2);
						TM2[k][j].setZero(c1,r2);
					}
				}
				for(int tid = 0; tid < D; tid++)
				{
					for(int j = 0; j < D; j++)
					{
						TM1[tid][j].noalias() = CRM[i+1][j] * LM[i][tid].transpose();
					}
				}
				for(int tid = 0; tid < D; tid++)
				{
					for(int j = 0; j < D; j++)
					{
						for(int l = 0; l < D; l++)
						{
							for(int k = 0; k < D; k++)
							{
								TM2[tid][j] +=  PEPS_[r][i][PhyC[r][i]][tid*D+k](j,l) * TM1[k][l];
							}
						}
					}
				}
				for(int j = 0; j < D; j++)
				{
					for(int k = 0; k < D; k++)
					{
						CRM[i][j].noalias() += UM[i][k] * TM2[k][j];
					}
				}
			}
			////////////////////////////////////////////////
			// Use UM[r-1], PEPS[r], LM[r+1] to build CLM
			// 1st site //
			r1=UM[0][0].rows();
			c1=UM[0][0].cols();
			r2=LM[0][0].rows();
			c2=LM[0][0].cols();
			for(int i = 0; i < D; i++)
			{
				CLM[0][i].setZero(c1,c2);
			}
			for(int j = 0; j < D; j++)
			{
				for(int k = 0; k < D; k++)
				{
					Mxd tp;
					tp.setZero(r2,c2);
					for(int l = 0; l < D; l++)
					{
						tp += PEPS_[r][0][PhyC[r][0]][k*D+l](0,j) * LM[0][l];
					}
					CLM[0][j].noalias() +=  UM[0][k].transpose() * tp;
				}
			}
			// 2nd to the last site //
			for (int i = 1; i < L-1; i++) {
				r1=UM[i][0].rows();
				c1=UM[i][0].cols();
				r2=LM[i][0].rows();
				c2=LM[i][0].cols();
				for(int j = 0; j < D; j++)
				{
					CLM[i][j].setZero(c1,c2);
					for(int k = 0; k < D; k++)
					{
						TM1[k][j].setZero(r1,c2);
						TM2[k][j].setZero(r1,c2);
					}
				}
				for(int tid = 0; tid < D; tid++)
				{
					for(int j = 0; j < D; j++)
					{
						TM1[tid][j].noalias() = CLM[i-1][j] * LM[i][tid];
					}
				}
				for(int tid = 0; tid < D; tid++)
				{
					for(int j = 0; j < D; j++)
					{
						for(int l = 0; l < D; l++)
						{
							for(int k = 0; k < D; k++)
							{
								TM2[tid][j] +=  TM1[k][l] * PEPS_[r][i][PhyC[r][i]][tid*D+k](l,j) ;
							}
						}
					}
				}
				for(int j = 0; j < D; j++)
				{
					for(int k = 0; k < D; k++)
					{
						CLM[i][j].noalias() += UM[i][k].transpose() * TM2[k][j];
					}
				}
			}
			////////////////////////////////////////////////
			// Build the derivative matrix
			for(int i = 0; i < L; i++)
			{
				if(i==0)
				{
					for(int j = 0; j < D; j++)
					{
						for(int k = 0; k < D; k++)
						{
							for(int l = 0; l < D; l++)
							{
								DiffT_[r][i][j*D+k](0,l) = (UM[i][j] * CRM[i+1][l] * LM[i][k].transpose())(0,0);
							}
						}
					}
				}else if(i==L-1)
				{
					for(int j = 0; j < D; j++)
					{
						for(int k = 0; k < D; k++)
						{
							for(int l = 0; l < D; l++)
							{
								DiffT_[r][i][j*D+k](l,0) = (UM[i][j].transpose() * CLM[i-1][l] * LM[i][k])(0,0);
							}
						}
					}
				}else
				{
					for(int j = 0; j < D; j++)
					{
						for(int k = 0; k < D; k++)
						{
							for(int l = 0; l < D; l++)
							{
								for(int s = 0; s < D; s++)
								{
									DiffT_[r][i][j*D+k](l,s) = (UM[i][j].transpose() * CLM[i-1][l] * LM[i][k] * CRM[i+1][s].transpose()).trace();
								}
							}
						}
					}
				}
			}
			////////////////////////////////////////////////
			for(int i = 0; i < L; i++)
			{
				delete [] CRM[i];
				delete [] CLM[i];
			}
			for(int i = 0; i < D; i++)
			{
				delete [] TM1[i];
				delete [] TM2[i];
			}
			delete [] CRM;
			delete [] CLM;
			delete [] TM1;
			delete [] TM2;
			////////////////////////////////////////////////
		}
	}
}

void PEPS_Base::setRandom()
{
	if(initted_)
	{
		for(int i = 0; i < W; i++)
		{
			for(int j = 0; j < L; j++)
			{
				for(int k = 0; k < phyD; k++)
				{
					if(i==0||i==W-1)
					{
						for(int l = 0; l < D; l++)
						{
							PEPS_[i][j][k][l].setRandom();
						}
					}
					else
					{
						for(int l = 0; l < D*D; l++)
						{
							PEPS_[i][j][k][l].setRandom();
						}
					}
				}
			}
		}
	}
}

void PEPS_Base::setUniform()
{
	if(initted_)
	{
		for(int i = 0; i < W; i++)
		{
			for(int j = 0; j < L; j++)
			{
				for(int k = 0; k < phyD; k++)
				{
					if(i==0||i==W-1)
					{
						for(int l = 0; l < D; l++)
						{
							Mxd tp(PEPS_[i][j][k][l].rows(),PEPS_[i][j][k][l].cols());
							tp.setIdentity();
							PEPS_[i][j][k][l] = 0.37 * tp;
						}
					}
					else
					{
						for(int l = 0; l < D*D; l++)
						{
							Mxd tp(PEPS_[i][j][k][l].rows(),PEPS_[i][j][k][l].cols());
							tp.setIdentity();
							PEPS_[i][j][k][l] = 0.37 * tp;
						}
					}
				}
			}
		}
	}
}

void PEPS_Base::setNearUniform()
{
	if(initted_)
	{
		for(int i = 0; i < W; i++)
		{
			for(int j = 0; j < L; j++)
			{
				for(int k = 0; k < phyD; k++)
				{
					if(i==0||i==W-1)
					{
						for(int l = 0; l < D; l++)
						{
							Mxd tp(PEPS_[i][j][k][l].rows(),PEPS_[i][j][k][l].cols());
							tp.setIdentity();
							PEPS_[i][j][k][l].setRandom();
							PEPS_[i][j][k][l] = 0.1 * PEPS_[i][j][k][l] + 0.37 * tp;
						}
					}
					else
					{
						for(int l = 0; l < D*D; l++)
						{
							Mxd tp(PEPS_[i][j][k][l].rows(),PEPS_[i][j][k][l].cols());
							tp.setIdentity();
							PEPS_[i][j][k][l].setRandom();
							PEPS_[i][j][k][l] = 0.1 * PEPS_[i][j][k][l] + 0.37 * tp;
						}
					}
				}
			}
		}
	}
}

void PEPS_Base::setProductState()
{
	if(initted_)
	{
		for(int i = 0; i < W; i++)
		{
			for(int j = 0; j < L; j++)
			{
				for(int k = 0; k < phyD; k++)
				{
					if(i==0||i==W-1)
					{
						for(int l = 0; l < D; l++)
						{
							Mxd tp(PEPS_[i][j][k][l].rows(),PEPS_[i][j][k][l].cols());
							if( (i+j)%2==0 && k==1 )
								tp.setIdentity();
							else if( (i+j)%2==1 && k==2 )
								tp.setIdentity();
							else
								tp.setZero();
							PEPS_[i][j][k][l] = 0.37 * tp;
						}
					}
					else
					{
						for(int l = 0; l < D*D; l++)
						{
							Mxd tp(PEPS_[i][j][k][l].rows(),PEPS_[i][j][k][l].cols());
							if( (i+j)%2==0 && k==1 )
								tp.setIdentity();
							else if( (i+j)%2==1 && k==2 )
								tp.setIdentity();
							else
								tp.setZero();
							PEPS_[i][j][k][l] = 0.37 * tp;
						}
					}
				}
			}
		}
	}
}

void PEPS_Base::setNearProductState()
{
	if(initted_)
	{
		for(int i = 0; i < W; i++)
		{
			for(int j = 0; j < L; j++)
			{
				for(int k = 0; k < phyD; k++)
				{
					if(i==0||i==W-1)
					{
						for(int l = 0; l < D; l++)
						{
							Mxd tp(PEPS_[i][j][k][l].rows(),PEPS_[i][j][k][l].cols());
							if( (i+j)%2==0 && k==1 )
								tp.setIdentity();
							else if( (i+j)%2==1 && k==2 )
								tp.setIdentity();
							else
								tp.setZero();
							PEPS_[i][j][k][l].setRandom();
							PEPS_[i][j][k][l] = 0.36 * PEPS_[i][j][k][l] + 0.37 * tp;
						}
					}
					else
					{
						for(int l = 0; l < D*D; l++)
						{
							Mxd tp(PEPS_[i][j][k][l].rows(),PEPS_[i][j][k][l].cols());
							if( (i+j)%2==0 && k==1 )
								tp.setIdentity();
							else if( (i+j)%2==1 && k==2 )
								tp.setIdentity();
							else
								tp.setZero();
							PEPS_[i][j][k][l].setRandom();
							PEPS_[i][j][k][l] = 0.36 * PEPS_[i][j][k][l] + 0.37 * tp;
						}
					}
				}
			}
		}
	}
}

void PEPS_Base::M_setZero(int phyDim, int Dim, Mxd ** MPS)
{
	for(int i = 0; i < phyDim; i++)
	{
		MPS[0][i].setZero(1,Dim);
		for(int j = 1; j < L-1; j++)
		{
			MPS[j][i].setZero(Dim,Dim);
		}
		MPS[L-1][i].setZero(Dim,1);
	}
}

void PEPS_Base::M_setRandom(int phyDim, int Dim, Mxd ** MPS)
{
	for(int i = 0; i < phyDim; i++)
	{
		MPS[0][i].setRandom(1,Dim);
		for(int j = 1; j < L-1; j++)
		{
			MPS[j][i].setRandom(Dim,Dim);
		}
		MPS[L-1][i].setRandom(Dim,1);
	}
}

double PEPS_Base::psiphi(int pD, Mxd ** Psi, Mxd ** Phi)
{
	int tid;
	int i, r1, r2, c1, c2;
	double EnergyExpec = 0.0;
	////////////////////////////////////////////
	Mxd * CRM = new Mxd [L];
	Mxd * TM1 = new Mxd [pD];
	////////////////////////////////////////////////
	// Last site //
	r1=Psi[L-1][0].rows();
	r2=Phi[L-1][0].rows();
	CRM[L-1].setZero(r1,r2);
	for(int i = 0; i < pD; i++)
	{
		CRM[L-1] += Psi[L-1][i] * Phi[L-1][i].transpose();
	}
	// Last-but-1 to the second site //
	for (i = L-2; i >= 0; i--) {
		r1=Psi[i][0].rows();
		c1=Psi[i][0].cols();
		r2=Phi[i][0].rows();
		c2=Phi[i][0].cols();
		CRM[i].setZero(r1,r2);
		for(int j = 0; j < pD; j++)
		{
			TM1[j].setZero(c1,r2);
		}
		for(tid = 0; tid < pD; tid++)
		{
			// tid = omp_get_thread_num();
			TM1[tid].noalias() = CRM[i+1] * Phi[i][tid].transpose();
		}
		for(int j = 0; j < pD; j++)
		{
			CRM[i] += Psi[i][j] * TM1[j];
		}
	}
	EnergyExpec = CRM[0](0,0);
	////////////////////////////////////////////////
	delete [] CRM;
	delete [] TM1;
	////////////////////////////////////////////////
	return EnergyExpec;
}

double PEPS_Base::RC(int Dim, int pD, Mxd ** MPS)
{
	// (QR factorization) //
	Mxd TM(Dim,pD), Q, R;
	Mxd tempM;
	////////////////////////////////////
	// Last site //
	for (int i = 0; i < pD; i++) {
		TM.block(0,i,Dim,1)=MPS[L-1][i];
	}
	TM.transposeInPlace();
	rQR(&TM,&Q);
	Q.transposeInPlace();
	if (Dim>=pD) {
		R = Q * TM;
	}else {
		tempM = Q * TM;
		R = tempM.block(0,0,Dim,Dim);
	}
	
	if (Dim>=pD) {
		for (int i = 0; i < pD; i++) {
			MPS[L-1][i].block(0,0,pD,1) = Q.block(0,i,pD,1);
		}
		for (int i = pD; i < Dim; i++) {
			for (int j = 0; j < pD; j++) {
				MPS[L-1][j](i,0)=0;
			}
		}
		for (int i = 0; i < pD; i++) {
			tempM = MPS[L-2][i] * R.transpose();
			MPS[L-2][i] = Mxd::Zero(Dim,Dim);
			MPS[L-2][i].block(0,0,Dim,pD)=tempM;
		}
	}else {
		for (int i = 0; i < pD; i++) {
			MPS[L-1][i].block(0,0,Dim,1) = Q.block(0,i,Dim,1);
		}
		for (int i = 0; i < pD; i++) {
			tempM = MPS[L-2][i] * R.transpose();
			MPS[L-2][i] = Mxd::Zero(Dim,Dim);
			MPS[L-2][i].block(0,0,Dim,Dim)=tempM;
		}
	}
	// Last-but-1 to the second site //
	for (int i = L-2; i >= 1; i--) {
		TM.resize(Dim,pD*Dim);
		for (int j = 0; j < pD; j++) {
			TM.block(0,j*Dim,Dim,Dim)=MPS[i][j];
		}
		
		TM.transposeInPlace();
		rQR(&TM,&Q);
		Q.transposeInPlace();
		R = Q * TM;
		for (int j = 0; j < pD; j++) {
			MPS[i][j] = Q.block(0,j*Dim,Dim,Dim);
		}
		
		tempM = R.block(0,0,Dim,Dim);
		R = Mxd::Zero(Dim,Dim);
		R = tempM;
		for (int j = 0; j < pD; j++) {
			tempM = MPS[i-1][j] * R.transpose();
			MPS[i-1][j] = tempM;
		}
	}
	
	// 1st site //
	TM.resize(1,pD*Dim);
	for (int i = 0; i < pD; i++) {
		TM.block(0,i*Dim,1,Dim)=MPS[0][i];
	}
	
	TM.transposeInPlace();
	rQR(&TM,&Q);
	Q.transposeInPlace();
	R = Q * TM;
	// cout<<"Norm(1): "<<R(0,0)<<endl;
	for (int i = 0; i < pD; i++) {
		MPS[0][i].block(0,0,1,Dim) = Q.block(0,i*Dim,1,Dim);
	}
	if (R(0,0)<0) {
		for (int i = 0; i < pD; i++) {
			MPS[0][i] = -MPS[0][i];
		}
		return -R(0,0);
	}else
	{
		return R(0,0);
	}
}

void PEPS_Base::CompressL(int Dim, int DT, int pD, Mxd ** VMPS, Mxd ** MPS)
{
	int tid;
	double * sv = new double [pD*Dim];
	int row, col;
	int ActualD = 1;
	int lD;
	/////////////////////////////////////////////
	Mxd TM, U, V;
	///////////////////////////////////////////
	for(int i = 0; i < L; i++)
	{
		row=VMPS[i][0].rows();
		col=VMPS[i][0].cols();
		TM.resize(pD*ActualD,col);
		for(tid = 0; tid < pD; tid++)
		{
			// tid = omp_get_thread_num();
			TM.block(tid*ActualD,0,ActualD,col)=VMPS[i][tid].block(0,0,ActualD,col);
		}
		lD=rDenMatDecomp(&TM,pD*Dim,sv,&U,&V,DT,SVD_Tolerance);
		for(tid = 0; tid < pD; tid++)
		{
			// tid = omp_get_thread_num();
			MPS[i][tid].setZero(row,min(DT,col));
			MPS[i][tid].real().block(0,0,ActualD,U.cols())=U.block(tid*ActualD,0,ActualD,U.cols());
			if(i!=L-1)
			{
				Mxd tempM;
				tempM.noalias()=V*VMPS[i+1][tid].real();
				VMPS[i+1][tid].setZero(DT,tempM.cols());
				VMPS[i+1][tid].real().block(0,0,tempM.rows(),tempM.cols())=tempM;
			}
		}
		ActualD=min(lD,int(U.cols()));
	}
	delete [] sv;
}

double PEPS_Base::IterCompress (int IfRandom, int MaxIter, double tol, int VDim, int Dim, int pD, Mxd ** VMPS, Mxd ** MPS)
{
	double dt; //Distance from target MPS
	int i, r1, r2, c1, c2;
	////////////////////////////////////////////
	Mxd * CRM = new Mxd [L];
	Mxd * CLM = new Mxd [L];
	Mxd * TM1 = new Mxd [pD];
	////////////////////////////////////////////////
	// Start from Random Compressed MPS //
	if (IfRandom!=0) {
		for (int i = 1; i < L-1; i++) {
			// Physical degrees of freedom: (0:+, 1:-)
			for (int j = 0; j < pD; j++) {
				MPS[i][j].setRandom();
			}
		}
		for (int i = 0; i < pD; i++) {
			MPS[0][i].setRandom();
			MPS[L-1][i].setRandom();
		}
	}
	//Right Canonicalize the MPS and Build CR to the left
	// (QR factorization) //
	Mxd TM(Dim,pD), Q, R;
	Mxd tempM;
	////////////////////////////////////
	// Last site //
	for (int i = 0; i < pD; i++) {
		TM.block(0,i,Dim,1)=MPS[L-1][i];
	}
	TM.transposeInPlace();
	rQR(&TM,&Q);
	Q.transposeInPlace();
	if (Dim>=pD) {
		R = Q * TM;
	}else {
		tempM = Q * TM;
		R = tempM.block(0,0,Dim,Dim);
	}
	if (pD>Dim) {
		for (int i = 0; i < pD; i++) {
			MPS[L-1][i].block(0,0,Dim,1) = Q.block(0,i,Dim,1);
		}
		for (int i = 0; i < pD; i++) {
			tempM = MPS[L-2][i] * R.transpose();
			MPS[L-2][i] = Mxd::Zero(Dim,Dim);
			MPS[L-2][i].block(0,0,Dim,Dim)=tempM;
		}
	}else {
		for (int i = 0; i < pD; i++) {
			MPS[L-1][i].block(0,0,pD,1) = Q.block(0,i,pD,1);
		}
		for (int i = pD; i < Dim; i++) {
			for (int j = 0; j < pD; j++) {
				MPS[L-1][j](i,0)=0;
			}
		}
		for (int i = 0; i < pD; i++) {
			tempM = MPS[L-2][i] * R.transpose();
			MPS[L-2][i] = Mxd::Zero(Dim,Dim);
			MPS[L-2][i].block(0,0,Dim,pD)=tempM;
		}
	}
	// Building the CR //
	r1=MPS[L-1][0].rows();
	r2=VMPS[L-1][0].rows();
	CRM[L-1].setZero(r2,r1);
	for(int i = 0; i < pD; i++)
	{
		CRM[L-1] += VMPS[L-1][i] * MPS[L-1][i].transpose();
	}
	// Last-but-1 to the second site //
	for (int i = L-2; i >= 1; i--) {
		// QR factorization //
		TM.resize(Dim,pD*Dim);
		for (int j = 0; j < pD; j++) {
			TM.block(0,j*Dim,Dim,Dim)=MPS[i][j];
		}
		
		TM.transposeInPlace();
		rQR(&TM,&Q);
		Q.transposeInPlace();
		R = Q * TM;
		
		for (int j = 0; j < pD; j++) {
			MPS[i][j] = Q.block(0,j*Dim,Dim,Dim);
		}
		
		tempM = R.block(0,0,Dim,Dim);
		R = Mxd::Zero(Dim,Dim);
		R = tempM;
		for (int j = 0; j < pD; j++) {
			MPS[i-1][j] = MPS[i-1][j] * R.transpose();
		}
		// Building the CR //
		r1=MPS[i][0].rows();
		c1=MPS[i][0].cols();
		r2=VMPS[i][0].rows();
		c2=VMPS[i][0].cols();
		CRM[i].setZero(r2,r1);
		for(int j = 0; j < pD; j++)
		{
			TM1[j].setZero(r2,c1);
		}
		for(int j = 0; j < pD; j++)
		{
			TM1[j].noalias() = VMPS[i][j] * CRM[i+1];
		}
		for(int j = 0; j < pD; j++)
		{
			CRM[i] += TM1[j] * MPS[i][j].transpose();
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////
	// Sweeping and improve the compressed MPS (QR factorization) //
	for (int l = 0; l < MaxIter; l++) {
		// Sweep from the left //
		// First site //
		// cout<<"Sweeping to the right."<<endl;
		// Variation //
		for(int i = 0; i < pD; i++)
		{
			MPS[0][i].setZero();
		}
		for(int i = 0; i < pD; i++)
		{
			MPS[0][i].noalias() = VMPS[0][i] * CRM[1];
		}
		/////
		dt = 1.0;
		for(int i = 0; i < pD; i++)
		{
			dt -= MPS[0][i].squaredNorm();
		}
		// cout<<"Distance:"<<dt<<endl;
		if (abs(dt)<=tol) {
//			cout<<"Distance:"<<dt<<endl;
			l=MaxIter;
			// break;
			return dt;
		}
		/////
		// QR factorization //
		TM.resize(pD,Dim);
		for (int i = 0; i < pD; i++) {
			TM.block(i,0,1,Dim)=MPS[0][i];
		}
		rQR(&TM,&Q);

		if (pD>Dim) {
			for (int i = 0; i < pD; i++) {
				MPS[0][i].block(0,0,1,Dim) = Q.block(i,0,1,Dim);
			}
		}else {
			for (int i = 0; i < pD; i++) {
				MPS[0][i].block(0,0,1,pD) = Q.block(i,0,1,pD);
			}
			for (int i = pD; i < Dim; i++) {
				for (int j = 0; j < pD; j++) {
					MPS[0][j](0,i)=0;
				}
			}
		}
		// Building the CL //
		c1=MPS[0][0].cols();
		c2=VMPS[0][0].cols();
		CLM[0].setZero(c1,c2);
		for(int i = 0; i < pD; i++)
		{
			CLM[0] += MPS[0][i].transpose() * VMPS[0][i];
		}
		// the 2nd to last-but 1 site //
		for (int i = 1; i < L-1; i++) {
			// Variation //
			for (int j = 0; j < pD; j++) {
				MPS[i][j].setZero();
			}
			for(int j = 0; j < pD; j++)
			{
				MPS[i][j] = CLM[i-1] * VMPS[i][j] * CRM[i+1];
			}
			
			dt = 1.0;
			for(int j = 0; j < pD; j++)
			{
				dt -= MPS[i][j].squaredNorm();
			}
			// cout<<"Distance:"<<dt<<endl;
			if (abs(dt)<=tol) {
//				cout<<"Distance:"<<dt<<endl;
				l=MaxIter;
				// break;
				return dt;
			}
			// QR factorization //
			TM.resize(pD*Dim,Dim);
			for (int j = 0; j < pD; j++) {
				TM.block(j*Dim,0,Dim,Dim)=MPS[i][j];
			}
			rQR(&TM,&Q);
			for (int j = 0; j < pD; j++) {
				MPS[i][j] = Q.block(j*Dim,0,Dim,Dim);
			}
			// Building the CL //
			r1=MPS[i][0].rows();
			c1=MPS[i][0].cols();
			r2=VMPS[i][0].rows();
			c2=VMPS[i][0].cols();
			CLM[i].setZero(c1,c2);
			for(int j = 0; j < pD; j++)
			{
				TM1[j].setZero(r1,c2);
			}
			for(int j = 0; j < pD; j++)
			{
				TM1[j].noalias() = CLM[i-1] * VMPS[i][j];
			}
			for(int j = 0; j < pD; j++)
			{
				CLM[i] += MPS[i][j].transpose() * TM1[j];
			}
		}
		if (l==MaxIter) {
			break;
		}
///////////////////////////////////////////////////////////////////////////////////////////
		// Sweep from the right //
//		cout<<"Sweeping to the left."<<endl;
		// Last site //
		// Variation //
		for (int i = 0; i < pD; i++) {
			MPS[L-1][i].setZero();
		}
		for(int i = 0; i < pD; i++)
		{
			MPS[L-1][i].noalias() = CLM[L-2] * VMPS[L-1][i];
		}
		
		dt = 1.0;
		for(int i = 0; i < pD; i++)
		{
			dt -= MPS[L-1][i].squaredNorm();
		}
		// cout<<"Distance:"<<dt<<endl;
		if (abs(dt)<=tol) {
//			cout<<"Distance:"<<dt<<endl;
			// break;
			return dt;
		}
		// QR factorization //
		TM.resize(Dim,pD);
		for (int i = 0; i < pD; i++) {
			TM.block(0,i,Dim,1)=MPS[L-1][i];
		}
		
		TM.transposeInPlace();
		rQR(&TM,&Q);
		Q.transposeInPlace();
		if (Dim<pD) {
			for (int i = 0; i < pD; i++) {
				MPS[L-1][i].block(0,0,Dim,1) = Q.block(0,i,Dim,1);
			}
		}else {
			for (int i = 0; i < pD; i++) {
				MPS[L-1][i].block(0,0,pD,1) = Q.block(0,i,pD,1);
			}
			for (int i = pD; i < Dim; i++) {
				for (int j = 0; j < pD; j++) {
					MPS[L-1][j](i,0)=0;
				}
			}
		}
		// Building the CR //
		r1=MPS[L-1][0].rows();
		r2=VMPS[L-1][0].rows();
		CRM[L-1].setZero(r2,r1);
		for(int i = 0; i < pD; i++)
		{
			CRM[L-1] += VMPS[L-1][i] * MPS[L-1][i].transpose();
		}
		// Sweep from the right -- last-but 1 to the 2nd site //
		for (int i = L-2; i >= 1; i--) {
			// Variation //
			for (int j = 0; j < pD; j++) {
				MPS[i][j].setZero();
			}
			for(int j = 0; j < pD; j++)
			{
				MPS[i][j] = CLM[i-1] * VMPS[i][j] * CRM[i+1];
			}
			
			dt = 1.0;
			for(int j = 0; j < pD; j++)
			{
				dt -= MPS[i][j].squaredNorm();
			}
			// cout<<"Distance:"<<dt<<endl;
			if (abs(dt)<=tol) {
//				cout<<"Distance:"<<dt<<endl;
				l=MaxIter;
				// break;
				return dt;
			}
			// QR factorization //
			TM.resize(Dim,pD*Dim);
			for (int j = 0; j < pD; j++) {
				TM.block(0,j*Dim,Dim,Dim)=MPS[i][j];
			}
			
			TM.transposeInPlace();
			rQR(&TM,&Q);
			Q.transposeInPlace();
			for (int j = 0; j < pD; j++) {
				MPS[i][j] = Q.block(0,j*Dim,Dim,Dim);
			}
			// Building the CR //
			r1=MPS[i][0].rows();
			c1=MPS[i][0].cols();
			r2=VMPS[i][0].rows();
			c2=VMPS[i][0].cols();
			CRM[i].setZero(r2,r1);
			for(int j = 0; j < pD; j++)
			{
				TM1[j].setZero(r2,c1);
			}
			for(int j = 0; j < pD; j++)
			{
				TM1[j].noalias() = VMPS[i][j] * CRM[i+1];
			}
			for(int j = 0; j < pD; j++)
			{
				CRM[i] += TM1[j] * MPS[i][j].transpose();
			}
		}
	}
	delete [] CRM;
	delete [] CLM;
	delete [] TM1;
}

#endif
