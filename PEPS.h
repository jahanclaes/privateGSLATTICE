#ifndef PEPS_H
#define PEPS_H
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>

#include "Random/Random.h"
#include "Blitz.h"
#include "SystemClass.h"
#include "WaveFunction.h"
#include "PairingFunctionMany.h"

#include "PEPS_Base_Class/PEPS_Base.h"
//#include "PEPS_Base_Class/PEPS_Base.cpp"

class PEPSClass  : public WaveFunctionClass
{
public:

  int L_;
  int W_;
  int D_;
  int VD_;
  PEPS_Base * peps;

  int NumSpinUp;
  PEPSClass() 
  {
  }

  complex<double> logevaluate(SystemClass &system,int &sign);

  //currently only the size matters
  void Swap(int i, int j)
  {
  }

  void Move(int site, int end_site, int spin)
  {
    
  }
  void Reject(SystemClass &system,int site,int end_site,int spin)
  {
    
  }

  void AllDerivs(SystemClass &system, Array<complex<double>,1>  &derivs);
  void AllDerivs(SystemClass &system, Array<complex<double>,1> &derivs,int start,int stop);
  void RealDerivs(SystemClass &system, Array<complex<double>,1> &derivs)  ;


  void  CheckDerivs(SystemClass &system, Array<complex<double>,1>  &derivs,int start, int stop);
    
  void Init(SystemClass &system,int tL_, int tW_, int tD_, int tVD_)
  {
    Name="PEPS";
    NeedFrequentReset=false;
    NumSpinUp=system.x.size()/2;
    bool ReadParams=false;
    if (ReadParams){
    }
    //SET ME CORRECTLY!  IMPLEMENT!
    //Need to find out how to read paramters
        L_ = tL_;
        W_ = tW_;
        D_ = tD_;
        VD_ = tVD_;//used to be 25

    //peps = new PEPS_Base(W, L, phyD, D, VD);
    peps = new PEPS_Base(W_, L_, 4, D_, VD_);
	
	// choose a starting configuration
    //	peps->setNearProductState();
    peps->setProductState();
	// peps->setNearUniform();
	// peps->setUniform();

	
    int N_edge = 4*(2*D_*D_+D_*D_*D_*(L_-2));
    int N_mid  = 4*(2*D_*D_*D_+D_*D_*D_*D_*(L_-2));
    NumParams=2*N_edge + (W_-2)*N_mid;
  }
  double GetParam_real(int i)
  {
    //IMPLEMENT ME!
    int N_edge = 4*(2*D_*D_+D_*D_*D_*(L_-2));
    int N_mid  = 4*(2*D_*D_*D_+D_*D_*D_*D_*(L_-2));
    int r_, c_, phy_, s_, l_;
    if( i<0 || i>=(2*N_edge+(W_-2)*N_mid) )
      {
	cout<<"GetParam_real(int i): "<<i <<" "<<(2*N_edge+(W_-2)*N_mid)<<"out of bound!"<<endl;
	exit(0);
      }
    else if(i<N_edge)
      {
	r_ = 1;
	if(i<4*D_*D_)
	  {
	    c_ = 1;
	    phy_ = i/(D_*D_);
	    s_ = (i%(D_*D_))/D_;
	    l_ = (i%(D_*D_))%D_;
	  }
	else if( i>=(N_edge-4*D_*D_) )
	  {
            c_ = L_ - 1;
            phy_ = (i+4*D_*D_-N_edge)/(D_*D_);
            s_ = ((i+4*D_*D_-N_edge)%(D_*D_))/D_;
            l_ = ((i+4*D_*D_-N_edge)%(D_*D_))%D_;
	  }
	else
	  {
            c_ = (i-4*D_*D_)/(4*D_*D_*D_) + 1;
            phy_ = (i-(c_-1)*4*D_*D_*D_-4*D_*D_)/(D_*D_*D_);
            s_ = ((i-(c_-1)*4*D_*D_*D_-4*D_*D_)%(D_*D_*D_))/(D_*D_);
            l_ = ((i-(c_-1)*4*D_*D_*D_-4*D_*D_)%(D_*D_*D_))%(D_*D_);
	  }
      }
    else if( i>=N_edge && i<(NumParams-N_edge) )
      {
	r_ = (i-N_edge)/N_mid + 1;
	int ii = (i-N_edge)%N_mid;
	if(ii<4*D_*D_*D_)
          {
            c_ = 1;
            phy_ = ii/(D_*D_*D_);
            s_ = (ii%(D_*D_*D_))/D_;
            l_ = (ii%(D_*D_*D_))%D_;
          }
	else if( ii>=(N_mid-4*D_*D_*D_) )
          {
            c_ = L_ - 1;
            phy_ = (ii+4*D_*D_*D_-N_mid)/(D_*D_*D_);
            s_ = ((ii+4*D_*D_*D_-N_mid)%(D_*D_*D_))/D_;
            l_ = ((ii+4*D_*D_*D_-N_mid)%(D_*D_*D_))%D_;
          }
        else
          {
            c_ = (ii-4*D_*D_*D_)/(4*D_*D_*D_*D_) + 1;
            phy_ = (ii-(c_-1)*4*D_*D_*D_*D_-4*D_*D_*D_)/(D_*D_*D_*D_);
            s_ = ((ii-(c_-1)*4*D_*D_*D_*D_-4*D_*D_*D_)%(D_*D_*D_*D_))/(D_*D_);
            l_ = ((ii-(c_-1)*4*D_*D_*D_*D_-4*D_*D_*D_)%(D_*D_*D_*D_))%(D_*D_);
          }
      }
    else if( i>=(NumParams-N_edge) )
      {
        r_ = W_ - 1;
        int ii = i+N_edge-NumParams;
        if(ii<4*D_*D_)
          {
            c_ = 1;
            phy_ = ii/(D_*D_);
            s_ = (ii%(D_*D_))/D_;
            l_ = (ii%(D_*D_))%D_;
          }
        else if( ii>=(N_edge-4*D_*D_) )
          {
            c_ = L_ - 1;
            phy_ = (ii+4*D_*D_-N_edge)/(D_*D_);
            s_ = ((ii+4*D_*D_-N_edge)%(D_*D_))/D_;
            l_ = ((ii+4*D_*D_-N_edge)%(D_*D_))%D_;
          }
        else
          {
            c_ = (ii-4*D_*D_)/(4*D_*D_*D_) + 1;
            phy_ = (ii-(c_-1)*4*D_*D_*D_-4*D_*D_)/(D_*D_*D_);
            s_ = ((ii-(c_-1)*4*D_*D_*D_-4*D_*D_)%(D_*D_*D_))/(D_*D_);
            l_ = ((ii-(c_-1)*4*D_*D_*D_-4*D_*D_)%(D_*D_*D_))%(D_*D_);
          }
      }

    return peps->PEPS_[r_][c_][phy_][s_](l_);
  }
  void SetParam_real(int i, double param)
  {
    //IMPLEMENT ME!
    int N_edge = 4*(2*D_*D_+D_*D_*D_*(L_-2));
    int N_mid  = 4*(2*D_*D_*D_+D_*D_*D_*D_*(L_-2));
    int r_, c_, phy_, s_, l_;
    if( i<0 || i>=(2*N_edge+(W_-2)*N_mid) )
    {
      cout<<"SetParam_real(int i): "<<i <<" "<<(2*N_edge+(W_-2)*N_mid)<<"out of bound!"<<endl;
      //      cout<<"GetParam_real(int i): i out of bound!"<<endl;
      exit(0);
    }
    else if(i<N_edge)
      {
        r_ = 1;
        if(i<4*D_*D_)
          {
            c_ = 1;
            phy_ = i/(D_*D_);
            s_ = (i%(D_*D_))/D_;
            l_ = (i%(D_*D_))%D_;
          }
        else if( i>=(N_edge-4*D_*D_) )
          {
            c_ = L_ - 1;
            phy_ = (i+4*D_*D_-N_edge)/(D_*D_);
            s_ = ((i+4*D_*D_-N_edge)%(D_*D_))/D_;
            l_ = ((i+4*D_*D_-N_edge)%(D_*D_))%D_;
          }
        else
          {
            c_ = (i-4*D_*D_)/(4*D_*D_*D_) + 1;
            phy_ = (i-(c_-1)*4*D_*D_*D_-4*D_*D_)/(D_*D_*D_);
            s_ = ((i-(c_-1)*4*D_*D_*D_-4*D_*D_)%(D_*D_*D_))/(D_*D_);
            l_ = ((i-(c_-1)*4*D_*D_*D_-4*D_*D_)%(D_*D_*D_))%(D_*D_);
          }
      }
    else if( i>=N_edge && i<(NumParams-N_edge) )
      {
        r_ = (i-N_edge)/N_mid + 1;
        int ii = (i-N_edge)%N_mid;
        if(ii<4*D_*D_*D_)
          {
            c_ = 1;
            phy_ = ii/(D_*D_*D_);
            s_ = (ii%(D_*D_*D_))/D_;
            l_ = (ii%(D_*D_*D_))%D_;
          }
        else if( ii>=(N_mid-4*D_*D_*D_) )
          {
            c_ = L_ - 1;
            phy_ = (ii+4*D_*D_*D_-N_mid)/(D_*D_*D_);
            s_ = ((ii+4*D_*D_*D_-N_mid)%(D_*D_*D_))/D_;
            l_ = ((ii+4*D_*D_*D_-N_mid)%(D_*D_*D_))%D_;
          }
        else
          {
            c_ = (ii-4*D_*D_*D_)/(4*D_*D_*D_*D_) + 1;
            phy_ = (ii-(c_-1)*4*D_*D_*D_*D_-4*D_*D_*D_)/(D_*D_*D_*D_);
            s_ = ((ii-(c_-1)*4*D_*D_*D_*D_-4*D_*D_*D_)%(D_*D_*D_*D_))/(D_*D_);
            l_ = ((ii-(c_-1)*4*D_*D_*D_*D_-4*D_*D_*D_)%(D_*D_*D_*D_))%(D_*D_);
          }
      }
    else if( i>=(NumParams-N_edge) )
      {
        r_ = W_ - 1;
        int ii = i+N_edge-NumParams;
        if(ii<4*D_*D_)
          {
            c_ = 1;
            phy_ = ii/(D_*D_);
            s_ = (ii%(D_*D_))/D_;
            l_ = (ii%(D_*D_))%D_;
          }
        else if( ii>=(N_edge-4*D_*D_) )
          {
            c_ = L_ - 1;
            phy_ = (ii+4*D_*D_-N_edge)/(D_*D_);
            s_ = ((ii+4*D_*D_-N_edge)%(D_*D_))/D_;
            l_ = ((ii+4*D_*D_-N_edge)%(D_*D_))%D_;
          }
        else
          {
            c_ = (ii-4*D_*D_)/(4*D_*D_*D_) + 1;
            phy_ = (ii-(c_-1)*4*D_*D_*D_-4*D_*D_)/(D_*D_*D_);
            s_ = ((ii-(c_-1)*4*D_*D_*D_-4*D_*D_)%(D_*D_*D_))/(D_*D_);
            l_ = ((ii-(c_-1)*4*D_*D_*D_-4*D_*D_)%(D_*D_*D_))%(D_*D_);
          }
      }
    peps->PEPS_[r_][c_][phy_][s_](l_) = param;
  }

  double GetParam_imag(int i)
  {
    //IMPLEMENT ME!
    return 0;
  }
  void SetParam_imag(int i, double param)
  {
    //IMPLEMENT ME!
  }


  complex<double> evaluateRatio(SystemClass &system,int start, int stop, int spin);
  

  complex<double> evaluate(SystemClass &system);
  complex<double> evaluateRatio(SystemClass &system,int swap1, int swap2);
  complex<double> evaluateRatio_check(SystemClass &system, int swap1, int swap2);
  void UpdateDets(SystemClass &system,int swap1, int swap2)
  {
  }

  void UpdateDets(SystemClass &system,int site, int end_site,int spin)
  {

  }


};

#endif
