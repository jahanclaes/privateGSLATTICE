#ifndef cPEPS_H
#define cPEPS_H

#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>

#include "Random/Random.h"
#include "Blitz.h"
#include "SystemClass.h"
#include "WaveFunction.h"
#include "PairingFunctionMany.h"

#include "cPEPS_Base_Class/cpeps.h"

class cPEPSClass  : public WaveFunctionClass
{
public:

  int xL;        // x-direction is periodic
  int yL;        // y-direction is open
  int pD;        // physical dimension -- 4 for Hubbard
  int xBD;       // bond dimnesion in x-direction
  int yBD;       // bond dimnesion in y-direction
  int max_yBD;   // max bond dimnesion in y-direction (beyond which trucation is performed)
  double tol;    // error tollerance of the truncation scheme
  cPEPS * peps;

  int NumSpinUp;
  cPEPSClass() 
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
    
  void Init(SystemClass &system, int xlen, int ylen, int phy, int xd, int yd, int maxd, double err_tol)
  {
    Name="cPEPS";
    NeedFrequentReset=false;
    NumSpinUp=system.x.size()/2;
    bool ReadParams=false;
    if (ReadParams){
    }
    //SET ME CORRECTLY!  IMPLEMENT!
    //Need to find out how to read paramters
	xL      = xlen;
	yL      = ylen;
	pD      = phy;
	xBD     = xd;
	yBD     = yd;
	max_yBD = maxd;
	tol     = err_tol;

    peps = new cPEPS(xL,yL,pD,xBD,yBD,max_yBD,tol);
	
	// choose a starting configuration
    //	peps->setNearProductState();
    // peps->setProductState();
	// peps->setNearUniform();
	peps->setUniform();

    NumParams=peps->numParams;
  }
  
  double GetParam_real(int i)
  {
	  //IMPLEMENT ME!
	  assert(i>=0);
	  
	  int x   = i / (peps->numParams/xL);
	  int idx = i - x * (peps->numParams/xL);
	  int y   = 0;
	  int pT  = 0;
	  int pF  = 0;
	  
	  for(int j = 0; j < yL; ++j)
	  {
		  if(idx >= pD*xBD*xBD*peps->TN[x].Dim[j]*peps->TN[x].Dim[j+1])
		  {
			  idx -= pD*xBD*xBD*peps->TN[x].Dim[j]*peps->TN[x].Dim[j+1];
		  }else
		  {
			  y    = j;
			  pT   = idx / (xBD*xBD*peps->TN[x].Dim[j]*peps->TN[x].Dim[j+1]);
			  idx -= pT * (xBD*xBD*peps->TN[x].Dim[j]*peps->TN[x].Dim[j+1]);
			  pF   = idx / (peps->TN[x].Dim[j]*peps->TN[x].Dim[j+1]);
			  idx -= pF * (peps->TN[x].Dim[j]*peps->TN[x].Dim[j+1]);
			  break;
		  }
	  }
	  
	  return peps->TN[x].T[y][pT][pF](idx);
  }
  
  void SetParam_real(int i, double param)
  {
	  //IMPLEMENT ME!
	  assert(i>=0);
	  
	  int x   = i / (peps->numParams/xL);
	  int idx = i - x * (peps->numParams/xL);
	  int y   = 0;
	  int pT  = 0;
	  int pF  = 0;
	  
	  for(int j = 0; j < yL; ++j)
	  {
		  if(idx >= pD*xBD*xBD*peps->TN[x].Dim[j]*peps->TN[x].Dim[j+1])
		  {
			  idx -= pD*xBD*xBD*peps->TN[x].Dim[j]*peps->TN[x].Dim[j+1];
		  }else
		  {
			  y    = j;
			  pT   = idx / (xBD*xBD*peps->TN[x].Dim[j]*peps->TN[x].Dim[j+1]);
			  idx -= pT * (xBD*xBD*peps->TN[x].Dim[j]*peps->TN[x].Dim[j+1]);
			  pF   = idx / (peps->TN[x].Dim[j]*peps->TN[x].Dim[j+1]);
			  idx -= pF * (peps->TN[x].Dim[j]*peps->TN[x].Dim[j+1]);
			  break;
		  }
	  }
	  peps->TN[x].T[y][pT][pF](idx) = param;
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
