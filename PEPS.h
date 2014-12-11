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

class PEPSClass  : public WaveFunctionClass
{
public:
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
    
  void Init(SystemClass &system)
  {
    Name="PEPS";
    NeedFrequentReset=false;
    NumSpinUp=system.x.size()/2;
    bool ReadParams=false;
    if (ReadParams){
    }
    //SET ME CORRECTLY!  IMPLEMENT!
    NumParams=0;
  }
  double GetParam_real(int i)
  {
    //IMPLEMENT ME!
  }
  void SetParam_real(int i,double param)
  {
    //IMPLEMENT ME!
  }

  double GetParam_imag(int i)
  {
    //IMPLEMENT ME!
  }
  void SetParam_imag(int i,double param)
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
