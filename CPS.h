#ifndef CPS_H
#define CPS_H
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>

#include "Random/Random.h"
#include "Blitz.h"
#include "SystemClass.h"
#include "WaveFunction.h"
#include "PairingFunctionMany.h"

class CPSClass  : public WaveFunctionClass
{
public:
  PairingFunctionMany &PF;
  int NumSpinUp;
  CPSClass(PairingFunctionMany &t_PF) : PF(t_PF)
  {
  }

  complex<double> logevaluate(SystemClass &system,int &sign);

  //currently only the size matters
  void Swap(int i, int j)
  {
  }
  void AllDerivs(SystemClass &system, Array<complex<double>,1>  &derivs);
  void AllDerivs(SystemClass &system, Array<complex<double>,1> &derivs,int start,int stop);
  void RealDerivs(SystemClass &system, Array<complex<double>,1> &derivs)  ;


  void  CheckDerivs(SystemClass &system, Array<complex<double>,1>  &derivs,int start, int stop);
    
  void Init(SystemClass &system)
  {
    Name="CPS";
    NeedFrequentReset=false;
    NumSpinUp=system.x.size()/2;
    PF.ProcessSystem(system);

    bool ReadParams=false;
    if (ReadParams){
    }
    NumParams=PF.NumParams; //Correlators;
  }
  double GetParam_real(int i)
  {
    pair<int,int> myLoc=PF.ParamLoc[i];
    return PF.f0[myLoc.first][myLoc.second].real();
  }
  void SetParam_real(int i,double param)
  {
    pair<int,int> myLoc=PF.ParamLoc[i];
    PF.f0[myLoc.first][myLoc.second].real()=param;
  }

  double GetParam_imag(int i)
  {
    pair<int,int> myLoc=PF.ParamLoc[i];
    return PF.f0[myLoc.first][myLoc.second].imag();
  }
  void SetParam_imag(int i,double param)
  {
    pair<int,int> myLoc=PF.ParamLoc[i];
    PF.f0[myLoc.first][myLoc.second].imag()=param;
  }

  
  complex<double> evaluate(SystemClass &system);
  complex<double> evaluateRatio(SystemClass &system,int swap1, int swap2);
  complex<double> evaluateRatioFlip(SystemClass &system,int swap1);
  complex<double> evaluateRatio_check(SystemClass &system, int swap1, int swap2);
  void UpdateDets(SystemClass &system,int swap1, int swap2)
  {
  }

};

#endif
