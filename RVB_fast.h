///TEST TEST TEST
#ifndef RVB_FAST_H
#define RVB_FAST_H
#include <iostream>
#include "Blitz.h"
#include "Random/Random.h"
#include <vector>
#include <algorithm>
#include <fstream>
#include "SystemClass.h"
#include "WaveFunction.h"
#include "SmartMatrix.h"
#include "PairingFunctionAllBin.h"


class RVBFastPsiClass  : public WaveFunctionClass
{
public:
  SmartMatrix mat;
  bool ReadParams;
  bool ReadPairingFunction;
  //  RVBFastPsiClass::
 RVBFastPsiClass(PairingFunctionAllBin &pf)  : PairingFunction(pf) {
    ReadParams=false;
    ReadPairingFunction=false;
    Name="RVBFastPsi";
  }
  vector<Array<complex<double>  ,1> > newCols;
  vector<Array<complex<double>  ,1> > newRows;
  vector<int> colIndices;
  vector<int> rowIndices;
  void Copy(WaveFunctionClass* wf)
  {
    RVBFastPsiClass &b(*(RVBFastPsiClass*)wf);
    mat=b.mat;
    newCols=b.newCols;
    newRows=b.newRows;
    colIndices=b.colIndices;
    rowIndices=b.rowIndices;
    NumSpinUp=b.NumSpinUp;
    rebuild=b.rebuild;
    

  }
  //  PairingFunctionVectorClass PairingFunction; 
  PairingFunctionAllBin &PairingFunction; 

  int NumSpinUp;
  bool rebuild;
  Array<complex<double>,1> u;  Array<complex<double>,1> up;
  void SetParams(int i,double delta3, SystemClass &system);

  void AllDerivs(SystemClass &system, Array<complex<double>,1> &derivs,int start,int stop);
  void AllDerivs(SystemClass &system, Array<complex<double>,1> &derivs);

  double GetParam_real(int i);
  void SetParam_real(int i, double param);

  double GetParam_imag(int i);
  void SetParam_imag(int i, double param);

  complex<double> Deriv(SystemClass &system,int bin); 

  void Init(SystemClass &system);
  void Swap(int i, int j);
  void SetParams(double delta3,SystemClass &system);
  
  complex<double>Phi(int i,int j,SystemClass &system);

  void FillDet(SystemClass &system, SmartMatrix &myMat);
  complex<double> evaluate(SystemClass &system);
  complex<double> evaluate_noInverse(SystemClass &system);
  complex<double> evaluateRatio(SystemClass &system,int swap1, int swap2);
  complex<double> evaluateRatio_energy(SystemClass &system,int swap1, int swap2);

  double Sign(SystemClass &system);

  complex<double> evaluateRatio_honeycomb(SystemClass &system,TinyVector<int,6>  &honeycomb_locs,
					  TinyVector<int,6> &honeycomb_backup,int amt);
  complex<double> evaluateRatio_swap(SystemClass &system, int swap1, int swap2);

  void Reject(SystemClass &system,int swap1,int swap2);
  void UpdateDets(SystemClass &system,int swap1, int swap2);



  //For debugging
  SmartMatrix mat_check;
  complex<double> logevaluate(SystemClass &system,int &sign);
  complex<double> evaluateRatio_check(SystemClass &system, int swap1, int swap2);

  double Psi_alpha_over_psi_FD(SystemClass &system);
  complex<double> evaluateRatio_noStore(SystemClass &system,int swap1, int swap2);
  


  void FillDet_check(SystemClass &system,SmartMatrix &myMat);

};

#endif
