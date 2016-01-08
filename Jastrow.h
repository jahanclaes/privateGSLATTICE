#ifndef JASTROW_H_
#define JASTROW_H_
#include "SystemClass.h"
#include "WaveFunction.h" 

class JastrowClass : public WaveFunctionClass
{
  public:

  Array<int,1> reduced_index;
  Array<double,1> vij;
  Array<double,1> Ti;
 
  void Init(SystemClass &system);

  complex<double> logevaluate(SystemClass &system,int &sign);
  complex<double> evaluate(SystemClass &system);
  complex<double> evaluateRatio(SystemClass &system,int swap1, int swap2);
  complex<double> evaluateRatio_check(SystemClass &system, int swap1, int swap2);
  complex<double> evaluateRatio(SystemClass &system,int start, int stop, int spin);
  complex<double> evaluateRatio_check(SystemClass &system, int start,int stop,int spin);

  void Swap(int i, int j);
  void Move(int site, int end_site, int s);
  void Reject(SystemClass &system,int site,int end_site,int spin);
 
  void AllDerivs(SystemClass &system, Array<complex<double>,1> &derivs,int start,int stop); //compute the derivative
  void CheckDerivs(SystemClass &system, Array<complex<double>,1> &derivs,int start, int stop);
  int  GetReducedIndex(SystemClass &system,int i,int j);
 
  double GetParam(int i);
  void SetParam(int i, double param);

  double GetParam_real(int i);
  void SetParam_real(int i, double param);
  double GetParam_imag(int i);
  void SetParam_imag(int i, double param);

  void InitTi(SystemClass &system);

  void GetParameters();

  void UpdateDets(SystemClass &system,int site, int end_site,int spin);
  void UpdateDets(SystemClass &system,int swap1, int swap2);
  complex<double> Sign(SystemClass &system); 
};

#endif
