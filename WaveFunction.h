#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H

#include "Blitz.h"
#include "SystemClass.h"

class WaveFunctionClass
{
public:
  string Name;
  int NumParams;
  bool NeedFrequentReset;
  virtual void RebuildParams();
  virtual void Init(SystemClass &system);
  virtual void AllDerivs(SystemClass &system, Array<complex<double>,1> &derivs,int start,int stop);
  virtual double GetParam(int i);
  virtual void SetParam(int i, double param);
  virtual double GetParam_real(int i);
  virtual void SetParam_real(int i, double param);

  virtual double GetParam_imag(int i);
  virtual void SetParam_imag(int i, double param);

  virtual void 
    CheckDerivs(SystemClass &system, 
		Array<complex<double>,1>  &derivs,int start, int stop);

  virtual void Copy(WaveFunctionClass* wf);
  virtual WaveFunctionClass* clone();

  virtual complex<double> logevaluate(SystemClass &system,int &sign)=0;
  virtual complex<double> evaluate(SystemClass &system)=0;
  virtual complex<double> evaluateRatio(SystemClass &system,int swap1, int swap2)=0;
  virtual complex<double> evaluateRatio_check(SystemClass &system, int swap1, int swap2)=0;
  virtual complex<double> evaluateRatio(SystemClass &system,int start, int stop, int spin);
  virtual void Swap(int i, int j)=0;
  virtual void Move(int site, int end_site, int spin);
  virtual void UpdateDets(SystemClass &system,int swap1, int swap2)=0;
  virtual void Reject(SystemClass &system,int swap1,int swap2);

  virtual void UpdateDets(SystemClass &system,int site, int end_site,int spin);
  virtual void Reject(SystemClass &system,int site,int end_site,int spin);


  virtual void SetParams(double delta3,SystemClass &system);
};

#endif
