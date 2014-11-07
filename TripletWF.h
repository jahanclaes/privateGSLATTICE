#ifndef TRIPLET_WAVE_FUNCTION_H
#define TRIPLET_WAVE_FUNCTION_H


#include <Eigen/Dense>
#include "Blitz.h" 
#include "SystemClass.h"
#include "SmartPfaffian.h"
#include "WaveFunction.h"
#include "PairingFunctionAllBin.h"


class TripletWF : public WaveFunctionClass
{
 public:
  SmartPfaffian UpM;
  PairingFunctionAllBin &PF;
  void ReadPairingFunction();
  TripletWF(PairingFunctionAllBin &t_PF): PF(t_PF)
  {

  }
  
  int SpinIndex;
  int UpIndex;
  Eigen::VectorXd NewUpCol;
  
  PairingFunctionAllBin pf_test;

  void FillUp(SystemClass &system);
  void FillUp(SystemClass &system,SmartPfaffian &UpM);
  void Init(SystemClass &system); //done
  complex<double> evaluate(SystemClass &system); //done
  
  complex<double> logevaluate(SystemClass &system,int &sign);
  complex<double> evaluateRatio(SystemClass &system,int swap1, int swap2);
  complex<double> evaluateRatio_check(SystemClass &system, int swap1, int swap2);
  
  void Swap(int i, int j);
  void UpdateDets(SystemClass &system,int swap1, int swap2);
  
  void Reject(SystemClass &system,int swap1,int swap2);
  complex<double> Phi(int i, int j);
  
  double Sign(SystemClass &system);
  
  void CheckPfaffian(SystemClass &system);
  

  void  FillDet_check(SystemClass &system,SmartPfaffian &myMat);
  
  void  CheckDerivs(SystemClass &system, 
		    Array<complex<double>,1>  &derivs,int start, int stop);
  
  
  void AllDerivs(SystemClass &system, Array<complex<double>,1> &derivs,
		 int start,int stop);
  
  /*   virtual double GetParam(int i); */
  /*   virtual void SetParam(int i, double param); */
  
  double GetParam_real(int i);
  double GetParam_imag(int i); 
  void SetParam_real(int i, double param); 
  void SetParam_imag(int i, double param); 






/*   virtual void SetParams(double delta3,SystemClass &system); */
};

#endif
