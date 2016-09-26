#include "WaveFunction.h"

void WaveFunctionClass::Copy(WaveFunctionClass* wf)
{
  assert(1==2);
}

double WaveFunctionClass::GetParam(int i)
{
  assert(1==2);
}

void WaveFunctionClass::SetParam(int i, double param)
{
  assert(1==2);

}

void WaveFunctionClass::CheckDerivs(SystemClass &system, 
				    Array<complex<double>,1>  &derivs,int start, int stop)
{
  assert(1==2);
}



double WaveFunctionClass::GetParam_real(int i)
{
  assert(1==2);
}

void WaveFunctionClass::SetParam_real(int i, double param)
{
  assert(1==2);

}


double WaveFunctionClass::GetParam_imag(int i)
{
  assert(1==2);
}

void WaveFunctionClass::SetParam_imag(int i, double param)
{
  assert(1==2);

}

void WaveFunctionClass::Init(SystemClass &system)

{
  assert(1==2);
}

void WaveFunctionClass::AllDerivs(SystemClass &system, Array<complex<double>,1> &derivs,int start,int stop)
{
  assert(1==2);

}
void WaveFunctionClass::Reject(SystemClass &system,int swap1,int swap2)
{

}

void 
WaveFunctionClass::SetParams(double delta3,SystemClass &system)
{

}
complex<double> 
WaveFunctionClass::evaluateRatioFlip(SystemClass &system, int swap1)
{

}
