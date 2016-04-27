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

void WaveFunctionClass::Move(int site, int end_site, int spin)
{
  assert(1==2);
}


void WaveFunctionClass::UpdateDets(SystemClass &system,int site, int end_site,int spin)
{
  assert(1==2);
}

void WaveFunctionClass::Reject(SystemClass &system,int site,int end_site,int spin)
{
  assert(1==2);
}

complex<double> WaveFunctionClass::evaluateRatio(SystemClass &system,int start, int stop, int spin)
{
  assert(1==2);

}

void WaveFunctionClass::RebuildParams()
{

}
WaveFunctionClass* WaveFunctionClass::clone()
{
  assert(1==2);
}

void WaveFunctionClass::MakeProductState(vector<int> &myState)
{
  assert(1==2);
}

void WaveFunctionClass::MakeUniformState()
{
  assert(1==2);
}
