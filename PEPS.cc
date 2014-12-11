#include "PEPS.h"
#include "SmartPfaffian.h"

void 
PEPSClass::AllDerivs(SystemClass &system, Array<complex<double>,1> &derivs)
{
  AllDerivs(system,derivs,0,derivs.size());
}

void 
PEPSClass::AllDerivs(SystemClass &system, Array<complex<double>,1>  &derivs,int start, int stop)
{
  for (int i=start;i<stop;i++)
    derivs(i)=0.0;

  //IMPLEMENT ME!!
  return;
}


void 
PEPSClass::RealDerivs(SystemClass &system, Array<complex<double>,1> &derivs)
{
  assert(1==2);
}





void 
PEPSClass::CheckDerivs(SystemClass &system, Array<complex<double>,1>  &derivs,int start, int stop)
{
  assert(1==2);
  return;
}



complex<double> 
PEPSClass::evaluate(SystemClass &system)
{
  //MIGHT IMPLEMENT
  assert(1==2);
  
}

complex<double> 
PEPSClass::evaluateRatio(SystemClass &system,int start, int stop, int spin)
{
  //IMLEMENT ME!
  //  complex<double> toCheck=evaluateRatio_check(system,swap1,swap2);
  double ratio=1.0;
  ///for each particle you need to know all the correlators you've upset
  double newVal=0.0;  //FIX ME !f(system.x);
  system.Move(start,stop,spin);
  double oldVal=0.0;  //FIX ME f(system.x);
  system.Move(start,stop,spin);
  ratio*=(newVal/oldVal);
  return ratio;
  
}

//assumes x(swap1) != x(swap2)
//When swap1 and swap2 exchange, the term that corresponds to them
//doesn't change at all. 
//Called after the particles have been swapped!
complex<double> 
PEPSClass::evaluateRatio(SystemClass &system, int swap1, int swap2)
{
  assert(1==2);
}


complex<double> 
PEPSClass::logevaluate(SystemClass &system,int &sign)
{
  assert(1==2);
  complex<double> val = evaluate(system);
  sign= val.real() >0 ? 1: -1;
  return log(abs(val.real()));
}


complex<double> 
PEPSClass::evaluateRatio_check(SystemClass &system, int swap1, int swap2)
{
  assert(1==2);
}

