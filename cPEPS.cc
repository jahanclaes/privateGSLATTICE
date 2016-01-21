#ifndef cPEPS_CC
#define cPEPS_CC

#include "cPEPS.h"
#include "SmartPfaffian.h"

void 
cPEPSClass::AllDerivs(SystemClass &system, Array<complex<double>,1> &derivs)
{
  AllDerivs(system,derivs,0,derivs.size());
}

void 
cPEPSClass::AllDerivs(SystemClass &system, Array<complex<double>,1>  &derivs,int start, int stop)
{
  complex<double> currentValue=evaluate(system);
  int** PhyCfg = new int* [xL];
  for(int i=0; i<xL; i++) PhyCfg[i] = new int [yL];
  // Specifically for Hubbard, change the following if using other models
  for(int i=0; i<xL*yL; i++)
    {
      if(system.x(i)==0 ) PhyCfg[i/yL][i%yL] = 0;
      if(system.x(i)==1 ) PhyCfg[i/yL][i%yL] = 1;
      if(system.x(i)==-1) PhyCfg[i/yL][i%yL] = 2;
      if(system.x(i)==2 ) PhyCfg[i/yL][i%yL] = 3;
   }
  ///////////////////////////////
  peps->diffPEPS(PhyCfg);
  ///////////////////////////////
  // all set to zero, for precaution
  for(int i=start; i<stop; i++) derivs(i)=0;
  ///////////////////////////////
  for(int ii=start; ii<stop; ii++)
  {
	  int i = ii-start;
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
	  if(PhyCfg[x][y]==pT) derivs(ii) = peps->dTN[x].T[y][pF](idx)/currentValue;
  }
  ///////////////////////////////
  for(int i=0; i<xL; i++) delete [] PhyCfg[i];
  delete [] PhyCfg;
  for(int i=start; i<stop; i++) 
    CheckDerivs(system,derivs,i);

  return;
}


void 
cPEPSClass::RealDerivs(SystemClass &system, Array<complex<double>,1> &derivs)
{
  assert(1==2);
}





double
cPEPSClass::CheckDerivs(SystemClass &system, Array<complex<double>,1>  &derivs,int derivInt)
{
  //  cerr<<"Testing "<<derivInt<<endl;
  double currParam=GetParam_real(derivInt);
//   for (double myStep=-0.01;myStep<0.01;myStep+=0.001){
//     double step_energy=0.0;
//     double countSteps=0.0;
//     SetParam_real(derivInt,currParam+myStep);
//     cerr<<myStep<<" "<<evaluate_noInverse(system)<<endl;
//   }

  SetParam_real(derivInt,currParam+0.001);
  double up=evaluate(system).real();
  SetParam_real(derivInt,currParam-0.001);
  double down=evaluate(system).real();

  SetParam_real(derivInt,currParam);
  double curr=evaluate(system).real();

  //uncomment next line if you want it to print
  //  cerr<<"Derivative of "<<derivInt<<" "<<derivs(derivInt)<<" "<<(up-down)/(0.002*curr)<<endl;
  derivs(derivInt)=(up-down)/(0.002*curr);
  //  cerr<<"Done Testing "<<derivInt<<endl;
  return  (up-down)/(0.002*curr);



}



complex<double> 
cPEPSClass::evaluate(SystemClass &system)
{
  int** PhyCfg = new int* [xL];
  for(int i=0; i<xL; i++) PhyCfg[i] = new int [yL];
  for(int i=0; i<xL*yL; i++)
    {
      if(system.x(i)==0 ) PhyCfg[i/yL][i%yL] = 0;
      if(system.x(i)==1 ) PhyCfg[i/yL][i%yL] = 1;
      if(system.x(i)==-1) PhyCfg[i/yL][i%yL] = 2;
      if(system.x(i)==2 ) PhyCfg[i/yL][i%yL] = 3;
    }
  complex<double> temp = peps->contractPEPS(PhyCfg);
  for(int i=0; i<xL; i++) delete [] PhyCfg[i];
  delete [] PhyCfg;
  return temp;
  //assert(1==2);
  
}

complex<double> 
cPEPSClass::evaluateRatio(SystemClass &system,int start, int stop, int spin)
{
  //IMPLEMENT ME!
  //  complex<double> toCheck=evaluateRatio_check(system,swap1,swap2);
  complex<double> temp;
  double ratio=1.0;
  ///for each particle you need to know all the correlators you've upset
  temp = evaluate(system);
  double newVal=temp.real();  //FIX ME !f(system.x);
  system.Move(stop,start,spin);
  temp = evaluate(system);
  double oldVal=temp.real();  //FIX ME f(system.x);
  system.Move(start,stop,spin);
  ratio*=(newVal/oldVal);
  return ratio;
  
}

//assumes x(swap1) != x(swap2)
//When swap1 and swap2 exchange, the term that corresponds to them
//doesn't change at all. 
//Called after the particles have been swapped!
complex<double> 
cPEPSClass::evaluateRatio(SystemClass &system, int swap1, int swap2)
{
  assert(1==2);
}


complex<double> 
cPEPSClass::logevaluate(SystemClass &system,int &sign)
{
  assert(1==2);
  complex<double> val = evaluate(system);
  sign= val.real() >0 ? 1: -1;
  return log(abs(val.real()));
}


complex<double> 
cPEPSClass::evaluateRatio_check(SystemClass &system, int swap1, int swap2)
{
  assert(1==2);
}



#endif
