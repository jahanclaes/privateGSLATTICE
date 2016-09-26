#include "CPS.h"
#include "SmartPfaffian.h"

void 
CPSClass::AllDerivs(SystemClass &system, Array<complex<double>,1> &derivs)
{
  AllDerivs(system,derivs,0,derivs.size());
}

void 
CPSClass::AllDerivs(SystemClass &system, Array<complex<double>,1>  &derivs,int start, int stop)
{
  for (int i=start;i<stop;i++)
    derivs(i)=0.0;

  for (int i=0;i<PF.NumCorrelators;i++){
    int theBin=PF.corr2Bin(system.x,i);
    assert(start+PF.binLoc[i][theBin]<derivs.size());
    derivs(start+PF.binLoc[i][theBin])+=1.0/(PF.corr2Val(system.x,i)); // f0[PairingFunction.FindBin(r)].real());

  }

  return;
}


void 
CPSClass::RealDerivs(SystemClass &system, Array<complex<double>,1> &derivs)
{
  int start=0;
  int stop=derivs.size();
  for (int i=start;i<stop;i++)
    derivs(i)=0.0;
  for (int j=0;j<PF.NumCorrelators;j++){
    int theBin=PF.corr2Bin(system.x,j);
    double total=1.0;
    for (int i=0;i<PF.NumCorrelators;i++){
      if (i!=theBin)
	total*=PF.corr2Val(system.x,i).real();
    }
    derivs(start+PF.binLoc[j][theBin])+=total;
  }

}





void 
CPSClass::CheckDerivs(SystemClass &system, Array<complex<double>,1>  &derivs,int start, int stop)
{
  Array<double,1> check_derivs(derivs.size());
  for (int i=start;i<stop;i++)
    check_derivs(i)=0.0;
  double eps=1e-10;
  for (int i=0;i<PF.NumCorrelators;i++){
    int theBin=PF.corr2Bin(system.x,i);
    PF.f0[i][theBin]+=eps;
    complex<double> up=evaluate(system);
    PF.f0[i][theBin]-=2*eps;
    complex<double> down=evaluate(system);
    PF.f0[i][theBin]+=eps;
    complex<double> thestart=evaluate(system);
    assert(start+PF.binLoc[i][theBin]<stop);
    check_derivs(start+PF.binLoc[i][theBin])=(up.real()-down.real())/(2*eps*thestart.real());

    //    derivs(start+PF.binLoc[i][theBin])+=1.0/(PF.corr2Val(system.x,i)); // f0[PairingFunction.FindBin(r)].real());
  }
  for (int i=start;i<stop;i++)
    cerr<<"DERIV CHECK: "<<derivs(i).real()<<" "<<check_derivs(i)<<endl;

  return;
}



complex<double> 
CPSClass::evaluate(SystemClass &system)
{
  double total=1.0;
  for (int i=0;i<PF.NumCorrelators;i++){
    total*=PF.corr2Val(system.x,i).real();
  }
  return total;
}


//assumes x(swap1) != x(swap2)
//When swap1 and swap2 exchange, the term that corresponds to them
//doesn't change at all. 
//Called after the particles have been swapped!
complex<double> 
CPSClass::evaluateRatio(SystemClass &system, int swap1, int swap2)
{
  //  complex<double> toCheck=evaluateRatio_check(system,swap1,swap2);
  //  cerr<<"ENDTERING"<<endl;
  double ratio=1.0;
  ///for each particle you need to know all the correlators you've upset
  set<int> correlatorsChanged;
  correlatorsChanged.insert(PF.correlatorsForSite[swap1].begin(),PF.correlatorsForSite[swap1].end());
  correlatorsChanged.insert(PF.correlatorsForSite[swap2].begin(),PF.correlatorsForSite[swap2].end());
  //  cerr<<"PRE FOR LOOP"<<endl;
  for (set<int>::iterator corrIter = correlatorsChanged.begin();corrIter!=correlatorsChanged.end();corrIter++){
    int corr=*corrIter;
    double newVal=PF.corr2Val(system.x,corr).real();
    //    cerr<<"swap1 swap2 "<<swap1<<" "<<swap2<<" "<<corr<<endl;
    //    for (int k=0;k<PF.myCorrs[corr].size();k++)
    //      cerr<<PF.myCorrs[corr][k]<<endl;
    //    cerr<<"newVal is "<<newVal<<endl;
    system.Swap(swap1,swap2);
    double oldVal=PF.corr2Val(system.x,corr).real();
    system.Swap(swap1,swap2);
    ratio*=(newVal/oldVal);
  }
  //  cerr<<"POST FOR LOOP"<<endl;
  //  cerr<<"CHECKING: "<<ratio<<" "<<toCheck<<endl;
  return ratio;
}


// Called AFTER spin flip
// No side effects
complex<double>
CPSClass::evaluateRatioFlip(SystemClass &system, int swap1)
{
  double ratio=1.0;
  set<int> correlatorsChanged;
  correlatorsChanged.insert(PF.correlatorsForSite[swap1].begin(),PF.correlatorsForSite[swap1].end());
  for (set<int>::iterator corrIter = correlatorsChanged.begin();corrIter!=correlatorsChanged.end();corrIter++){
    int corr=*corrIter;
    double newVal=PF.corr2Val(system.x,corr).real();
    system.Flip(swap1);
    double oldVal=PF.corr2Val(system.x,corr).real();
    system.Flip(swap1);
    ratio*=(newVal/oldVal);
  }
  return ratio;
}




complex<double> 
CPSClass::logevaluate(SystemClass &system,int &sign)
{
  complex<double> val = evaluate(system);
  sign= val.real() >0 ? 1: -1;
  return log(abs(val.real()));
}


complex<double> 
CPSClass::evaluateRatio_check(SystemClass &system, int swap1, int swap2)
{
  system.Swap(swap1,swap2);
  int signa;
  complex<double> before=logevaluate(system,signa);
  system.Swap(swap1,swap2);
  int signb;
  complex<double> after=logevaluate(system,signb);
  return exp(after.real()-before.real())*signa*signb;
}

