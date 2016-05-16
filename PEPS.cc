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
  complex<double> currentValue=evaluate(system);
  int** PhyCfg = new int* [W_];
  for(int i=0; i<W_; i++) PhyCfg[i] = new int [L_];
  for(int i=0; i<L_*W_; i++)
    {
      if(system.x(i)==0 ) PhyCfg[i/L_][i%L_] = 0;
      if(system.x(i)==1 ) PhyCfg[i/L_][i%L_] = 1;
      if(system.x(i)==-1) PhyCfg[i/L_][i%L_] = 2;
      if(system.x(i)==2 ) PhyCfg[i/L_][i%L_] = 3;
    }

  ///////////////////////////////
  peps->diff(PhyCfg);
  ///////////////////////////////
  for(int i=start; i<stop; i++) derivs(i)=0;
  ///////////////////////////////
  int N_edge = 4*(2*D_*D_+D_*D_*D_*(L_-2));
  int N_mid  = 4*(2*D_*D_*D_+D_*D_*D_*D_*(L_-2));
  int r_, c_, phy_, s_, l_;
  for (int i=start;i<stop;i++)
    {
      if( i<0 || i>=(2*N_edge+(W_-2)*N_mid) )
      {
	cout<<"GetParam_real(int i): i out of bound!"<<endl;
	exit(0);
      }
      else if(i<N_edge)
	{
	  r_ = 1;
	  if(i<4*D_*D_)
	    {
	      c_ = 1;
	      phy_ = i/(D_*D_);
	      s_ = (i%(D_*D_))/D_;
	      l_ = (i%(D_*D_))%D_;
	    }
	  else if( i>=(N_edge-4*D_*D_) )
	    {
	      c_ = L_ - 1;
	      phy_ = (i+4*D_*D_-N_edge)/(D_*D_);
	      s_ = ((i+4*D_*D_-N_edge)%(D_*D_))/D_;
	      l_ = ((i+4*D_*D_-N_edge)%(D_*D_))%D_;
	    }
	  else
	    {
	      c_ = (i-4*D_*D_)/(4*D_*D_*D_) + 1;
	      phy_ = (i-(c_-1)*4*D_*D_*D_-4*D_*D_)/(D_*D_*D_);
	      s_ = ((i-(c_-1)*4*D_*D_*D_-4*D_*D_)%(D_*D_*D_))/(D_*D_);
	      l_ = ((i-(c_-1)*4*D_*D_*D_-4*D_*D_)%(D_*D_*D_))%(D_*D_);
	    }
	}
      else if( i>=N_edge && i<(NumParams-N_edge) )
	{
	  r_ = (i-N_edge)/N_mid + 1;
	  int ii = (i-N_edge)%N_mid;
	  if(ii<4*D_*D_*D_)
	    {
	      c_ = 1;
	      phy_ = ii/(D_*D_*D_);
	      s_ = (ii%(D_*D_*D_))/D_;
	      l_ = (ii%(D_*D_*D_))%D_;
	    }
	  else if( ii>=(N_mid-4*D_*D_*D_) )
	    {
	      c_ = L_ - 1;
	      phy_ = (ii+4*D_*D_*D_-N_mid)/(D_*D_*D_);
	      s_ = ((ii+4*D_*D_*D_-N_mid)%(D_*D_*D_))/D_;
	      l_ = ((ii+4*D_*D_*D_-N_mid)%(D_*D_*D_))%D_;
	    }
	  else
	    {
	      c_ = (ii-4*D_*D_*D_)/(4*D_*D_*D_*D_) + 1;
	      phy_ = (ii-(c_-1)*4*D_*D_*D_*D_-4*D_*D_*D_)/(D_*D_*D_*D_);
	      s_ = ((ii-(c_-1)*4*D_*D_*D_*D_-4*D_*D_*D_)%(D_*D_*D_*D_))/(D_*D_);
	      l_ = ((ii-(c_-1)*4*D_*D_*D_*D_-4*D_*D_*D_)%(D_*D_*D_*D_))%(D_*D_);
	    }
	}
      else if( i>=(NumParams-N_edge) )
	{
	  r_ = W_ - 1;
	  int ii = i+N_edge-NumParams;
	  if(ii<4*D_*D_)
	    {
	      c_ = 1;
	      phy_ = ii/(D_*D_);
	      s_ = (ii%(D_*D_))/D_;
	      l_ = (ii%(D_*D_))%D_;
	    }
	  else if( ii>=(N_edge-4*D_*D_) )
	    {
	      c_ = L_ - 1;
	      phy_ = (ii+4*D_*D_-N_edge)/(D_*D_);
	      s_ = ((ii+4*D_*D_-N_edge)%(D_*D_))/D_;
	      l_ = ((ii+4*D_*D_-N_edge)%(D_*D_))%D_;
	    }
	  else
	    {
	      c_ = (ii-4*D_*D_)/(4*D_*D_*D_) + 1;
	      phy_ = (ii-(c_-1)*4*D_*D_*D_-4*D_*D_)/(D_*D_*D_);
	      s_ = ((ii-(c_-1)*4*D_*D_*D_-4*D_*D_)%(D_*D_*D_))/(D_*D_);
	      l_ = ((ii-(c_-1)*4*D_*D_*D_-4*D_*D_)%(D_*D_*D_))%(D_*D_);
	    }
	}

      if(PhyCfg[r_][c_]==phy_) derivs(i) = peps->DiffT_[r_][c_][s_](l_)/currentValue;
    }
  ///////////////////////////////
  for(int i=0; i<W_; i++) delete [] PhyCfg[i];
  delete [] PhyCfg;
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
  int** PhyCfg = new int* [W_];
  for(int i=0; i<W_; i++) PhyCfg[i] = new int [L_];
  for(int i=0; i<L_*W_; i++)
    {
      if(system.x(i)==0 ) PhyCfg[i/L_][i%L_] = 0;
      if(system.x(i)==1 ) PhyCfg[i/L_][i%L_] = 1;
      if(system.x(i)==-1) PhyCfg[i/L_][i%L_] = 2;
      if(system.x(i)==2 ) PhyCfg[i/L_][i%L_] = 3;
    }
  complex<double> temp = peps->contract(PhyCfg,W_/2);
  for(int i=0; i<W_; i++) delete [] PhyCfg[i];
  delete [] PhyCfg;
  return temp;
  //assert(1==2);
  
}

complex<double> 
PEPSClass::evaluateRatio(SystemClass &system,int start, int stop, int spin)
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

