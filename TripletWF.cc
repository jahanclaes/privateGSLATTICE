#include "TripletWF.h"
#include <fstream>
using namespace std;
void TripletWF::ReadPairingFunction()
{
  int bin=0;
  ifstream infile;
  infile.open("PairingFunction.dat");
  assert(infile);
  while (!infile.eof()){
    int i;      int j; double tr; double ti;
    assert(!infile.eof());
    infile>>i;
    if (!infile.eof()){
      assert(!infile.eof());
      infile>>j;
      bin=PF.FindBin(i,j);
      assert(bin<PF.f0.size());
      assert(!infile.eof());
      infile>>PF.f0[bin].real();
      PF.f0[bin].real()=PF.f0[bin].real(); ///100.0;
      assert(!infile.eof());
      infile>>PF.f0[bin].imag();
      PF.f0[bin].imag()=PF.f0[bin].imag(); ///100.0;
    }
    cerr<<i<<" "<<j<<" "<<PF.f0[bin]<<" "<<bin<<endl;
  }
  infile.close();
  cerr<<"READ IN PAIRING FUNCTIONS"<<endl;
}

//Pairing Function should be initialized somewhere!!
void TripletWF::Init(SystemClass &system)
{
  Name="TripletWF";
  NeedFrequentReset=true;

  PF.Init(system);

  pf_test.Init(system);
  pf_test.GetEigs(system);
  pf_test.gutz_build(system);

  UpM.Init(system.x.size()/2);
  NewUpCol.resize(system.x.size()/2);
  NumParams=PF.f0.size();

  
}

double 
TripletWF::Sign(SystemClass &system)
{
  vector<int> detPos;
  for (int i=0;i<UpM.DetPos.size();i++)
    if (system.x(i)==SpinIndex)
      detPos.push_back(UpM.DetPos[i]);
  int numSwap=0;

  for (int i=0;i<detPos.size();i++){
    while (detPos[i]!=i){
      numSwap++;
      int temp=detPos[i];
      swap(detPos[i],detPos[temp]);

    }
    
  }
  if (numSwap %2==0)
    return 1;
  else 
    return -1;
}


void TripletWF::CheckPfaffian(SystemClass &system)
{
 SmartPfaffian UpM_check;
 UpM_check.Init(system.x.size()/2);
 FillUp(system,UpM_check);
 cerr<<"PFaffian should be "<<UpM_check.Pfaffian()<<" "<<Sign(system)<<endl;
}

void TripletWF::FillUp(SystemClass &system,SmartPfaffian &UpM)
{


  UpM.M=Eigen::MatrixXd::Zero(system.x.size()/2,system.x.size()/2);
  int first=-1;
  for (int i=0;i<system.x.size();i++){
    if (system.x(i)==SpinIndex){
      first++;
      int second=first;
      for (int j=i+1;j<system.x.size();j++){
	if (system.x(j)==SpinIndex){
	  second++;
	  UpM.M(first,second)=PF.Phi(i,j).real();
	  UpM.M(second,first)=-UpM.M(first,second);
	  UpM.DetPos[i]=first;
	  UpM.DetPos[j]=second;
	}
	
      }
      
    }
  }
  
}
    

void TripletWF::FillUp(SystemClass &system)
{
  FillUp(system,UpM);
}



complex<double> TripletWF::evaluate(SystemClass &system)
{

  FillUp(system);
  complex<double> myAns=UpM.Pfaffian();
  UpM.CalcAndStoreInverse();

  return myAns;
}



complex<double> TripletWF::logevaluate(SystemClass &system,int &sign)
{
  assert(1==2);

}
//need a checkM
//need to set up Pairing 

complex<double> TripletWF::Phi(int i,int j)
{

  //  cerr<<"CHECKING: "<<pf_test.Phi(i,j)<<" "<<PF.Phi(i,j)<<" "<<PF.Phi(i,j)/pf_test.Phi(i,j)<<endl;
  return PF.Phi(i,j);
}


void TripletWF::SetParam_real(int i,double param)
{
  PF.f0[i].real()=param;
}

void TripletWF::SetParam_imag(int i,double param)
{
  PF.f0[i].imag()=param;
}


double TripletWF::GetParam_real(int i)
{
  return PF.f0[i].real();
}



double TripletWF::GetParam_imag(int i)
{
  return PF.f0[i].imag();
}


void 
TripletWF::FillDet_check(SystemClass &system,SmartPfaffian &myMat)
{
  cerr<<"The pfaffian I am checking is "<<endl<<myMat.M<<endl;
  int upDet=-1;
  int downDet=-1;
  bool ok=true;
  for (int i=0;i<system.x.size();i++){
    if (system.x(i)==SpinIndex){
      upDet++;
      downDet=-1;
      for (int j=0;j<system.x.size();j++){
	if (system.x(j)==SpinIndex){
	  downDet++;
	  
	  if ( (i==j)  && (abs(myMat.M(myMat.DetPos[i],myMat.DetPos[j]))>1e-10) ){
	    cerr<<"zero broke: "<<i<<" "<<j<<" "<<myMat.DetPos[i]<<" "<<myMat.DetPos[j]<<" "<<myMat.M(myMat.DetPos[i],myMat.DetPos[j])<<" "<<Phi(i,j)<<" "<<endl;
	    ok=false;
	  }
	  
	  if  ( (i > j) && (!(myMat.M(myMat.DetPos[i],myMat.DetPos[j])==-Phi(i,j)))){
	    cerr<<"pos broke: "<<i<<" "<<j<<" "<<myMat.DetPos[i]<<" "<<myMat.DetPos[j]<<" "<<myMat.M(myMat.DetPos[i],myMat.DetPos[j])<<" "<<-Phi(i,j)<<" "<<endl;
	    ok=false;
	  }

	  
	  if  ( (i < j) && (!(myMat.M(myMat.DetPos[i],myMat.DetPos[j])==Phi(i,j)))){
	    cerr<<"neg broke: "<<i<<" "<<j<<" "<<myMat.DetPos[i]<<" "<<myMat.DetPos[j]<<" "<<myMat.M(myMat.DetPos[i],myMat.DetPos[j])<<" "<<Phi(i,j)<<" "<<endl;
	    ok=false;
	  }


	}

      }
    }
  }
  for (int i=0;i<system.x.size();i++)
    cerr<<system.x(i)<<" ";
  cerr<<endl;
  assert(ok);
}



complex<double> TripletWF::evaluateRatio(SystemClass &system,int swap1, int swap2)
{

//assumes swap1 and swap2 are of different spin
  int count=0;
  for (int i=min(swap1,swap2)+1;i<max(swap1,swap2);i++){
    if (system.x(i)==SpinIndex){
      count++;
    }
  }
  double mySign = count % 2==0 ? 1: -1 ;
  ///let's define spin up as swap1
  if (system.x(swap1)!=SpinIndex){
    swap(swap1,swap2);
  }


  //evaluate the new spin up column
  UpIndex=UpM.DetPos[swap1];


  //swap1 has been set to be the spin up value
  //loop over the other spin up particles
  for (int j=0;j<system.x.size();j++){
    if (system.x(j)==SpinIndex && swap1!=j){
      //      if (UpM.DetPos[j]<UpIndex)
      if (j<swap1)
	NewUpCol(UpM.DetPos[j])=Phi(swap1,j).real()*-1; //-1 * 
      else 
	NewUpCol(UpM.DetPos[j])=Phi(swap1,j).real();
    }
    else if (system.x(j)==SpinIndex && swap1==j)
      NewUpCol(UpM.DetPos[j])=0.0;
  }
  complex<double> ratio=UpM.RowColRatio(UpIndex,NewUpCol)*(mySign*-1);
  
// complex<double> check_ratio = evaluateRatio_check(system,swap1,swap2);
//     complex<double> diff=ratio-check_ratio;
//      cerr<<"Diff: "<<diff<<" "<<ratio<<" "<<check_ratio<<" "<<mySign<<" "<<Sign(system)<<endl;
  //HACK!      assert((diff*conj(diff)).real()<1e-10);
  //  cerr<<"test ratio is "<<test_ratio<<" "<<endl; //evaluateRatio_check(system,swap1,swap2)<<endl;
  //  rebuild=false;
  return ratio;
}
  

complex<double> TripletWF::evaluateRatio_check(SystemClass &system, int swap1, int swap2)
{
  SmartPfaffian UpM_check;
  //  cerr<<"Swapping from "<<swap1<<" to "<<swap2<<endl;
  UpM_check.Init(system.x.size()/2);

  system.Swap(swap1,swap2);
  FillUp(system,UpM_check);
  Eigen::MatrixXd oldM=UpM_check.M;
  complex<double> pre=UpM_check.Pfaffian();
  system.Swap(swap1,swap2);
  FillUp(system,UpM_check);
  complex<double> post=UpM_check.Pfaffian();
  //  cerr<<"Diff of matrices: "<<endl;
  //  cerr<<oldM-UpM_check.M<<endl;
  //  cerr<<endl<<oldM<<endl;
  //  cerr<<endl<<UpM_check.M<<endl;
  cerr<<"Pre and Post: "<<pre<<" "<<post<<endl;

  return post/pre;


}


void TripletWF::Swap(int i, int j)
{
  swap(UpM.DetPos[i],UpM.DetPos[j]);
}


void TripletWF::UpdateDets(SystemClass &system,int swap1, int swap2)
{


  ///let's define spin up as swap1
  if (system.x(swap1)!=SpinIndex){
    swap(swap1,swap2);
  }


  UpM.InverseUpdate(UpIndex,NewUpCol);

  //  CheckPfaffian(system);
  //  FillDet_check(system,UpM);
}


void TripletWF::Reject(SystemClass &system, int swap1, int swap2)
{
  //do nothing for now
  

}

void 
TripletWF::AllDerivs(SystemClass &system, Array<complex<double>,1> &derivs,
		     int start,int stop)
{
  for (int i=start;i<stop;i++)
    derivs(i)=0.0;
  for (int i=0;i<system.x.size();i++){
    for (int j=0;j<system.x.size();j++){
      if (system.x(i)==SpinIndex && system.x(j)==SpinIndex){
	int bin=start+PF.PhiMap(i,j);
	if (i>j)
	  derivs(bin)+=UpM.MInverse(UpM.DetPos[i],UpM.DetPos[j])/2.0;
	else 
	  derivs(bin)-=UpM.MInverse(UpM.DetPos[i],UpM.DetPos[j])/2.0;
      }
      else {
	
      }
    }
  }



}


void 
TripletWF::CheckDerivs(SystemClass &system, 
		       Array<complex<double>,1>  &derivs,int start, int stop)
{
  Array<double,1> check_derivs(derivs.size());
  for (int i=start;i<stop;i++)
    check_derivs(i)=0.0;
  double eps=1e-10;
  
  for (int i=0;i<PF.f0.size();i++){
    PF.f0[i]+=eps;
    complex<double> up=evaluate(system);
    PF.f0[i]-=2*eps;
    complex<double> down=evaluate(system);
    PF.f0[i]+=eps;
    complex<double> thestart=evaluate(system);
    assert(start+i<stop);
    check_derivs(start+i)=(up.real()-down.real())/(2*eps*thestart.real());
  }
  for (int i=start;i<stop;i++)
    cerr<<"DERIV CHECK: "<<derivs(i).real()<<" "<<check_derivs(i)<<endl;
  return;
}
