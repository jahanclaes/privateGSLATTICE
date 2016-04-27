#ifndef  OPTIMIZE_BOTH_CLASS_H
#define  OPTIMIZE_BOTH_CLASS_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include "HuseElser_pair2.h"
#include "RVB_fast.h"

#include "CSP.h"

using namespace blitz;
using namespace std;


class OptimizeBothClass
{
public:
  RandomClass Random;
  SystemClass System;
  int VMC_equilSweeps;
  int VMC_SampleSweeps;
  int opt_equilSweeps;
  int opt_SampleSweeps;

  int sign;
  double J2;


  HuseElserPair2Class Psi;
  //  HuseElserPair2Class Psi_HE_2;
  //  HuseElserPair2Class Psi_HE_3;
  //  HuseElserPair2Class Psi_HE_4;
  //  HuseElserPair2Class Psi_HE_5;
  //RVBFastPsiClass Psi_HE_2;
  RVBFastPsiClass Psip;
  CSPClass Psi_CSP;

  void Copy(OptimizeBothClass &b)
  {
    sign=b.sign;
    System=b.System;
    Psi=b.Psi;
    Psip=b.Psip;
    Psi_CSP=b.Psi_CSP;
    //    Psi_HE_2=b.Psi_HE_2;
    //    Psi_HE_3=b.Psi_HE_3;
    //    Psi_HE_4=b.Psi_HE_4;
    //    Psi_HE_5=b.Psi_HE_5;
    wf_list.clear();
    //    wf_list.push_back(&Psi);
    wf_list.push_back(&Psip);
    //a    wf_list.push_back(&Psi_CSP);
    //    wf_list.push_back(&Psi_HE_2);
    //    wf_list.push_back(&Psi_HE_3);
    //    wf_list.push_back(&Psi_HE_4);
    //    wf_list.push_back(&Psi_HE_5);
  }


  list<WaveFunctionClass*> wf_list;


  NNHeisenberg H;
  NNNHeisenberg Hp;
  int makeNeg;
  int OptimizeCalled;
  int NumberOfParams;
  void Init()
  {
    VMC_equilSweeps=1000;
    //    VMC_SampleSweeps=10000;
    VMC_SampleSweeps=10000;
    opt_equilSweeps=1000;
   opt_SampleSweeps=100000;
    //    opt_SampleSweeps=10000;
    //    opt_equilSweeps=2;
    //    opt_SampleSweeps=2;


    sign=1;
    System.Init();
    System.Stagger();

    wf_list.clear();
    //    wf_list.push_back(&Psi);
    wf_list.push_back(&Psip);
    //a    wf_list.push_back(&Psi_CSP);
    //    wf_list.push_back(&Psi_HE_2);

    //    wf_list.push_back(&Psi_HE_3);
    //    wf_list.push_back(&Psi_HE_4);
    //    wf_list.push_back(&Psi_HE_5);
    for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++)
      (*wf_iter)->Init(System); 

    //HACK!
    Psi.Init(System);
    //    cerr<<"READING CSP"<<endl;
    //    Psi_CSP.Init(System);
    //    cerr<<"READING CSP DONE"<<endl;
    ////    Psi.PairingFunction.f0[makeNeg]=-1;
    /*     ifstream infile; */
/*     infile.open("toConnect.dat"); */
/*     assert(infile); */
    
/*     while (!infile.eof()){ */
/*       int i1; */
/*       int j1; */
/*       infile>>i1; */
/*       infile>>j1; */
/*       Psip.PairingFunction.f0[Psip.PairingFunction.PhiMap(i1,j1)].real()=0.0; //10.0; //Psip.PairingFunction.f0[Psip.PairingFunction.PhiMap(i1,j1)].real()*100; */
/*     } */
/*     infile.close(); */

    NumberOfParams=0;
    for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      cerr<<"THE NUMBER OF PARAMS IS "<<(*wf_iter)->NumParams<<endl;
      NumberOfParams+=(*wf_iter)->NumParams;
    }
    ParamsOld.resize(NumberOfParams);
    InitDeriv();
    OptimizeCalled=0;
    H.Init(System);
    Hp.Init(System,J2);
    if (1==2){
     for (int i=0;i<System.x.size();i++) 
       for (int j=i+1;j<System.x.size();j++){ 
 	set<int> theSet; 
 	std::set_intersection( Psi_CSP.PF.correlatorsForSite[i].begin(),  Psi_CSP.PF.correlatorsForSite[i].end(), 
 			       Psi_CSP.PF.correlatorsForSite[j].begin(),  Psi_CSP.PF.correlatorsForSite[j].end(),  
 			       std::inserter( theSet, theSet.begin() ) ); 
 	int theBin=Psi.PairingFunction.FindBin(i,j); 
 	int CSP_bin=(*(theSet.begin())); 
 	if (Psi.PairingFunction.f0[theBin].real()==0 || Psi.PairingFunction_diff.f0[theBin].real()==0) 
 	  assert(1==2); 
 	Psi_CSP.PF.f0[CSP_bin][0]=Psi.PairingFunction.f0[theBin]; 
 	Psi_CSP.PF.f0[CSP_bin][3]=Psi.PairingFunction.f0[theBin]; 
 	Psi_CSP.PF.f0[CSP_bin][1]=Psi.PairingFunction_diff.f0[theBin]; 
 	Psi_CSP.PF.f0[CSP_bin][2]=Psi.PairingFunction_diff.f0[theBin]; 
       } 
     for (int i=0;i<Psi_CSP.PF.f0.size();i++) 
       if (Psi_CSP.PF.f0[i][0].real()==0 || Psi_CSP.PF.f0[i][1].real()==0) 
 	assert(1==2); 
    }
  }


  double Sweep()
  {
    int numAccepted=0;
    int numAttempted=0;
    //HACK!

    for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      if ((*wf_iter)->NeedFrequentReset)
	(*wf_iter)->evaluate(System);
    }

      //     ((RVBFastPsiClass*)(&Psip))->evaluate(System);
     //     ((RVBFastPsiClass*)(&Psi_HE_2))->evaluate(System);
     //END HACK!

    for (int step=0;step<System.x.size();step++){
      int swap1=step;
      int swap2=Random.randInt(System.x.size());
      while (System.x(swap1)==System.x(swap2))
	swap2=Random.randInt(System.x.size());
      System.Swap(swap1,swap2);
      
      for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++)
	(*wf_iter)->Swap(swap1,swap2);

      complex<double> quick_ratio=1.0;
      for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
	complex<double> myRatioIs=(*wf_iter)->evaluateRatio(System,swap1,swap2);
	//	cerr<<"The ratio is "<<myRatioIs<<endl;
	quick_ratio*=myRatioIs;
      }
      //      exit(1);

      //      complex<double> quick_ratio_check=Psi.evaluateRatio_check(System,swap1,swap2);
      double ranNum=Random.ranf();
      numAttempted++;
      if (2*log(abs(quick_ratio.real())) >log(ranNum)){
	numAccepted++;
	for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++)
	  (*wf_iter)->UpdateDets(System,swap1,swap2);
	//accept
      }
      else {
	System.Swap(swap1,swap2);
	for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++)
	  (*wf_iter)->Swap(swap1,swap2);

	for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++)
	  (*wf_iter)->Reject(System,swap1,swap2) ;
      }
    }
    return (double)numAccepted/(double)numAttempted;
    //cerr<<"Accepted: "<<(double)numAccepted/(double)numAttempted<<endl;
  }

  double VMC(double J2)
   {
     double numAttempt=0;
     double numAccept=0;
     //     ((HuseElserPair2Class*)(&Psi))->evaluate(System);
     //     ((RVBFastPsiClass*)(&Psip))->evaluate(System);
    if (OptimizeCalled==0){
      OptimizeCalled++;
      for (int sweeps=0;sweeps<VMC_equilSweeps;sweeps++){
	Sweep();
      }
    }
    double energy=0.0;
    int NumCounts=0;
    for (int sweeps=0;sweeps<VMC_SampleSweeps;sweeps++){
      numAccept+=Sweep();
      numAttempt+=1;
      double te=H.Energy(System,wf_list)+Hp.Energy(System,wf_list);
      energy+=te;
      NumCounts++;
    }
    double oldEnergy=energy/(double)NumCounts;
    if ((1==1)) { 
      cerr<<"VMCENERGY: "<<energy/(double)NumCounts<<" "<<endl;
      cerr<<"ACCEPT: "<<numAccept/numAttempt<<endl;
      energy=0.0;
      NumCounts=0;
    }
    return oldEnergy;
   }
  

  Array<complex<double>,1> Derivs;
  Array<complex<double>,1> EnergyDerivs;
  Array<complex<double>,1> derivs_times_energy;
  Array<complex<double>,1> derivs;
  Array<complex<double>,1> FullDerivs;
  Array<complex<double >,1> ParamsOld;
  void InitDeriv()
  {
    int NumParams=NumberOfParams;
    Derivs.resize(NumParams);
    EnergyDerivs.resize(NumParams);
    derivs_times_energy.resize(NumParams);
    derivs.resize(NumParams);
    FullDerivs.resize(NumParams);
  }

  void ResetParam()
  {

    int currStart=0;
    for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      WaveFunctionClass &Psi =**wf_iter;
      for (int i=0;i<Psi.NumParams;i++){
	Psi.SetParam_real(i,ParamsOld(currStart+i).real());
	Psi.SetParam_imag(i,ParamsOld(currStart+i).imag());
      }
      currStart=currStart+Psi.NumParams;
    }
  }
  void OptimizeRVB(double J2)
  {

    if (OptimizeCalled==0){

      H.Init(System);
      Hp.Init(System,J2);
      OptimizeCalled++;
    }
    OptimizeCalled++;

    int currStart=0;
    for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      WaveFunctionClass &Psi =**wf_iter;
      for (int i=0;i<Psi.NumParams;i++){
	ParamsOld(currStart+i).real()=Psi.GetParam_real(i);
	ParamsOld(currStart+i).imag()=Psi.GetParam_imag(i);
      }
      currStart=currStart+Psi.NumParams;
    }
    

    DerivClass VarDeriv;
    int NumParams=NumberOfParams;
    VarDeriv.Init(NumParams);
    VarDeriv.Clear();
    derivs=0;
     
    double energyPsi=0;
    double energyPsip=0;
    double energy=0.0;
    int NumCounts=0;
    for (int sweeps=0;sweeps<opt_equilSweeps;sweeps++)
      Sweep();
    for (int sweeps=0;sweeps<opt_SampleSweeps;sweeps++){
      Sweep();

      double te=H.Energy(System,wf_list)+Hp.Energy(System,wf_list);
      int currStart=0;
      for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
	(*wf_iter)->AllDerivs(System,derivs,currStart,currStart+(*wf_iter)->NumParams);
	currStart=currStart+(*wf_iter)->NumParams;
      }
      VarDeriv.Add(te,derivs);
      energy+=te; 
      NumCounts++;
    }

    if (1==1)
    {
      energy=energy/NumCounts;
      cerr<<"ENERGY: "<<" "<<energy<<endl; 
      int currStart=0;
      for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
	WaveFunctionClass &Psi =**wf_iter;
	cerr<<"NEW WAVEFUNCTION"<<endl;
	for (int i=0;i<Psi.NumParams;i++){
	//	for (int i=3;i<4;i++){
	  complex<double> myDeriv=VarDeriv.ComputeDerivp(currStart+i);
	  cerr<<"I IS "<<i<<" "<<myDeriv.real()<<endl;
	  cerr<<"PAIRING: "<<currStart+i<<" "<<" "<<Psi.GetParam_real(i)<<" "<<energy<<" "<<myDeriv.real()<<" "<<myDeriv.imag()<<endl;
	  //	  ParamsOld(currStart+i)=Psi.GetParam(i);
	  //	  if (myDeriv.real()!=0)
	  //	    Psi.SetParam(i,Psi.GetParam(i)+-0.01*abs(myDeriv.real())/myDeriv.real()*Random.ranf());
	  //	  if (Psi.Name!="RVBFastPsi" && myDeriv.real()!=0){
	    //	    Psi.SetParam(i,Psi.GetParam(i)+-0.01*myDeriv.real()); ///abs(myDeriv.real())*Random.ranf());
	    //	    Psi.SetParam(i,Psi.GetParam(i)+-10.0*myDeriv.real()); ///abs(myDeriv.real())*Random.ranf());

	    //	    Psi.SetParam(i,Psi.GetParam(i)+-0.1*(myDeriv.real())/abs(myDeriv.real())*Random.ranf());
	  if (myDeriv.real()!=0)
	    Psi.SetParam_real(i,Psi.GetParam_real(i)+-0.1*(myDeriv.real())/abs(myDeriv.real())*Random.ranf());
	    //	    Psi.SetParam_real(i,Psi.GetParam_real(i)-0.1*(myDeriv.real())); // *(myDeriv.real()));

	  myDeriv=VarDeriv.ComputeDerivp_imag(currStart+i);
	  cerr<<"I imag IS "<<i<<" "<<myDeriv.real()<<endl;
	  cerr<<"PAIRING imag: "<<currStart+i<<" "<<" "<<Psi.GetParam_imag(i)<<" "<<energy<<" "<<myDeriv.real()<<" "<<myDeriv.imag()<<endl;
	  //	  ParamsOld(currStart+i)=Psi.GetParam(i);
	  //	  if (myDeriv.real()!=0)
	  //	    Psi.SetParam(i,Psi.GetParam(i)+-0.01*myDeriv.real()); //*abs(myDeriv.real())/myDeriv.real()*Random.ranf());
	  //	  if (Psi.Name!="RVBFastPsi" && myDeriv.real()!=0){
	    //	    Psi.SetParam(i,Psi.GetParam(i)+-0.01*myDeriv.real()); ///abs(myDeriv.real())*Random.ranf());
	    //	    Psi.SetParam(i,Psi.GetParam(i)+-10.0*myDeriv.real()); ///abs(myDeriv.real())*Random.ranf());

	    //	    Psi.SetParam(i,Psi.GetParam(i)+-0.1*(myDeriv.real())/abs(myDeriv.real())*Random.ranf());
	  if (myDeriv.real()!=0)
	    //	    Psi.SetParam_imag(i,Psi.GetParam_imag(i)-0.1); //*(myDeriv.real()));
	    //asdf	    Psi.SetParam_imag(i,Psi.GetParam_imag(i)+0.1*(myDeriv.real()));
	    //	    Psi.SetParam_imag(i,Psi.GetParam_imag(i)+0.1*(myDeriv.real()));
	    Psi.SetParam_imag(i,Psi.GetParam_imag(i)+0.1*(myDeriv.real())/abs(myDeriv.real())*Random.ranf());

	    //	  else{
	    //	    double newVal=exp(log(Psi.GetParam(i))+-0.01*(myDeriv.real())/abs(myDeriv.real())*Random.ranf());
	    //	    Psi.SetParam(i,newVal);
	    //	  }
	    //	  }
	    //	  else if (myDeriv.real()!=0)
	    //	    Psi.SetParam(i,Psi.GetParam(i)+-0.01*myDeriv.real()); ///abs(myDeriv.real())*Random.ranf());
	}
	currStart=currStart+Psi.NumParams;
      }

      for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
	(*wf_iter)->evaluate(System);
      }
      cerr<<"Variance: "<<VarDeriv.ComputeVariance()<<endl; 
	//      ((HuseElserPair2Class*)(&Psi))->evaluate(System);
	//      ((RVBFastPsiClass*)(&Psip))->evaluate(System);

      energy=0.0;
      NumCounts=0;
      
    
    }


    if ((1==2)) { 
/*       energy=energy/NumCounts; */
/*       cerr<<"ENERGY: "<<" "<<energy<<endl;  */
/*       for (int i=0;i<Psi.PairingFunction.binVals.size();i++){ */
/* 	  complex<double> myDeriv=VarDeriv.ComputeDerivp(i); */
/*      	  cerr<<"PAIRING: "<<i<<" "<<" "<<Psi.PairingFunction.f0[i].real()<<" "<<energy<<" "<<myDeriv.real()<<" "<<endl; */
/* 	  Psi.PairingFunction.f0[i]-=0.01*myDeriv.real(); */
/*       } */
/*       for (int i=0;i<Psi.PairingFunction_diff.binVals.size();i++){ */
/* 	int j=i+Psi.PairingFunction.binVals.size(); */
/* 	complex<double> myDeriv=VarDeriv.ComputeDerivp(j); */
/* 	cerr<<"PAIRING: "<<j<<" "<<" "<<Psi.PairingFunction_diff.f0[i].real()<<" "<<energy<<" "<<myDeriv.real()<<" "<<endl; */
/* 	if (myDeriv.real()!=0) */
/* 	  Psi.PairingFunction_diff.f0[i]-=0.01*myDeriv.real(); ///abs(myDeriv.real()); //0.01*myDeriv.real(); ///abs(myDeriv.real()); */
/*       } */
/*       for (int i=0;i<Psip.PairingFunction.binVals.size();i++){ */
/* 	int j=i+2*Psi.PairingFunction.binVals.size(); */
/*        	complex<double> myDeriv=VarDeriv.ComputeDerivp(j); */
/* 	cerr<<"PAIRING: "<<j<<" "<<" "<<Psip.PairingFunction.f0[i].real()<<" "<<energy<<" "<<myDeriv.real()<<" "<<endl; */
/* 	if (myDeriv.real()!=0 ){ */
/* 	  Psip.PairingFunction.f0[i]-=0.0*0.01*myDeriv.real(); // /abs(myDeriv.real())*Random.ranf(); */
/* 	} */

/*       } */
/*       cerr<<"Variance: "<<VarDeriv.ComputeVariance()<<endl; */
/*       ((HuseElserPair2Class*)(&Psi))->evaluate(System); */
/*       ((RVBFastPsiClass*)(&Psip))->evaluate(System); */

/*       energy=0.0; */
/*       NumCounts=0; */
    }

  }


  
  int Sample(list<pair<SpinSwap,double> > &vals,
	     RandomClass &Random,
	     double &w)
  {
    vector<double> cumulant;
    //    double w=0.0;
    w=0.0;
    for (list<pair<SpinSwap,double> >::iterator iter=vals.begin();iter!=vals.end();iter++){
      //      cerr<<"My vales are "<<(*iter).second<<" "<<(*iter).first.spin1<<" "<<(*iter).first.spin2<<endl;
      if (((*iter).second<-1e-10))
	cerr<<"My vales are "<<(*iter).second<<" "<<(*iter).first.spin1<<" "<<(*iter).first.spin2<<endl;
      assert((*iter).second>-1e-10);
      w+=(*iter).second;
      cumulant.push_back(w);
    }
    //really could just multiply waht we are lookign for by w
    for (int i=0;i<cumulant.size();i++)
      cumulant[i]=cumulant[i]/w;
    double toFind=Random.ranf();
    int i=0;
    while (toFind>cumulant[i])
      i++;
    return i;
  }

 bool Apply(list<pair<SpinSwap,double> > &vals,
	     SystemClass &system,
	    //	     WaveFunctionClass &wf,
	    list<WaveFunctionClass*> wf_list,
	    int num)
    
  {
    //    assert(&Psi==&wf);
    int count=0;
    for (list<pair<SpinSwap,double> >::iterator iter=vals.begin();iter!=vals.end();iter++){
      if (count==num){

	if (iter->first.spin1==0 && iter->first.spin2==0)
	  return false;
	else{
	  system.Swap(iter->first.spin1,iter->first.spin2);
	  
	  int s1=iter->first.spin1;
	  int s2=iter->first.spin2;
	  //	  cerr<<"I am swapping "<<s1<<" "<<s2<<endl;
	  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
	    (*wf)->Swap(s1,s2);
	  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
	    (*wf)->evaluateRatio(System,s1,s2);
	  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
	    (*wf)->UpdateDets(System,s1,s2);
	  return true;
	}
      }
      count++;
    }
  }


 bool Apply_sign(list<pair<SpinSwap,double> > &vals,
		 SystemClass &system,
		 //	     WaveFunctionClass &wf,
		 list<WaveFunctionClass*> &wf_list,
		 vector<int> &signVec,
		 int num)
    
  {
    //    assert(&Psi==&wf);
    int count=0;
    for (list<pair<SpinSwap,double> >::iterator iter=vals.begin();iter!=vals.end();iter++){
      if (count==num){

	if (iter->first.spin1==0 && iter->first.spin2==0)
	  return false;
	else{
	  system.Swap(iter->first.spin1,iter->first.spin2);
	  sign=sign*signVec[count];
	  int s1=iter->first.spin1;
	  int s2=iter->first.spin2;
	  //	  cerr<<"I am swapping "<<s1<<" "<<s2<<endl;
	  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
	    (*wf)->Swap(s1,s2);
	  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
	    (*wf)->evaluateRatio(System,s1,s2);
	  for (list<WaveFunctionClass*>::iterator wf=wf_list.begin();wf!=wf_list.end();wf++)
	    (*wf)->UpdateDets(System,s1,s2);
	  return true;
	}
      }
      count++;
    }
  }





};
#endif
