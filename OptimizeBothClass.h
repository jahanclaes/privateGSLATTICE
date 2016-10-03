#ifndef  OPTIMIZE_BOTH_CLASS_H
#define  OPTIMIZE_BOTH_CLASS_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <map>
using namespace blitz;
using namespace std;
#include "SystemClass.h"
#include "WaveFunction.h"
#include "Hamiltonian.h"
#include "CPS.h"
#include "TripletWF.h"
#include "DerivClass.h"
#include "ProjGutz.h"
#include "SharedWaveFunctionData.h"
#include "input.h"
#include "Timer.h"

enum OptType {GRADIENT, TIMEEVOLUTION, SR};

class OptimizeBothClass
{
public:
  SystemClass System;
  RandomClass &Random;
  list<WaveFunctionClass*> wf_list;  
  list<HamiltonianClass*> Ham;
  Array<complex<double >,1> ParamsOld;
  Array<complex<double>,1> derivs;
  bool computeObservables;
  
  OptType opt;
  double StepSize;
  

  int NumberOfParams;
 
  int VMC_equilSweeps;
  int VMC_SampleSweeps;
  int opt_equilSweeps;
  int opt_SampleSweeps;

  int sign;

  DerivClass VarDeriv;
  
  vector<TimerClass*> myTimers;
  TimerClass InitTimer;


  void Copy(OptimizeBothClass &b)
  {
    ParamsOld=b.ParamsOld;
    derivs=b.derivs;
    sign=b.sign;
    VarDeriv=b.VarDeriv;
    System=b.System;
    list<WaveFunctionClass*>::iterator iter=wf_list.begin();
    for (list<WaveFunctionClass*>::iterator iterb= b.wf_list.begin();
	 iterb!=b.wf_list.end();
	 iterb++){
      (**iter).Copy(*iterb);
      //      (*(RVBFastPsiClass*)(*iter))=(*(RVBFastPsiClass*)(*iterb));
      iter++;
    }
  }
  
  void Combine(vector<OptimizeBothClass*> &vec)
  {
    vector<DerivClass*> vec_Derivs;
    for (int i=0;i<vec.size();i++){
      vec_Derivs.push_back(&(vec[i]->VarDeriv));
    }
    VarDeriv.Combine(vec_Derivs);
    cerr<<"COMBINED: "<<VarDeriv.NumTimes<<endl;
  }




 OptimizeBothClass(RandomClass &t_random) : Random(t_random) , InitTimer("InitTimer")
  {

    
  }

  




   void Init(list<pair<string, SharedWaveFunctionDataClass*> > &wf_init)
   {
     assert(1==2);
     //     std::map<string,string> input;
     //     Init(wf_init,input);
   }

   void Init(list<pair<string, SharedWaveFunctionDataClass*> > &wf_init,std::map<string,string> &input)
   {
     assert(1==2);
   }

   void Init(list<pair<string, SharedWaveFunctionDataClass*> > &wf_init,InputClass &myInput)
   {
     computeObservables=false;
    VMC_equilSweeps=100;
    VMC_SampleSweeps=1000; // 10000*1000;
    opt_equilSweeps=100; //times 40
    opt_SampleSweeps=1000; //times 40

    sign=1;

    cerr<<"Initializing the system"<<endl;
    System.Init();
    System.Stagger();
    

    cerr<<"Initializing the wavefunction"<<endl;
    wf_list.clear();
  
    for (list<pair<string,SharedWaveFunctionDataClass*> >::iterator iter=wf_init.begin();
	 iter!=wf_init.end();iter++){
      string wf_type_string=(*iter).first;
      cerr<<"The wf_type_string is "<<wf_type_string<<endl;
      if (wf_type_string=="CPS"){
	cerr<<"ADDING CPS"<<endl;
	CPSClass *t_CPS=new CPSClass(*((PairingFunctionMany*)((*iter).second)));
	t_CPS->PF.Init(System);
	t_CPS->Init(System);
	wf_list.push_back(t_CPS);
      }
      
    }
    
    cerr<<"Initializing the hamiltonian"<<endl;
    Ham.clear();

    int i=1;


    std::stringstream ss;
    ss << "Heisenberg"<<i;
    string check=ss.str();
    
    cerr<<"Checking Hamiltonian "<<check<<" and getting "<<myInput.IsVariable(check)<<endl;
    //    while (input.count(check)!=0){
    while (myInput.IsVariable(check)){
      Heisenberg *H_temp=new Heisenberg(check);
      H_temp->Init(System,myInput.GetVariable(check));
      std::stringstream ssJ;
      ssJ<<"Heisenberg"<<i<<"_J";
      assert(myInput.IsVariable(ssJ.str()));
      H_temp->Set_J(myInput.toDouble(myInput.GetVariable(ssJ.str())));
      //      H_temp->Set_J(atof(input[ssJ.str()].c_str()));
      cerr<<"The J that is set for "<<check<<" is "<<H_temp->J<<endl;
      Ham.push_back(H_temp);
      std::stringstream ss;
      i++;
      ss << "Heisenberg"<<i;
      check=ss.str();
      cerr<<"Checking: "<<myInput.IsVariable(check)<<endl;
    }
    assert(myInput.IsVariable("OptType"));
    if (myInput.GetVariable("OptType")=="GRADIENT"){
      opt=GRADIENT;
      assert(myInput.IsVariable("StepSize"));
      StepSize=myInput.toDouble(myInput.GetVariable("StepSize"));
      cerr<<"I am optimizing with Gradient and a stepsize of "<<StepSize<<endl;
    }
    else if (myInput.GetVariable("OptType")=="SR"){
      opt=SR;
      assert(myInput.IsVariable("StepSize"));
      StepSize=myInput.toDouble(myInput.GetVariable("StepSize"));
      cerr<<"I am optimizing with SR and a stepsize of "<<StepSize<<endl;
    }
    else if (myInput.GetVariable("OptType")=="TimeEvolution"){
      opt=TIMEEVOLUTION;
      assert(myInput.IsVariable("StepSize"));
      StepSize=myInput.toDouble(myInput.GetVariable("StepSize"));
      cerr<<"I am optimizing with Time Evolution and a stepsize of "<<StepSize<<endl;
    }
    else
      assert(1==2);
    
    cerr<<"Initializing the parameters"<<endl;
    NumberOfParams=0;
    for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      cerr<<"THE NUMBER OF PARAMS IS "<<(*wf_iter)->NumParams<<endl;
      NumberOfParams+=(*wf_iter)->NumParams;
    }
    if (opt==TIMEEVOLUTION)
      VarDeriv.Init(NumberOfParams,true);
    else
      VarDeriv.Init(NumberOfParams,false);

    derivs.resize(NumberOfParams);
    VarDeriv.Clear();
    derivs=0;
    ParamsOld.resize(NumberOfParams);
    cerr<<"Done with init"<<endl;
  }


  double Sweep()
  {
    int numAccepted=0;
    int numAttempted=0;
    for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      if ((*wf_iter)->NeedFrequentReset){
	complex<double> ans=(*wf_iter)->evaluate(System);
      }
    }
    //random swaps
    for (int step=0;step<System.x.size();step++){ 
      int swap1=step;
      int swap2=Random.randInt(System.x.size());
      int reset=0;
      while (System.x(swap1)==System.x(swap2) and reset<50){
	    swap2=Random.randInt(System.x.size());
        reset++;
      }
      if (reset!=50){
          System.Swap(swap1,swap2);
          for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++)
    	    (*wf_iter)->Swap(swap1,swap2);
          complex<double> quick_ratio=1.0;
          for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
    	    complex<double> myRatioIs=(*wf_iter)->evaluateRatio(System,swap1,swap2);
    	    quick_ratio*=myRatioIs;
          }
          double ranNum=Random.ranf();
          
          numAttempted++;
          if ( (2*log(abs(quick_ratio.real())) >log(ranNum))){ //HACK FOR HONEYCOMB
    	    numAccepted++;
    	    for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++)
    	        (*wf_iter)->UpdateDets(System,swap1,swap2);
          }
          else {
    	    System.Swap(swap1,swap2);
    	    for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++)
    	    (*wf_iter)->Swap(swap1,swap2);
    	    for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++)
    	        (*wf_iter)->Reject(System,swap1,swap2) ;
          }
      }
    }
    // Random flips
    for (int step=0;step<System.x.size();step++){
        int flipSpin = Random.randInt(2);
        if (flipSpin==1){
            System.Flip(step);
            complex<double> quick_ratio=1.;
            for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
                complex<double> myRatioIs=(*wf_iter)->evaluateRatioFlip(System,step);
                quick_ratio*=myRatioIs;
            }
            double ranNum=Random.ranf();
            numAttempted++;
            if ( (2*log(abs(quick_ratio.real())) >log(ranNum))){
                numAccepted++;
            }
            else{
                System.Flip(step);
            }
        }
    }
    // Random global flip
    vector<int> flips;
    flips.resize(System.x.size());
    complex<double> quick_ratio = 1.;
    numAttempted++;
    for (int step=0;step<System.x.size();step++){
        int flipSpin = Random.randInt(2);
        if (flipSpin==1){
            flips[step]=1;
            System.Flip(step);
            for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
                complex<double> myRatioIs=(*wf_iter)->evaluateRatioFlip(System,step);
                quick_ratio*=myRatioIs;
            }
        }
        else
            flips[step]=0;
    }
    double ranNum=Random.ranf();
    if ( (2*log(abs(quick_ratio.real())) >log(ranNum))){
        numAccepted++;
    }
    else{
        for (int step=0;step<System.x.size();step++){
            if (flips[step]==1)
                System.Flip(step);
        }
    }

    for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      if ((*wf_iter)->NeedFrequentReset){
	    complex<double> ans=(*wf_iter)->evaluate(System);
      }
    }

    return (double)numAccepted/(double)numAttempted;
  }


void BroadcastParams(CommunicatorClass &myComm)
{
  vector<complex<double> > myParams; 
  for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      WaveFunctionClass &Psi =**wf_iter;
      for (int i=0;i<Psi.NumParams;i++){
	complex<double> myVal(Psi.GetParam_real(i),Psi.GetParam_imag(i));
	myParams.push_back(myVal);
      }
  }
  myComm.Broadcast(0,myParams);
  int countMe=0;
  for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      WaveFunctionClass &Psi =**wf_iter;
      for (int i=0;i<Psi.NumParams;i++){
/* 	if (myComm.MyProc()==0){ */
/* 	  cerr<<"SETTING PARAMETER "<< */
	  
/* 	} */
	
	Psi.SetParam_real(i,myParams[countMe].real());
	Psi.SetParam_imag(i,myParams[countMe].imag());
	countMe++;
      }
  }
  
  

}

 void SaveParams(string fileName)
 {
   ofstream outfile;
   outfile.open(fileName.c_str());
   for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
     WaveFunctionClass &Psi =**wf_iter;
     for (int i=0;i<Psi.NumParams;i++){
       outfile<<Psi.GetParam_real(i)<<" "<<Psi.GetParam_imag(i)<<" "<<endl;
     }
   }
   outfile.close();
 }
 void GetParams(string fileName)
 {
   ifstream infile;
   infile.open(fileName.c_str());
   assert(infile);
   for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
     WaveFunctionClass &Psi =**wf_iter;
     for (int i=0;i<Psi.NumParams;i++){
       double data;
       infile>>data;
       if (!infile.eof()){
	 Psi.SetParam_real(i,data);
	 infile>>data;
	 Psi.SetParam_imag(i,data);
       }
       else{
	 cerr<<"Not enough parameters in file to set parameter number "<<i<<endl;
	 cerr<<"Setting parameter to 1.0"<<endl;
	 Psi.SetParam_real(i,1.0);
	 Psi.SetParam_imag(i,0.0);
       }
     }
   }
   infile.close();
   
 }
 
 void TakeStep(CommunicatorClass &myComm,double cutoff)
 {
    int currStart=0;
    for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
        WaveFunctionClass &Psi =**wf_iter;
        for (int i=0;i<Psi.NumParams;i++){
            ParamsOld(currStart+i).real()=Psi.GetParam_real(i);
	        ParamsOld(currStart+i).imag()=Psi.GetParam_imag(i);
        }
        currStart=currStart+Psi.NumParams;
    }
    currStart=0;
   if (opt==TIMEEVOLUTION)
     VarDeriv.GetSInverse(cutoff);

   for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
     WaveFunctionClass &Psi =**wf_iter;
     if (myComm.MyProc()==0){
       cerr<<"NEW WAVEFUNCTION"<<endl;
     for (int i=0;i<Psi.NumParams;i++){
       complex<double> myDeriv;
       if (opt==TIMEEVOLUTION)
	 myDeriv=VarDeriv.ComputeDerivSR(currStart+i);
       else 
	 myDeriv=VarDeriv.ComputeDerivp(currStart+i);
       if (myDeriv.real()!=0)
	 cerr<<"My parameter deriv is "<<myDeriv<<" "<<(myDeriv.real())/abs(myDeriv.real())<<" "<<Psi.GetParam_real(i)<<endl;
       if (opt==GRADIENT || opt==TIMEEVOLUTION)
	 Psi.SetParam_real(i,Psi.GetParam_real(i)+-StepSize*(myDeriv.real())); 
       else if (opt==SR && myDeriv.real()!=0)
	 Psi.SetParam_real(i,Psi.GetParam_real(i)+-StepSize*(myDeriv.real())/abs(myDeriv.real())*Random.ranf());
       else if (myDeriv.real()==0){
       }
       else 
	 assert(1==2);
       /* 	  myDeriv=VarDeriv.ComputeDerivp_imag(currStart+i); */
       /* 	  //HACK!	  cerr<<"I imag IS "<<i<<" "<<Psi.GetParam_imag(i)<<" "<<myDeriv.real()<<endl; */
       /* 	  //	  cerr<<"PAIRING imag: "<<currStart+i<<" "<<" "<<Psi.GetParam_imag(i)<<" "<<energy<<" "<<myDeriv.real()<<" "<<myDeriv.imag()<<endl; */
       /* 	  if (myDeriv.real()!=0) */
       /* 	    Psi.SetParam_imag(i,Psi.GetParam_imag(i)+0.01*(myDeriv.real())); // /abs(myDeriv.real())*Random.ranf()); */
       /* 	    //	    Psi.SetParam_imag(i,Psi.GetParam_imag(i)+0.001*(myDeriv.real())/abs(myDeriv.real())*Random.ranf()); */
       /* 	} */
     }
     currStart=currStart+Psi.NumParams;
     }
   }
   EvaluateAll();
   if (myComm.MyProc()==0)
     cerr<<"Variance: "<<VarDeriv.ComputeVariance()<<endl; 
   //      energy=0.0;
   //      NumCounts=0;
   
   
 }
 void EvaluateAll()
 {
   //HACK!
   //    return; 
   for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
     (*wf_iter)->evaluate(System);
   }
   
 }
 
 
 /*   void SetParam(OptimizeBothClass &opt) */
 /*   { */
 
 /*     int currStart=0; */
 /*     for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){ */
/*       WaveFunctionClass &Psi =**wf_iter; */
/*       for (int i=0;i<Psi.NumParams;i++){ */
/* 	Psi.SetParam_real(i,ParamsOld(currStart+i).real()); */
/* 	Psi.SetParam_imag(i,ParamsOld(currStart+i).imag()); */
/*       } */
/*       currStart=currStart+Psi.NumParams; */
/*     } */
/*   } */

  int toInt(){
    int result=0;
    for (int i=0;i<System.x.size();i++){
        if( System.x(i)==1){
            result += (1 << (System.x.size()-i-1));
        }
    }
    return result;
    }

  void toSystem(int sysInt){
    complex<double> result=0;
    int current = sysInt;
    for (int i=0;i<System.x.size();i++){
        if (current >=(1 << (System.x.size()-i-1))){
            System.x(i)=1;
            current-=(1 << (System.x.size()-i-1));
        }
        else
            System.x(i)=0;
    }
  }
  

  double Optimize()
  {
       
    
    VarDeriv.Clear();
    derivs=0;
    double acceptanceRatio = 0.0;
    
    ///FIXME!
    
    double energy=0.0;
    int NumCounts=0;
    for (int sweeps=0;sweeps<opt_equilSweeps;sweeps++)
      Sweep();
    for (int sweeps=0;sweeps<opt_SampleSweeps;sweeps++){
      acceptanceRatio+=Sweep();
      double te=0;
      for (list<HamiltonianClass*>::iterator iter=Ham.begin();iter!=Ham.end();iter++){
	double tte=(*iter)->Energy(System,wf_list);
	te+=tte;
      }
      int currStart=0;
      for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
	(*wf_iter)->AllDerivs(System,derivs,currStart,currStart+(*wf_iter)->NumParams);
	currStart=currStart+(*wf_iter)->NumParams;
      }
      VarDeriv.Add(te,derivs);
      energy+=te; 
      NumCounts++;
    }

    if (1==2)
    {
      energy=energy/NumCounts;
      cerr<<"ENERGY: "<<" "<<energy<<endl; 
      int currStart=0;
      for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
	WaveFunctionClass &Psi =**wf_iter;
	cerr<<"NEW WAVEFUNCTION"<<endl;
	for (int i=0;i<Psi.NumParams;i++){
	  complex<double> myDeriv=VarDeriv.ComputeDerivp(currStart+i);
	  cerr<<"I IS "<<i<<" "<<myDeriv.real()<<endl;
	  cerr<<"PAIRING: "<<currStart+i<<" "<<" "<<Psi.GetParam_real(i)<<" "<<energy<<" "<<myDeriv.real()<<" "<<myDeriv.imag()<<endl;
	  if (myDeriv.real()!=0)
	    Psi.SetParam_real(i,Psi.GetParam_real(i)+-0.01*(myDeriv.real())/abs(myDeriv.real())*Random.ranf());
	  myDeriv=VarDeriv.ComputeDerivp_imag(currStart+i);
	  cerr<<"I imag IS "<<i<<" "<<myDeriv.real()<<endl;
	  cerr<<"PAIRING imag: "<<currStart+i<<" "<<" "<<Psi.GetParam_imag(i)<<" "<<energy<<" "<<myDeriv.real()<<" "<<myDeriv.imag()<<endl;
	  if (myDeriv.real()!=0)
	    Psi.SetParam_imag(i,Psi.GetParam_imag(i)+0.01*(myDeriv.real())/abs(myDeriv.real())*Random.ranf());
	}
	currStart=currStart+Psi.NumParams;
      }

      for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
	(*wf_iter)->evaluate(System);
      }
      cerr<<"Variance: "<<VarDeriv.ComputeVariance()<<endl; 
      energy=0.0;
      NumCounts=0;
      
    
    }
    ///DONE FIXME
    return acceptanceRatio/opt_SampleSweeps;



    
  }



 bool Apply(list<pair<SpinSwap,double> > &vals,
	    SystemClass &system,
	    list<WaveFunctionClass*> wf_list,
	    int num)
    
  {
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


 double VMC(bool equilibrate,
	    ofstream* outfile=NULL)
  {
    //    cerr<<"Started VMC"<<endl;
    double numAttempt=0;
    double numAccept=0;
    cerr << "Starting VMC "<<endl;
    if (equilibrate){
      for (int sweeps=0;sweeps<VMC_equilSweeps;sweeps++){
	    Sweep();
      }
    }
    vector<double> energy_terms(Ham.size(),0);
    double energy=0.0;
    int NumCounts=0;
    double observable=0.;
    for (int sweeps=0;sweeps<VMC_SampleSweeps;sweeps++){
      numAccept+=Sweep();
      if ((outfile!=NULL) && computeObservables){
        observable+=StaggeredMagnetization();
      }
      numAttempt+=1;
      NumCounts++;
      int i=0;
      for (list<HamiltonianClass*>::iterator iter=Ham.begin();iter!=Ham.end();iter++){
	double te=(*iter)->Energy(System,wf_list);
	energy_terms[i]+=te;
	i++;
	energy+=te;

      }

    }
    int NumParams=(*(wf_list.begin()))->NumParams;

    double oldEnergy=energy/(double)NumCounts;
    if ((1==1)) { 
      if (outfile!=NULL && (computeObservables==false)){
	    (*outfile)<< (energy/(double)NumCounts) <<" "<<endl;
      }
      if (outfile!=NULL && computeObservables){
	    (*outfile)<< (observable/(double)NumCounts) <<" "<<endl;
        (*outfile).flush();
      }
      (*(Ham.begin()))->term1=0; (*(Ham.begin()))->term2=0;

      for (int i=0;i<energy_terms.size();i++)
      energy=0.0;
      NumCounts=0;
    }
    return oldEnergy;
  }
  
double StaggeredMagnetization(){
    vector<int> s(System.x.size());
    s[0] = 1;
    for (int i=0; i<System.neighbors.size(); i++){
        for (int j=0; j<System.neighbors[i].size(); j++){
            s[System.neighbors[i][j]] = (-1)*s[i];
        }
    }
    double sl_local = 0.0;
    for (int i=0;i<System.x.size();i++){
        sl_local+=s[i]*(System.x(i)==0 ? 1: -1);
    }
    return sl_local*sl_local;
 }
 
 void RestoreParamsOld(){
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


 





};

/*     } */

#endif
