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
#include "PEPS.h"
#include "cPEPS.h"
#include "TripletWF.h"
#include "DerivClass.h"
#include "ProjGutz.h"
#include "RVB_fast.h"
#include "RVBp.h"
#include "SharedWaveFunctionData.h"
#include "input.h"
#include "Timer.h"
#include "BackFlow.h"
#include "SlaterDet.h"
#include "Jastrow.h"
#include "HopMove.h"
#include "kondoHelp.h"

enum OptType {GRADIENT, TIMEEVOLUTION, SR};

enum SweepType {HOP, EXCHANGE,KONDO};

class OptimizeBothClass
{
public:
  SystemClass System;
  RandomClass &Random;
  list<WaveFunctionClass*> wf_list;  
  KondoHelp kondo_help;  
  SweepType sweep;
  //  vector<list<WaveFunctionClass*> > wf_list_;  

  list<HamiltonianClass*> Ham;
  Array<complex<double >,1> ParamsOld;
  Array<complex<double>,1> derivs;
  
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


  void TestDerivs(int derivInt,WaveFunctionClass &Psi)
  {
    cerr<<"Testing the derivative "<<derivInt<<endl;
    double currParam=Psi.GetParam_real(derivInt);
    for (double myStep=-0.01;myStep<0.011;myStep+=0.01){
      double step_energy=0.0;
      double countSteps=0.0;
      Psi.SetParam_real(derivInt,currParam+myStep);
      for (int sweeps=0;sweeps<opt_equilSweeps;sweeps++)
	Sweep_hop();

      for (int sweeps=0;sweeps<opt_SampleSweeps*1000;sweeps++){
	Sweep_hop();
      
	double te=0;
	for (list<HamiltonianClass*>::iterator iter=Ham.begin();iter!=Ham.end();iter++){
	  double tte=(*iter)->Energy(System,wf_list);
	  te+=tte;
	}
	countSteps+=1;
	step_energy+=te;
      }
      
      cerr<<"The energy for step "<<myStep<<" is "<<step_energy/countSteps<<endl;
    }
    
  }

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
    //    cerr<<"COMBINED: "<<VarDeriv.NumTimes<<endl;
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

    VMC_equilSweeps=100;
    VMC_SampleSweeps=1000; // 10000*1000;
    opt_equilSweeps=100; //times 40
    opt_SampleSweeps=1000; //times 40

    sign=1;

    cerr<<"Initializing the system"<<endl;
    System.Init(myInput);
    //    System.Stagger();
    

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
      else if (wf_type_string=="PEPS"){
	cerr<<"ADDING PEPS"<<endl;
	cPEPSClass *t_PEPS=new cPEPSClass();

	int L,Dx,Dy,W,chi;
	double tol;

	int check=0;
	while (myInput.OpenSection("WaveFunction",check)){
	  string waveFunction=myInput.GetVariable("name");
	  cerr<<"Found the "<<check<<" wavefunction called"<<waveFunction<<endl;

	  if (waveFunction=="PEPS"){
	    Dx=myInput.toInteger(myInput.GetVariable("Dx"));
	    Dy=myInput.toInteger(myInput.GetVariable("Dy"));
	    chi=myInput.toInteger(myInput.GetVariable("chi"));
	    L=myInput.toInteger(myInput.GetVariable("L"));
	    W=myInput.toInteger(myInput.GetVariable("W"));
	    tol=myInput.toDouble(myInput.GetVariable("tol"));
	  }
	  myInput.CloseSection();
	  check++;
	}
	
	t_PEPS->Init(System,L,W,4,Dx,Dy,chi,tol);

      	wf_list.push_back(t_PEPS);
      }




      else if (wf_type_string=="RVB"){
	cerr<<"ADDING RVB"<<endl;
	//	RVBFastPsiClass *RVB=new RVBFastPsiClass(*((PairingFunctionAllBin*)((*iter).second)));
	{
	  RVBpPsiClass *RVB=new RVBpPsiClass(*((PairingFunctionAllBin*)((*iter).second)));
	  RVB->Init(System);
	  wf_list.push_back(RVB); 
	}

	{
	  //	  RVBFastPsiClass *RVB=new RVBFastPsiClass(*((PairingFunctionAllBin*)((*iter).second)));
	  //	  RVB->Init(System);
	  //	  wf_list.push_back(RVB); 
	}


      }
      else if (wf_type_string=="JASTROW"){
	cerr<<"ADDING JASTROW"<<endl;
	JastrowClass  *t_JASTROW  = new JastrowClass();
	t_JASTROW->Init(System);
	wf_list.push_back(t_JASTROW);
      }

      else if (wf_type_string=="BACKFLOW"){
        cerr<<"ADDING BACKFLOW"<<endl;
        BackFlowClass *t_BACKFLOW=new BackFlowClass();
        t_BACKFLOW->Init(System);
        wf_list.push_back(t_BACKFLOW);
      }
      
      else if (wf_type_string=="SLATERDET"){
        cerr<<"ADDING Slater Det"<<endl;
	{
	  SlaterDetPsiClass *t_SLATERDET=new SlaterDetPsiClass(*((SharedEigsClass*)((*iter).second)));
	  t_SLATERDET->Init(System,1);
	  wf_list.push_back(t_SLATERDET);
	}
	{
	  SlaterDetPsiClass *t_SLATERDET=new SlaterDetPsiClass(*((SharedEigsClass*)((*iter).second)));
	  t_SLATERDET->Init(System,-1);
	  wf_list.push_back(t_SLATERDET);
	}
      }
      else if (wf_type_string=="SLATERDETUP"){
        cerr<<"ADDING Slater Det up"<<endl;
	{
	  SlaterDetPsiClass *t_SLATERDET=new SlaterDetPsiClass(*((SharedEigsClass*)((*iter).second)));
	  t_SLATERDET->Init(System,1);
	  wf_list.push_back(t_SLATERDET);
	}

      }
      else if (wf_type_string=="SLATERDETDOWN"){
        cerr<<"ADDING Slater Det down"<<endl;
	{
	  SlaterDetPsiClass *t_SLATERDET=new SlaterDetPsiClass(*((SharedEigsClass*)((*iter).second)));
	  t_SLATERDET->Init(System,-1);
	  wf_list.push_back(t_SLATERDET);
	}
      }

      

    }
    
    cerr<<"Initializing the hamiltonian"<<endl;
    Ham.clear();

    int check=0;
    while (myInput.OpenSection("Hamiltonian",check)){
      string hamiltonian_Name=myInput.GetVariable("name");
      cerr<<"Found the "<<check<<" hamiltonian called "<<hamiltonian_Name<<endl;
      if (hamiltonian_Name=="Hubbard"){
	std::stringstream ss; 
	ss << "Hubbard"<<check; 
	string check=ss.str(); 
	Hubbard *H_temp=new Hubbard(check);
	H_temp->Init(System,myInput.GetVariable("bondFile"));
	//	std::stringstreamz ssJ;
	//	ssJ<<"Hubbard"<<i<<"_U";
	//	assert(myInput.IsVariable(ssJ.str()));
	assert(myInput.IsVariable("U"));
	assert(myInput.IsVariable("t"));
	//	H_temp->Set_U(myInput.toDouble(myInput.GetVariable(ssJ.str())));
	H_temp->Set_U(myInput.toDouble(myInput.GetVariable("U")));
	H_temp->t=myInput.toDouble(myInput.GetVariable("t")); //Set_U(myInput.toDouble(myInput.GetVariable("U")));
	//      H_temp->Set_J(atof(input[ssJ.str()].c_str()));
	cerr<<"The U that is set for "<<check<<" is "<<H_temp->U<<endl;
	Ham.push_back(H_temp);
      } 
      else if (hamiltonian_Name=="Heisenberg"){
	std::stringstream ss; 
	ss << "Heisenberg"<<check; 
	string check=ss.str(); 
	Heisenberg *H_temp=new Heisenberg(check);
	H_temp->Init(System,myInput.GetVariable("bondFile"));
	H_temp->Set_J(myInput.toDouble(myInput.GetVariable("J")));
	Ham.push_back(H_temp);
      }
      myInput.CloseSection();
      check++;
    }


    /* int i=1; */


    /* std::stringstream ss; */
    /* ss << "Hubbard"<<i; */
    /* string check=ss.str(); */
    
    /* cerr<<"Checking Hamiltonian "<<check<<" and getting "<<myInput.IsVariable(check)<<endl; */
    /* //    while (input.count(check)!=0){ */
    /* while (myInput.IsVariable(check)){ */
    /*   Hubbard *H_temp=new Hubbard(check); */
    /*   H_temp->Init(System,myInput.GetVariable(check)); */
    /*   std::stringstream ssJ; */
    /*   ssJ<<"Hubbard"<<i<<"_U"; */
    /*   assert(myInput.IsVariable(ssJ.str())); */
    /*   H_temp->Set_U(myInput.toDouble(myInput.GetVariable(ssJ.str()))); */
    /*   //      H_temp->Set_J(atof(input[ssJ.str()].c_str())); */
    /*   cerr<<"The U that is set for "<<check<<" is "<<H_temp->U<<endl; */
    /*   Ham.push_back(H_temp); */
    /*   std::stringstream ss; */
    /*   i++; */
    /*   ss << "Hubbard"<<i; */
    /*   check=ss.str(); */
    /*   cerr<<"Checking: "<<myInput.IsVariable(check)<<endl; */
    /* } */



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
   //Have to have doS set correctly to call the reset
   void ResetParams()
   {
     cerr<<"Resetting the parameters"<<endl;
     NumberOfParams=0;
     for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
       cerr<<"THE NUMBER OF PARAMS IS "<<(*wf_iter)->NumParams<<endl;
       NumberOfParams+=(*wf_iter)->NumParams;
     }
     VarDeriv.Init(NumberOfParams,VarDeriv.doS);
     derivs.resize(NumberOfParams);
     VarDeriv.Clear();
     derivs=0;
     ParamsOld.resize(NumberOfParams);
     cerr<<"Done with reset"<<endl;
   }
   

  double Sweep()
  {
    //    cerr<<"Start sweep "<<endl;
    int numAccepted=0;
    int numAttempted=0;
    for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      if ((*wf_iter)->NeedFrequentReset){
	complex<double> ans=(*wf_iter)->evaluate(System);
      }
    }
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
	quick_ratio*=myRatioIs;
      }
      double ranNum=Random.ranf();
      
      numAttempted++;
      if ( (2*log(abs(quick_ratio.real())) >log(ranNum))){ //HACK FOR HONEYCOMB
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
    //    cerr<<"Accepted: "<<(double)numAccepted/(double)numAttempted<<endl;
    //    cerr<<"end sweep "<<endl;
    return (double)numAccepted/(double)numAttempted;
    //
  }

  complex<double> EvaluateWF_ifNeedReset(list<WaveFunctionClass*> &wf_list)
  {
    complex<double> ans=1.0;
    for (auto wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      if ((*wf_iter)->NeedFrequentReset){
	complex<double> val=(*wf_iter)->evaluate(System);
	ans*=val;
      }
    }
    return ans;
  }

  HopMove ChooseHop(vector<int> &positions,SystemClass &system)
  {
    int spin = (Random.randInt(2)==0 ? -1: 1);
    int site;
    do { 
      site=positions[Random.randInt(positions.size())];
    } while ((system.x(site)!=spin) && (system.x(site)!=2)) ;
    int end_site;
    do {
      end_site=positions[Random.randInt(positions.size())];
    } while ((system.x(end_site)==spin) || (system.x(end_site)==2));
    //    HopMove move(site,end_site,spin);
    return HopMove(site,end_site,spin);
  }
  void MakeMove(HopMove hop, SystemClass &system, list<WaveFunctionClass*> &wf)
  {
    System.Move(hop.start,hop.end,hop.spin);
    for (auto wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      (*wf_iter)->Move(hop.start,hop.end,hop.spin);
    }
  }

  complex<double> Ratio(HopMove &hop, SystemClass &system, list<WaveFunctionClass*> &wf)
  {
    complex<double> quick_ratio=1.0;
    for (auto wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      complex<double> myRatioIs=(*wf_iter)->evaluateRatio(System,hop.start,hop.end,hop.spin);
      quick_ratio*=myRatioIs;
    }
    return quick_ratio;
  }

  void Accept(HopMove &hop, SystemClass &system, list<WaveFunctionClass*> &wf_list)
  {
    for (auto wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      (*wf_iter)->UpdateDets(System,hop.start,hop.end,hop.spin);
    }
  }
  
  void Reject(HopMove &hop, SystemClass &system, list<WaveFunctionClass*> &wf_list)
  {
    for (auto wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      (*wf_iter)->Reject(System,hop.start,hop.end,hop.spin);
    }
  }

  HopMove ReverseHop(HopMove &hop){
    return HopMove(hop.end,hop.start,hop.spin);
  }

  double Sweep_kondo_hop()
  {
    bool noHop=true;
    for (int i=1;i<kondo_help.layer1Sites.size();i++)
      if (System.x(kondo_help.layer1Sites[0])!=System.x(kondo_help.layer1Sites[i]))
	noHop=false;
    if (noHop)
      return 0.0;
    int numAccepted=0;
    int numAttempted=0;
    EvaluateWF_ifNeedReset(wf_list);
    HopMove hop_move;
    for (int step=0;step<System.x.size();step++){ 
      hop_move=ChooseHop(kondo_help.layer1Sites,System);
      //      if (Random.randInt(2) ==0)
      //	hop_move=ChooseHop(kondo_help.layer1Sites,System);
      //      else 
      //	hop_move=ChooseHop(kondo_help.layer2Sites,System);
      //      cerr<<"I am moving from "<<hop_move.start<<" to "<<hop_move.end<<" with "<<hop_move.spin<<endl;
      MakeMove(hop_move,System,wf_list);
      complex<double> quick_ratio=Ratio(hop_move,System,wf_list);
      //      cerr<<"My ratio is "<<quick_ratio<<endl;
      double ranNum=Random.ranf();
      numAttempted++;
      if ( (2*log(abs(quick_ratio.real())) >log(ranNum))){ 
	numAccepted++;
	Accept(hop_move,System,wf_list);

	//accept
      }
      else {
	MakeMove(ReverseHop(hop_move),System,wf_list);
	Reject(hop_move,System,wf_list);

      }
    }
    //    cerr<<"Accepted: "<<(double)numAccepted/(double)numAttempted<<endl;
    return (double)numAccepted/(double)numAttempted;
    //
  }
  bool LegalExchange(ExchangeMove &e,SystemClass &system)
  {
    return ((e.site1 != e.site2) &&
      ( 
      ( (system.x(e.site1)==1) && (system.x(e.site2)==-1)) || 
      ( (system.x(e.site1)==-1) && (system.x(e.site2)==1)) 
	)
	    );
  }

  ExchangeMove ChooseExchange(vector<pair<int, int > > &bonds)
  {
    
    ExchangeMove e;
    //    do {
      int bondIndex=Random.randInt(bonds.size());
      e.site1=bonds[bondIndex].first;
      e.site2=bonds[bondIndex].second;
      //    } while (!(LegalExchange(e,System)));
    e.spin1=System.x(e.site1);
    e.spin2=System.x(e.site2);
    return e;
  }
  
  void MakeMove(ExchangeMove &e, SystemClass &system,list<WaveFunctionClass*> &wf_list)
  {
    System.Swap(e.site1,e.site2);
    for (auto wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      (*wf_iter)->Swap(e.site1,e.site2);
    }
  }
  
  complex<double> Ratio(ExchangeMove &e, SystemClass &system, list<WaveFunctionClass*> &wf_list)
  {
    complex<double> quick_ratio=1.0;
    for (auto wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      quick_ratio*=(*wf_iter)->evaluateRatio(system,e.site1,e.site2);
    }
    return quick_ratio;
  }
  void Accept(ExchangeMove &e, SystemClass &system, list<WaveFunctionClass*> &wf_list)
  {
    for (auto wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      (*wf_iter)->UpdateDets(system,e.site1,e.site2);
    }
  }

  void Reject(ExchangeMove &e, SystemClass &system, list<WaveFunctionClass*> &wf_list)
  {
    for (auto wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      (*wf_iter)->Reject(system,e.site1,e.site2);
    }
  }
    
  
  double Sweep_kondo_exchange()
  {


    int numAccepted=0;
    int numAttempted=0;
    EvaluateWF_ifNeedReset(wf_list);
    ExchangeMove exchange;

    for (int step=0;step<System.x.size()*System.x.size();step++){ 
      ExchangeMove exchange_move=ChooseExchange(kondo_help.heisenbergBonds);
      if (LegalExchange(exchange_move,System)){
	//	cerr<<"I am going to try "<<exchange_move.site1<<" to "<<exchange_move.site2<<" with spin "<<exchange_move.spin1<<" to "<<exchange_move.spin2<<endl;
	MakeMove(exchange_move,System,wf_list);
	complex<double> quick_ratio=Ratio(exchange_move,System,wf_list);
	double ranNum=Random.ranf();
	numAttempted++;
	if ( (2*log(abs(quick_ratio.real())) >log(ranNum))){ 
	  numAccepted++;
	  Accept(exchange_move,System,wf_list);

	  //accept
	}
	else {
	  MakeMove(exchange_move,System,wf_list);
	  Reject(exchange_move,System,wf_list);

	}
      }

    }

    //    cerr<<"Accepted: "<<(double)numAccepted/(double)numAttempted<<endl;

    return (double)numAccepted/(double)numAttempted;
    //
  }



  double Sweep_hop()
  {
    if (sweep==EXCHANGE)
      return Sweep();
    else if (sweep==KONDO){
      double acceptProb_hop=Sweep_kondo_hop();
      //      cerr<<"Hop prob is "<<acceptProb_hop<<endl;
      double acceptProb_exchange=Sweep_kondo_exchange();
      //      cerr<<"exchange prob is "<<acceptProb_exchange<<endl;
      //      double acceptProb_exchange=0.0;

      return (acceptProb_hop+acceptProb_exchange)/2.0;
    }
      //    cerr<<endl;
    int numAccepted=0;
    int numAttempted=0;
    for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      if ((*wf_iter)->NeedFrequentReset){
	complex<double> ans=(*wf_iter)->evaluate(System);
      }
    }
    for (int step=0;step<System.x.size();step++){ 
      //      cerr<<"On step "<<step<<endl;

      int spin = (Random.randInt(2) == 0 ? -1: 1);
      int site=Random.randInt(System.x.size());
      while ( (System.x(site)!=spin && System.x(site)!=2))
        site=Random.randInt(System.x.size());

      int end_site=Random.randInt(System.x.size());
      while ( (System.x(end_site)==spin || System.x(end_site)==2))
        end_site=Random.randInt(System.x.size());

      System.Move(site,end_site,spin);
      for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++)
	(*wf_iter)->Move(site,end_site,spin);
      complex<double> quick_ratio=1.0;
      for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
	complex<double> myRatioIs=(*wf_iter)->evaluateRatio(System,site,end_site,spin);
	quick_ratio*=myRatioIs;
      }
      //      cerr<<"My ratio is "<<quick_ratio.real()<<endl;
      double ranNum=Random.ranf();
      numAttempted++;
      if ( (2*log(abs(quick_ratio.real())) >log(ranNum))){ 
	numAccepted++;
	for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++)
	  (*wf_iter)->UpdateDets(System,site,end_site,spin);
	
	//accept
	
      }
      else {
	System.Move(end_site,site,spin);
	for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++)
	  (*wf_iter)->Move(end_site,site,spin);
	
 	for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++)
	  (*wf_iter)->Reject(System,site,end_site,spin);
	
      }
    }
    //    cerr<<"Accepted: "<<(double)numAccepted/(double)numAttempted<<endl;
    return (double)numAccepted/(double)numAttempted;
    //
  }

  double VMC_correlated(bool equilibrate)
  {
    cerr<<"Started VMC"<<endl;
    double numAttempt=0;
    double numAccept=0;
    if (equilibrate){
      for (int sweeps=0;sweeps<VMC_equilSweeps;sweeps++){
	Sweep_hop();
      }
    }


    vector<list<WaveFunctionClass*> > wf_list_corr;      
    wf_list_corr.resize(3);
    for (int i=0;i<3;i++){
      for (list<WaveFunctionClass*>::iterator iter=wf_list.begin();iter!=wf_list.end();iter++){
	wf_list_corr[i].push_back((*iter)->clone());
      iter++;
      }
    }
    
    vector<double> energy_terms(Ham.size(),0);
    double energy=0.0;
    vector<double> energy_corr(wf_list_corr.size(),0.0);
    int NumCounts=0;
    for (int sweeps=0;sweeps<VMC_SampleSweeps;sweeps++){
      numAccept+=Sweep_hop_correlated(wf_list_corr);
      numAttempt+=1;
      NumCounts++;
      int i=0;
      for (list<HamiltonianClass*>::iterator iter=Ham.begin();iter!=Ham.end();iter++){
	double te=(*iter)->Energy(System,wf_list);
	energy_terms[i]+=te;
	i++;
	energy+=te;
	for (int j=0;j<energy_corr.size();j++){
	  energy_corr[j]+=(*iter)->Energy(System,wf_list_corr[j]);
	}
	

      }
    }
    int NumParams=(*(wf_list.begin()))->NumParams;

    double oldEnergy=energy/(double)NumCounts;
    if ((1==1)) { 
      cerr<<"VMCENERGY: "<<energy/(double)NumCounts<<" "<<endl;
      cerr<<"J1 J2: "<<((*(Ham.begin()))->term1)/(double)NumCounts<<" "<<((*(Ham.begin()))->term2)/(double)NumCounts<<endl;
      (*(Ham.begin()))->term1=0; (*(Ham.begin()))->term2=0;

      cerr<<"ENERGY TERMS: ";
      for (int i=0;i<energy_terms.size();i++)
	cerr<<energy_terms[i]/(double)NumCounts<<" ";
      cerr<<endl;

      cerr<<"ACCEPT: "<<numAccept/numAttempt<<endl;
      energy=0.0;
      NumCounts=0;
    }
    cerr<<"Ended VMC"<<endl;
    return oldEnergy;
    


  }
 
  double Sweep_hop_correlated(vector<list<WaveFunctionClass*> > wf_list_corr)
  {
/*     int numAccepted=0; */
/*     int numAttempted=0; */
    
    
/*     for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){ */
/*       if ((*wf_iter)->NeedFrequentReset){ */
/* 	complex<double> ans=(*wf_iter)->evaluate(System); */
/*       } */
/*     } */
/*     for (int step=0;step<System.x.size();step++){  */
/*       //      cerr<<"On step "<<step<<endl; */

/*       int spin = (Random.randInt(2) == 0 ? -1: 1); */
/*       int site=Random.randInt(System.x.size()); */
/*       while ( (System.x(site)!=spin && System.x(site)!=2)) */
/*         site=Random.randInt(System.x.size()); */

/*       int end_site=Random.randInt(System.x.size()); */
/*       while ( (System.x(end_site)==spin || System.x(end_site)==2)) */
/*         end_site=Random.randInt(System.x.size()); */

/*       System.Move(site,end_site,spin); */
/*       for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++) */
/* 	(*wf_iter)->Move(site,end_site,spin); */
/*       complex<double> quick_ratio=1.0; */
/*       for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){ */
/* 	complex<double> myRatioIs=(*wf_iter)->evaluateRatio(System,site,end_site,spin); */
/* 	quick_ratio*=myRatioIs; */
/*       } */
/*       //      cerr<<"My ratio is "<<quick_ratio.real()<<endl; */
/*       double ranNum=Random.ranf(); */
/*       numAttempted++; */
/*       if ( (2*log(abs(quick_ratio.real())) >log(ranNum))){  */
/* 	numAccepted++; */
/* 	for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++) */
/* 	  (*wf_iter)->UpdateDets(System,site,end_site,spin); */
	
/* 	//accept */
	
/*       } */
/*       else { */
/* 	System.Move(end_site,site,spin); */
/* 	for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++) */
/* 	  (*wf_iter)->Move(end_site,site,spin); */
	
/*  	for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++) */
/* 	  (*wf_iter)->Reject(System,site,end_site,spin); */
	
/*       } */
/*     } */
/*     //    cerr<<"Accepted: "<<(double)numAccepted/(double)numAttempted<<endl; */
/*     return (double)numAccepted/(double)numAttempted; */
    //
  }




  void MatchParams(OptimizeBothClass &vmc)
  {
  vector<complex<double> > myParams; 
  for (list<WaveFunctionClass*>::iterator wf_iter=vmc.wf_list.begin();wf_iter!=vmc.wf_list.end();wf_iter++){
      WaveFunctionClass &Psi =**wf_iter;
      for (int i=0;i<Psi.NumParams;i++){
	complex<double> myVal(Psi.GetParam_real(i),Psi.GetParam_imag(i));
	myParams.push_back(myVal);
      }
  }
  int countMe=0;
  for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      WaveFunctionClass &Psi =**wf_iter;
      for (int i=0;i<Psi.NumParams;i++){
	Psi.SetParam_real(i,myParams[countMe].real());
	Psi.SetParam_imag(i,myParams[countMe].imag());
	countMe++;
      }
  }
  


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
   outfile.precision(std::numeric_limits<double>::digits10 + 1);
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

 
 void GetGradient(CommunicatorClass &myComm,
		  vector<complex<double> > &gradient)
 {
   gradient.resize(NumberOfParams);

   for (int i=0;i<gradient.size();i++){
     gradient[i]=0.0;
   }
   
   int currStart=0;
   if (opt==TIMEEVOLUTION)
     VarDeriv.GetSInverse();
   for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
     WaveFunctionClass &Psi =**wf_iter;
     for (int i=0;i<Psi.NumParams;i++){
       complex<double> myDeriv;
       if (opt==TIMEEVOLUTION)
	 myDeriv=VarDeriv.ComputeDerivSR(currStart+i);
       else 
	 myDeriv=VarDeriv.ComputeDerivp(currStart+i);
       if (opt==SR && myDeriv!=0.0)
	 myDeriv=myDeriv/abs(myDeriv);
       gradient[currStart+i]=myDeriv;
     }
     currStart=currStart+Psi.NumParams;
   }
 }
 void TakeStep_old(CommunicatorClass &myComm)
 {
   vector<complex<double> > gradient;
   GetGradient(myComm,gradient);
   
   

 }
 
 void TakeStep(CommunicatorClass &myComm)
 {
   int currStart=0;
   if (opt==TIMEEVOLUTION)
     VarDeriv.GetSInverse();

   for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
     WaveFunctionClass &Psi =**wf_iter;
     if (myComm.MyProc()==0){
     for (int i=0;i<Psi.NumParams;i++){
       complex<double> myDeriv;
       if (opt==TIMEEVOLUTION)
	 myDeriv=VarDeriv.ComputeDerivSR(currStart+i);
       else 
	 myDeriv=VarDeriv.ComputeDerivp(currStart+i);
       if (myDeriv.real()!=0)
	 cerr<<"My parameter deriv is "<<i<< " of "<<Psi.NumParams<<" "<<myDeriv<<" "<<(myDeriv.real())/abs(myDeriv.real())<<" "<<Psi.GetParam_real(i)<<endl;
       else 
	 cerr<<"My parameter deriv is "<<i<< " of "<<Psi.NumParams<<" "<<myDeriv<<endl;
       //	 cerr<<"My parameter deriv is "<<myDeriv<<" "<<(myDeriv.real())/abs(myDeriv.real())<<" "<<Psi.GetParam_real(i)<<endl;

       if (opt==GRADIENT || opt==TIMEEVOLUTION){
	 Psi.SetParam_real(i,Psi.GetParam_real(i)+-StepSize*(myDeriv.real())); //Han-Yi Chou
	 //	 if (i==0)
	 //	   Psi.SetParam_real(i, Psi.GetParam_real(i)+0.2);
       }
       else if (opt==SR && myDeriv.real()!=0){
	 double newParam=Psi.GetParam_real(i)+-StepSize*(myDeriv.real())/abs(myDeriv.real())*Random.ranf();
	 Psi.SetParam_real(i,newParam); //Psi.GetParam_real(i)+-StepSize*(myDeriv.real())/abs(myDeriv.real())*Random.ranf());
       }
       else if (myDeriv.real()==0){
	 double newParam=Psi.GetParam_real(i)+-StepSize*0.1*Random.ranf()* (Random.ranf()>0.5 ? -1: 1);
	 Psi.SetParam_real(i,newParam); //Psi.GetParam_real(i)+-StepSize*(myDeriv.real())/abs(myDeriv.real())*Random.ranf());

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
   for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
     WaveFunctionClass &Psi =**wf_iter;
     Psi.RebuildParams();
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


  

  void Optimize()
  {
    cerr<<"IN optimize"<<endl;
    int currStart=0;
    for (list<WaveFunctionClass*>::iterator wf_iter=wf_list.begin();wf_iter!=wf_list.end();wf_iter++){
      WaveFunctionClass &Psi =**wf_iter;
      for (int i=0;i<Psi.NumParams;i++){
	ParamsOld(currStart+i).real(Psi.GetParam_real(i));
	ParamsOld(currStart+i).imag(Psi.GetParam_imag(i));
      }
      currStart=currStart+Psi.NumParams;
    }
        
    //cout << "I am here" << endl; 
    VarDeriv.Clear();
    derivs=0;
    
    ///FIXME!
    
    double energy=0.0;
    int NumCounts=0;
    for (int sweeps=0;sweeps<opt_equilSweeps;sweeps++)
      Sweep_hop();
    for (int sweeps=0;sweeps<opt_SampleSweeps;sweeps++){
      Sweep_hop();
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


  double VMC(bool equilibrate)
  {
    cerr<<"Started VMC"<<endl;
    double numAttempt=0;
    double numAccept=0;
    if (equilibrate){
      for (int sweeps=0;sweeps<VMC_equilSweeps;sweeps++){
	//	cerr<<"Doing step "<<sweeps<<endl;
	Sweep_hop();
      }
    }
    vector<double> energy_terms(Ham.size(),0);
    double energy=0.0;
    int NumCounts=0;
    for (int sweeps=0;sweeps<VMC_SampleSweeps;sweeps++){
      //      cerr<<"Doing step "<<sweeps<<endl;
      numAccept+=Sweep_hop();
      numAttempt+=1;
      NumCounts++;
      int i=0;
      for (list<HamiltonianClass*>::iterator iter=Ham.begin();iter!=Ham.end();iter++){
	//       	cerr<<"J1 J2: "<<(*iter)->term1<<" "<<(*iter)->term2<<endl;
	//	(*iter)->term1=0; (*iter)->term2=0;
	double te=(*iter)->Energy(System,wf_list);
	energy_terms[i]+=te;
	i++;
	energy+=te;

      }

    }
    int NumParams=(*(wf_list.begin()))->NumParams;

    //    cerr<<(*(wf_list.begin()))->NumParams<<endl;
    //    Array<complex<double>,1> derivs(NumParams);
    //    (*(wf_list.begin()))->AllDerivs(System,derivs,0,derivs.size());
    //    (*(wf_list.begin()))->CheckDerivs(System,derivs,0,derivs.size());

    double oldEnergy=energy/(double)NumCounts;
    if ((1==1)) { 
      cerr<<"VMCENERGY: "<<energy/(double)NumCounts<<" "<<endl;
      cerr<<"J1 J2: "<<((*(Ham.begin()))->term1)/(double)NumCounts<<" "<<((*(Ham.begin()))->term2)/(double)NumCounts<<endl;
      (*(Ham.begin()))->term1=0; (*(Ham.begin()))->term2=0;

      cerr<<"ENERGY TERMS: ";
      for (int i=0;i<energy_terms.size();i++)
	cerr<<energy_terms[i]/(double)NumCounts<<" ";
      cerr<<endl;

      cerr<<"ACCEPT: "<<numAccept/numAttempt<<endl;
      energy=0.0;
      NumCounts=0;
    }
    cerr<<"Ended VMC"<<endl;
    return oldEnergy;
  }
  

 





};

//Old CPS initializiation
/*     InitDeriv(); */
/*     OptimizeCalled=0; */
/*     H.Init(System); */
/*     Hp.Init(System,J2); */
/*     if (1==2){ */
/*      for (int i=0;i<System.x.size();i++)  */
/*        for (int j=i+1;j<System.x.size();j++){  */
/*  	set<int> theSet;  */
/*  	std::set_intersection( Psi_CSP.PF.correlatorsForSite[i].begin(),  Psi_CSP.PF.correlatorsForSite[i].end(),  */
/*  			       Psi_CSP.PF.correlatorsForSite[j].begin(),  Psi_CSP.PF.correlatorsForSite[j].end(),   */
/*  			       std::inserter( theSet, theSet.begin() ) );  */
/*  	int theBin=Psi.PairingFunction.FindBin(i,j);  */
/*  	int CSP_bin=(*(theSet.begin()));  */
/*  	if (Psi.PairingFunction.f0[theBin].real()==0 || Psi.PairingFunction_diff.f0[theBin].real()==0)  */
/*  	  assert(1==2);  */
/*  	Psi_CSP.PF.f0[CSP_bin][0]=Psi.PairingFunction.f0[theBin];  */
/*  	Psi_CSP.PF.f0[CSP_bin][3]=Psi.PairingFunction.f0[theBin];  */
/*  	Psi_CSP.PF.f0[CSP_bin][1]=Psi.PairingFunction_diff.f0[theBin];  */
/*  	Psi_CSP.PF.f0[CSP_bin][2]=Psi.PairingFunction_diff.f0[theBin];  */
/*        }  */
/*      for (int i=0;i<Psi_CSP.PF.f0.size();i++)  */
/*        if (Psi_CSP.PF.f0[i][0].real()==0 || Psi_CSP.PF.f0[i][1].real()==0)  */
/*  	assert(1==2);  */
/*     } */

#endif
