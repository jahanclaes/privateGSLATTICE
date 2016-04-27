#ifndef VMC_H
#define VMC_H

#include <iostream>
#include "Communication/Communication.h"
#include "Random/Random.h"
#include "OptimizeBothClass.h"
#include <omp.h>
#include <bitset>
#include <map>
#include <utility>
#include <vector>
#include "Hamiltonian.h"
#include "exact.h"
#include "Timer.h"
#include "SharedEigs.h"
#include "kondoHelp.h"
class VMCDriverClass
{
 public:



  void RunSingleVMC(InputClass &myInput)
  {
    CommunicatorClass myComm;
    cerr<<"My processor is "<<myComm.MyProc()<<endl;
    RandomClass Random(myComm);
    Random.Init();
    list<pair<string,SharedWaveFunctionDataClass* > > wf_list;


    ReadWaveFunction(myInput,wf_list);


    OptimizeBothClass VMC(Random);
    VMC.Init(wf_list,myInput);
    VMC.VMC_equilSweeps=myInput.toInteger(myInput.GetVariable("EquilSweeps"));
    VMC.VMC_SampleSweeps=myInput.toInteger(myInput.GetVariable("SampleSweeps"));

    /////    VMC.GetParams("params.dat"); 
    string mySweep=myInput.GetVariable("MoveType");
    if (mySweep=="HOP")
      VMC.sweep=HOP;
    else if (mySweep=="EXCHANGE")
      VMC.sweep=EXCHANGE;
    ///GETTING PARAMETERS
    if (myInput.IsVariable("ReadParams")){
      bool toRead=myInput.toBool(myInput.GetVariable("ReadParams"));
      if (toRead){
	string param_fileName=myInput.GetVariable("ParamFileName");
	VMC.GetParams(param_fileName); 
      }
    }
    //DONE GETTING PARAMETERS



    VMC.EvaluateAll();

    //    VMC.SaveParams("params.dat");
<<<<<<< HEAD

=======
>>>>>>> master
    int step=0;
    VMC.VMC(true);
    while (1==1){
      step++;
      VMC.VMC(false);
    }
  }

  void RunSingleOpt(ifstream &infile)
  {
    std::map<string,string> input;
    ReadInput(infile,input);

    CommunicatorClass myComm;
    cerr<<"My processor is "<<myComm.MyProc()<<endl;
    RandomClass Random(myComm);
    Random.Init();
    list<pair<string,SharedWaveFunctionDataClass* > > wf_list;
    //need to delete these eventually if we don't want memory to leak
    //    wf_list.push_back(make_pair("RVB",new PairingFunctionAllBin()));
    //wf_list.push_back(make_pair("CPS",new PairingFunctionMany()));
    wf_list.push_back(make_pair("BACKFLOW",new PairingFunctionMany()));
    
    OptimizeBothClass VMC(Random);
    //Must send the input so that it can initilalize Hamiltonians
    //and things!
    VMC.Init(wf_list,input);
    VMC.VMC_equilSweeps=atoi(input["EquilSweeps"].c_str());
    VMC.VMC_SampleSweeps=atoi(input["SampleSweeps"].c_str());



    VMC.EvaluateAll();


    VMC.SaveParams("params.dat");
    int step=0;

    for  (int i=0;i<100;i++){
      step++;
      VMC.VMC(true);
    }

    //If you want to optimize:    
    for (int step=0;step<10000;step++){
      VMC.Optimize();
      VMC.TakeStep(myComm);
      cerr<<"ACENERGY: "<<VMC.VarDeriv.ComputeEnergy()<<endl;;
      VMC.VarDeriv.Clear();
      VMC.EvaluateAll();
      VMC.SaveParams("params.dat");
    }
  }
  
  

  void ReadInput(ifstream &infile,std::map<string,string> &input)
  {
    //    std::map<string,string> input;
    while (!infile.eof()){
      string line;
      infile>>line;
      int eqLoc=line.find('=');
      string token=line.substr(0,eqLoc);
      string outVals=line.substr(eqLoc+1);
      input[token]=outVals;
      cerr<<"INPUT: "<<token<<" "<<input[token]<<endl;
    }
    
  }


double MeasureStaggered(OptimizeBothClass &vmc)
{
  
  vector<int> s(16);
  int currSign=1;
  for (int i=0;i<4;i++){
    for (int j=0;j<4;j++){
      s[4*i+j]=currSign;
      currSign*=-1;
    }
    currSign*=-1;
  }
  double sl_local=0.0;
  for (int i=0;i<vmc.System.x.size();i++)
    sl_local+=s[(i % 16)]*(vmc.System.x(i)==0 ? 1: -1);
  return sl_local*sl_local; 
}

 void ReadWaveFunction(InputClass &myInput,
		       list<pair<string,SharedWaveFunctionDataClass* > > &wf_list)
 {
   int check=0;
   while (myInput.OpenSection("WaveFunction",check)){
     string waveFunction=myInput.GetVariable("name");
     cerr<<"Found the "<<check<<" wavefunction called"<<waveFunction<<endl;
     if (waveFunction=="PEPS"){
       wf_list.push_back(make_pair("PEPS",new PairingFunctionMany()));
       
     }
     else if  (waveFunction=="CPS"){
       wf_list.push_back(make_pair("CPS",new PairingFunctionMany()));	
     }
     else if (waveFunction=="RVB"){
       wf_list.push_back(make_pair("RVB",new PairingFunctionAllBin()));
     }
     else if (waveFunction=="BACKFLOW")
       wf_list.push_back(make_pair("BACKFLOW",new PairingFunctionAllBin()));
     else if (waveFunction=="JASTROW"){
       wf_list.push_back(make_pair("JASTROW",new PairingFunctionAllBin()));
     }
     else if (waveFunction=="SLATERDET"){
       wf_list.push_back(make_pair("SLATERDETUP",new SharedEigsClass()));
       wf_list.push_back(make_pair("SLATERDETDOWN",new SharedEigsClass()));
     }
     else{
       cerr<<"Could not find wavefunction: "<<waveFunction <<endl;
       exit(1);
     }
     myInput.CloseSection();
     check++;
   }

 }

 void RunFiniteT_new(InputClass &myInput)
 {
   CommunicatorClass myComm; 
   cerr<<"My processor is "<<myComm.MyProc()<<endl; 
   RandomClass Random(myComm); 
   Random.Init(); 
   vector<OptimizeBothClass*> VMC_vec; 
   OptimizeBothClass VMC_combine(Random); 
   list<pair<string,SharedWaveFunctionDataClass* > > wf_list;
   ReadWaveFunction(myInput,wf_list);

   VMC_combine.Init(wf_list,myInput);
   VMC_combine.opt_equilSweeps=myInput.toInteger(myInput.GetVariable("EquilSweeps"));
   VMC_combine.opt_SampleSweeps=myInput.toInteger(myInput.GetVariable("SampleSweeps"));
   int NumWalkers=myInput.toInteger(myInput.GetVariable("NumWalkers"));
   cerr<<"Number of Walkers: "<<NumWalkers<<endl;
   cerr<<"Equil Sweeps: "<<VMC_combine.opt_equilSweeps<<endl;
   cerr<<"Sample Sweeps: "<<VMC_combine.opt_SampleSweeps<<endl;
   VMC_vec.resize(NumWalkers); 
   cerr<<"Running VMC"<<endl; 


   //Initializing must be done in serial since
   //everyone is setting up the pairing function
   for (int i=0;i<NumWalkers;i++){ 
     VMC_vec[i]=new OptimizeBothClass(Random); 
     VMC_vec[i]->Init(wf_list,myInput); 
     VMC_vec[i]->opt_equilSweeps=VMC_combine.opt_equilSweeps;
     VMC_vec[i]->opt_SampleSweeps=VMC_combine.opt_SampleSweeps;
   }
   //    VMC_combine.GetParams("params.dat"); 
   VMC_combine.EvaluateAll();
#pragma omp parallel for 
   for (int i=0;i<NumWalkers;i++){ 
     cerr<<"On walker "<<i<<endl; 
     VMC_vec[i]->EvaluateAll(); 
     VMC_vec[i]->VMC(true); 
   } 
   ofstream energyFile;
   energyFile.open("Energy.dat");
   
   ///INITIALIZATION UP TO HERE SHOULD BE SAME FOR FINITE TEMPERATURE!
   for (int markovSteps=0;markovSteps<1000;markovSteps++){
     int numSteps=10;
     for (int step=0;step<numSteps;step++){ 
       if (myComm.MyProc()==0)
	 cerr<<"Step number: "<<step<<endl; 
#pragma omp parallel for 
       for (int i=0;i<NumWalkers;i++) { 
	 VMC_vec[i]->Optimize(); 
       } 
       VMC_combine.Combine(VMC_vec);
       VMC_combine.VarDeriv.ParallelCombine(myComm); 
       if (myComm.MyProc()==0){ 
	 VMC_combine.TakeStep(myComm); 
       } 
       VMC_combine.BroadcastParams(myComm); 
       
       if (myComm.MyProc()==0){
	 string fileName="params.dat.";
	 ostringstream ss;
	 ss<<fileName<<(step+1);
	 VMC_combine.SaveParams(ss.str()); 
	 complex<double> theEnergy=VMC_combine.VarDeriv.ComputeEnergy();
	 complex<double> theVariance=VMC_combine.VarDeriv.ComputeVariance();
	 energyFile<<step<<" "<<theEnergy.real()<<" "<<theEnergy.imag()<<" "<<theVariance.real()<<" "<<theVariance.imag()<<endl;
       }
       VMC_combine.VarDeriv.Clear(); 
       for (int i=0;i<NumWalkers;i++)  
	 VMC_vec[i]->VarDeriv.Clear(); 
#pragma omp parallel for 
       for (int i=0;i<NumWalkers;i++)  
	 VMC_vec[i]->EvaluateAll(); 
       //I don't think VMC_Combine needs to evaluate all but I'm a bit worried?
     }

     VMC_vec[0]->VMC(true); 
     cerr<<"The particles are at ";
     for (int i=0;i<VMC_vec[0]->System.x.size();i++)
       cerr<<VMC_vec[0]->System.x(i)<<" ";
     cerr<<endl;

       //       cerr<<"The "<<i<<" particle is at "<<VMC_vec[0]->System.x(i)<<endl;
     //Measure something here!


     
     
     VMC_vec[0]->BroadcastParams(myComm); 
     ///WONt BROADCAST RIGHT!!!! BUG:
     for (int i=0;i<NumWalkers;i++) {
       for (int j=0;j<VMC_vec[0]->System.x.size();j++)
       	 VMC_vec[i]->System.x(j)=VMC_vec[0]->System.x(j);
       VMC_vec[i]->EvaluateAll(); 
       VMC_vec[i]->VMC(true); 
     }
     //I don't think VMC_Combine needs to evaluate all but I'm a bit worried?
     }
     
     energyFile.close(); 
     exit(1);


 }


  
  double RunMultipleVMC(vector<OptimizeBothClass*> &VMC_vec)
  {

    double energy=0.0;
    double count=0.0;
    double energy2=0.0;
#pragma omp parallel for 
    for (int i=0;i<VMC_vec.size();i++){ 
	VMC_vec[i]->EvaluateAll(); 
	VMC_vec[i]->VMC(true); 
      } 


#pragma omp parallel for reduction(+:energy,count,energy2)
      for (int i=0;i<VMC_vec.size();i++){ 
	  VMC_vec[i]->EvaluateAll(); 
	  count=count+1;
	  double myEnergy=VMC_vec[i]->VMC(true); 
	  energy+=myEnergy;
	  energy2+=myEnergy*myEnergy;
      }
      cerr<<"Variance: "<<(energy/count)*(energy/count)-(energy2/count)<<endl;
      cerr<<"std error: "<<((energy/count)*(energy/count)-(energy2/count))/sqrt(count)<<endl;
      return energy/count;
  }

  //  void RunMultipleOpt(ifstream &infile)
  void RunMultipleOpt(InputClass &myInput)
  {
    cerr<<"Running multiple opt "<<endl;
    //    std::map<string,string> input;
    //    ReadInput(infile,input);
    vector<TimerClass*> myTimers;
    TimerClass InitTimer("InitTimer");
    myTimers.push_back(&InitTimer);
    TimerClass OptimizeTimer("OptimizeTimer");
    myTimers.push_back(&OptimizeTimer);
    TimerClass CombineTimer("Combine");
    myTimers.push_back(&CombineTimer);
    TimerClass TakeStepTimer("TakeStep");
    myTimers.push_back(&TakeStepTimer);
    TimerClass BroadcastStepTimer("BroadcastStep");
    myTimers.push_back(&BroadcastStepTimer);

    InitTimer.Start();
    CommunicatorClass myComm; 
    cerr<<"My processor is "<<myComm.MyProc()<<endl; 
    RandomClass Random(myComm); 
    Random.Init(); 
    
    vector<OptimizeBothClass*> VMC_vec; 
    OptimizeBothClass VMC_combine(Random); 
    
    list<pair<string,SharedWaveFunctionDataClass* > > wf_list;
    ReadWaveFunction(myInput,wf_list);

    VMC_combine.Init(wf_list,myInput);
    
    VMC_combine.opt_equilSweeps=myInput.toInteger(myInput.GetVariable("EquilSweeps"));
    VMC_combine.opt_SampleSweeps=myInput.toInteger(myInput.GetVariable("SampleSweeps"));
    
    VMC_combine.VMC_equilSweeps=VMC_combine.opt_equilSweeps;
    VMC_combine.VMC_SampleSweeps=VMC_combine.opt_SampleSweeps;

    int NumWalkers=myInput.toInteger(myInput.GetVariable("NumWalkers"));
    cerr<<"Number of Walkers: "<<NumWalkers<<endl;
    cerr<<"Equil Sweeps: "<<VMC_combine.opt_equilSweeps<<endl;
    cerr<<"Sample Sweeps: "<<VMC_combine.opt_SampleSweeps<<endl;

    VMC_vec.resize(NumWalkers); 
    cerr<<"Running VMC"<<endl; 
    //Initializing must be done in serial since
    //everyone is setting up the pairing function
    for (int i=0;i<NumWalkers;i++){ 
       VMC_vec[i]=new OptimizeBothClass(Random); 
       VMC_vec[i]->Init(wf_list,myInput); 
       VMC_vec[i]->opt_equilSweeps=VMC_combine.opt_equilSweeps;
       VMC_vec[i]->opt_SampleSweeps=VMC_combine.opt_SampleSweeps;

       VMC_vec[i]->VMC_equilSweeps=VMC_combine.VMC_equilSweeps;
       VMC_vec[i]->VMC_SampleSweeps=VMC_combine.VMC_SampleSweeps;


    }


    ///GETTING PARAMETERS
    if (myInput.IsVariable("ReadParams")){
      bool toRead=myInput.toBool(myInput.GetVariable("ReadParams"));
      if (toRead){
	string param_fileName=myInput.GetVariable("ParamFileName");
	//	VMC_combine.GetParams("params.dat"); 
	VMC_combine.GetParams(param_fileName); 
	for (int i=0;i<NumWalkers;i++){ 
	  VMC_vec[i]->MatchParams(VMC_combine);
	}
      }
    }
    //DONE GETTING PARAMETERS

    //    (((cPEPSClass*)(*(VMC_combine.wf_list.begin())))->boostBD(1,1)); 
    //     VMC_combine.ResetParams(); 
    //     for (int i=0;i<NumWalkers;i++){ 
    //       (((cPEPSClass*)(*(VMC_vec[i]->wf_list.begin())))->boostBD(1,1)); 
    //       VMC_vec[i]->ResetParams(); 
    //     } 
     //     VMC_combine.SaveParams("params.dat.boosted");  
    VMC_combine.EvaluateAll();
#pragma omp parallel for 
     for (int i=0;i<NumWalkers;i++){ 
       cerr<<"On walker "<<i<<endl; 
       VMC_vec[i]->EvaluateAll(); 
       VMC_vec[i]->VMC(true); 
     } 

     ofstream energyFile;
     energyFile.open("Energy.dat");
     InitTimer.Stop();
     for (int step=0;step<10000*100;step++){ 
       if (myComm.MyProc()==0)
	 cerr<<"Step number: "<<step<<endl; 
       OptimizeTimer.Start();
#pragma omp parallel for 
       for (int i=0;i<NumWalkers;i++) { 
	 VMC_vec[i]->Optimize(); 
       } 
       OptimizeTimer.Stop();
       CombineTimer.Start();
       VMC_combine.Combine(VMC_vec);
       VMC_combine.VarDeriv.ParallelCombine(myComm); 
       CombineTimer.Stop();

       TakeStepTimer.Start();
       //       cerr<<"Outside "<<endl;
       if (myComm.MyProc()==0){ 
	 VMC_combine.TakeStep(myComm); 
       } 
       TakeStepTimer.Stop();
       BroadcastStepTimer.Start();
       VMC_combine.BroadcastParams(myComm); 
       for (int i=0;i<NumWalkers;i++){ 
	 VMC_vec[i]->MatchParams(VMC_combine);
       }
       if (myComm.MyProc()==0){
	 string fileName="params.dat.";
	 ostringstream ss;
	 ss<<fileName<<(step+1);
	 VMC_combine.SaveParams(ss.str()); 
	 //	 VMC_combine.SaveParams("params.dat"); 
	 cerr<<"ENERGY: "<<VMC_combine.VarDeriv.ComputeEnergy()<<endl;
	 complex<double> theEnergy=VMC_combine.VarDeriv.ComputeEnergy();
	 complex<double> theVariance=VMC_combine.VarDeriv.ComputeVariance();
	 energyFile<<step<<" "<<theEnergy.real()<<" "<<theEnergy.imag()<<" "<<theVariance.real()<<" "<<theVariance.imag()<<endl;
       }
       VMC_combine.VarDeriv.Clear(); 
       for (int i=0;i<NumWalkers;i++)  
	 VMC_vec[i]->VarDeriv.Clear(); 
#pragma omp parallel for 
       for (int i=0;i<NumWalkers;i++)  
	 VMC_vec[i]->EvaluateAll(); 
       BroadcastStepTimer.Stop();
       if (myComm.MyProc()==0){
       for (int i=0;i<myTimers.size();i++){
	 cerr<<"Timer: "<<myTimers[i]->Name<<" "<<myTimers[i]->Time()<<endl;
       }
       
       cerr<<"Step ENDING: "<<step<<endl; 
       }
     }
     
     energyFile.close();
  }
};
#endif
