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
class VMCDriverClass
{
 public:



  void RunSingleVMC_old(ifstream &infile)
  {

    std::map<string,string> input;

    ReadInput(infile,input);

    CommunicatorClass myComm;
    cerr<<"My processor is "<<myComm.MyProc()<<endl;
    RandomClass Random(myComm);
    Random.Init();
    list<pair<string,SharedWaveFunctionDataClass* > > wf_list;
    //need to delete these eventually if we don't want memory to leak
    wf_list.push_back(make_pair("RVB",new PairingFunctionAllBin()));
    //    wf_list.push_back(make_pair("CPS",new PairingFunctionMany()));
    OptimizeBothClass VMC(Random);
    VMC.Init(wf_list,input);
    VMC.VMC_equilSweeps=atoi(input["EquilSweeps"].c_str());
    VMC.VMC_SampleSweeps=atoi(input["SampleSweeps"].c_str());

    cerr<<"POST INIT"<<endl;
    VMC.EvaluateAll();
    cerr<<"POST evaluate all"<<endl;
    VMC.SaveParams("params.dat");
    cerr<<"C"<<endl;
    int step=0;
    VMC.VMC(true);
    while (1==1){
      step++;
      VMC.VMC(false);
      
      //      if (step % 5==0){
      //	((Heisenberg*)(*VMC.Ham.begin()))->PrintHistogram(VMC.System);
      //	((Heisenberg*)(*VMC.Ham.begin()))->ClearHistogram();
      //      }
    }
  }


  void RunSingleVMC(InputClass &myInput)
  {

    //    std::map<string,string> input;

    //    ReadInput(infile,input);

    CommunicatorClass myComm;
    cerr<<"My processor is "<<myComm.MyProc()<<endl;
    std::ostringstream filename;
    string outFileBase=myInput.GetVariable("OutFileBase");
    filename<<outFileBase<<"energy."<<myComm.MyProc();
    //    unique_ptr<ofstream> outfile(new ofstream);
    ofstream *outfile = new ofstream();
    outfile->open(filename.str().c_str());

    RandomClass Random(myComm);
    Random.Init();
    list<pair<string,SharedWaveFunctionDataClass* > > wf_list;
    //need to delete these eventually if we don't want memory to leak
    wf_list.push_back(make_pair("RVB",new PairingFunctionAllBin()));
    //    wf_list.push_back(make_pair("CPS",new PairingFunctionMany()));
    OptimizeBothClass VMC(Random);
    VMC.Init(wf_list,myInput);
    VMC.VMC_equilSweeps=myInput.toInteger(myInput.GetVariable("EquilSweeps"));
    VMC.VMC_SampleSweeps=myInput.toInteger(myInput.GetVariable("SampleSweeps"));

    //    cerr<<"POST INIT"<<endl;
    VMC.EvaluateAll();
    //    cerr<<"POST evaluate all"<<endl;
    if (myComm.MyProc()==0)
      VMC.SaveParams("params.dat");
    //    cerr<<"C"<<endl;
    int step=0;
    VMC.VMC(true,outfile);
    while (1==1){
    //    for (int aa=0;aa<100;aa++){
      step++;
      VMC.VMC(false,outfile);
      
      //      if (step % 5==0){
      //	((Heisenberg*)(*VMC.Ham.begin()))->PrintHistogram(VMC.System);
      //	((Heisenberg*)(*VMC.Ham.begin()))->ClearHistogram();
      //      }
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
    //wf_list.push_back(make_pair("RVB",new PairingFunctionAllBin()));
    wf_list.push_back(make_pair("CPS",new PairingFunctionMany()));
    
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



  void RunFiniteT(InputClass &input)
  {
    //    std::map<string,string> input;
    //    ReadInput(infile,input);
    CommunicatorClass myComm; 
    cerr<<"My processor is "<<myComm.MyProc()<<endl; 
    RandomClass Random(myComm); 
    Random.Init(); 

    vector<OptimizeBothClass*> VMC_vec; 
    OptimizeBothClass VMC_combine(Random); 

    
    list<pair<string,SharedWaveFunctionDataClass* > > wf_list;
    //    wf_list.push_back(make_pair("RVB",new PairingFunctionAllBin()));
    wf_list.push_back(make_pair("CPS",new PairingFunctionMany()));
    VMC_combine.Init(wf_list,input);

    //    VMC_combine.opt_equilSweeps=atoi(input["EquilSweeps"].c_str());
    //    VMC_combine.opt_SampleSweeps=atoi(input["SampleSweeps"].c_str()); 
    VMC_combine.opt_equilSweeps=input.toInteger(input.GetVariable("EquilSweeps"));  //atoi(input["EquilSweeps"].c_str());
    VMC_combine.opt_SampleSweeps=input.toInteger(input.GetVariable("SampleSweeps")); //atoi(input["SampleSweeps"].c_str()); 

    int NumWalkers=input.toInteger(input.GetVariable("NumWalkers")); // atoi(input["NumWalkers"].c_str()); 
    cerr<<"Number of Walkers: "<<NumWalkers<<endl;
    cerr<<"Equil Sweeps: "<<VMC_combine.opt_equilSweeps<<endl;
    cerr<<"Sample Sweeps: "<<VMC_combine.opt_SampleSweeps<<endl;
    VMC_vec.resize(NumWalkers); 

    cerr<<"Running VMC"<<endl; 
    //Initializing must be done in serial since
    //everyone is setting up the pairing function
    for (int i=0;i<NumWalkers;i++){ 
       VMC_vec[i]=new OptimizeBothClass(Random); 
       VMC_vec[i]->Init(wf_list,input); 
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


     int max_markovSteps= input.toInteger(input.GetVariable("markovSteps"));
     int numSteps=input.toInteger(input.GetVariable("TimeStepsTaken"));
     string outFileBase=input.GetVariable("OutFileBase");
     vector<ofstream*> myFiles(NumWalkers);
     for (int i=0;i<NumWalkers;i++){
       std::ostringstream filename;
       filename<<outFileBase<<"observables."<<myComm.MyProc()<<"_"<<i;
       myFiles[i]= new ofstream();
       myFiles[i]->open(filename.str().c_str());
     }
     ///INITIALIZATION UP TO HERE SHOULD BE SAME FOR FINITE TEMPERATURE!
     for (int markovSteps=0;markovSteps<max_markovSteps;markovSteps++){
       //     int numSteps=10;
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
     //Here is where you want to compute the observables

     ///This will be slow when you have many walkers per node!
     //BUG: Slow 
     for (int i=0;i<NumWalkers;i++)  {
       VMC_vec[i]->computeObservables=true;
       cout <<"computingObservables"<<endl;
       VMC_vec[i]->VMC(true,myFiles[i]);
       VMC_vec[i]->computeObservables=false;
     }
	 
     VMC_vec[0]->VMC(true); 
     cerr<<"The particles are at ";
     for (int i=0;i<VMC_vec[0]->System.x.size();i++)
       cerr<<VMC_vec[0]->System.x(i)<<" ";
     cerr<<endl;
       //       cerr<<"The "<<i<<" particle is at "<<VMC_vec[0]->System.x(i)<<endl;
     cerr<<"My staggered magnetization is "<<MeasureStaggered(*VMC_vec[0])<<endl;
     CPSClass *cps =((CPSClass*)(*(VMC_vec[0]->wf_list.begin())));
     for (int i=0;i<cps->PF.f0.size();i++){
       int theBin=cps->PF.corr2Bin(VMC_vec[0]->System.x,i);
       for (int j=0;j<cps->PF.f0[i].size();j++){
	 if (theBin==j || i>=4)
	   cps->PF.f0[i][j]=1+0.001*i+0.002*j;
	 else
       	   cps->PF.f0[i][j]=0.001+0.001*i+0.002*j;
       }
     }
     for (int i=0;i<cps->PF.f0.size();i++){
       for (int j=0;j<cps->PF.f0[i].size();j++){
	 cerr<<i<<" "<<j<<" "<<cps->PF.f0[i][j]<<endl;
       }
     }
     
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


  void RunFiniteT_exact(ifstream &infile)
  {
    Hamiltonian2 H;
    H.Init();
    H.GenerateBasis();

    std::map<string,string> input;
    ReadInput(infile,input);
    CommunicatorClass myComm; 
    cerr<<"My processor is "<<myComm.MyProc()<<endl;  
    RandomClass Random(myComm);  
    Random.Init();  
    
    OptimizeBothClass VMC(Random); 
    list<pair<string,SharedWaveFunctionDataClass* > > wf_list;
/*     //    wf_list.push_back(make_pair("RVB",new PairingFunctionAllBin())); */
     wf_list.push_back(make_pair("CPS",new PairingFunctionMany())); 
     
     VMC.Init(wf_list,input);
     VMC.opt_equilSweeps=atoi(input["EquilSweeps"].c_str()); 
     VMC.opt_SampleSweeps=atoi(input["SampleSweeps"].c_str());  
     
     int NumWalkers=atoi(input["NumWalkers"].c_str());  
     
     cerr<<"Number of Walkers: "<<NumWalkers<<endl; 
     cerr<<"Equil Sweeps: "<<VMC.opt_equilSweeps<<endl; 
     cerr<<"Sample Sweeps: "<<VMC.opt_SampleSweeps<<endl; 
     
     VMC.EvaluateAll(); 
     VMC.VMC(true); 
     
     CPSClass *cps =((CPSClass*)(*(VMC.wf_list.begin()))); 


     vector<double> Psi(H.basis.size());


     Array<complex<double>,1 > derivs;
     derivs.resize(VMC.NumberOfParams);


     vector<vector<double> > PsiDerivs;
     PsiDerivs.resize(derivs.size());
     for (int i=0;i<PsiDerivs.size();i++)
       PsiDerivs[i].resize(H.basis.size());



     for (int step=0;step<1000;step++){
       VMC.EvaluateAll(); 
       VMC.VMC(true); 

       for (int i=0;i<H.basis.size();i++){
	 for (int j=0;j<VMC.System.x.size();j++)
	   VMC.System.x(j)=H.basis[i][j];
	 Psi[i]=(cps->evaluate(VMC.System)).real();
	 cps->RealDerivs(VMC.System,derivs);
	 
	 for (int d=0;d<derivs.size();d++){
	   PsiDerivs[d][i]=derivs(d).real();
	 }
       }
       
       Eigen::MatrixXd S;
       Eigen::MatrixXd Hp;
       S.resize(derivs.size(),derivs.size());
       for (int i=0;i<S.cols();i++){
	 for (int j=0;j<S.rows();j++){
	   S(i,j)=0.0;
	   for (int k=0;k<PsiDerivs[i].size();k++){
	     S(i,j)+=PsiDerivs[i][k]*PsiDerivs[j][k];
	   }
	 }
       }
       vector<double> scratch(H.basis.size()); 
       Hp.resize(derivs.size(),derivs.size());     
       for (int i=0;i<S.cols();i++){
	 H.mult(PsiDerivs[i],scratch);
	 for (int j=0;j<S.rows();j++){
	   Hp(i,j)=0.0;
	   for (int k=0;k<PsiDerivs[i].size();k++){
	     Hp(i,j)+=scratch[k]*PsiDerivs[j][k];
	   }
	 }
       }

       for (int i=0;i<S.cols();i++)
	   S(i,i)+=1e-6;

       //       cerr<<"S is "<<S<<endl;
       //       cerr<<"Hp is "<<Hp<<endl;
       
       double tau=0.01;
       Eigen::VectorXd psi_v;
       psi_v.resize(derivs.size());
       for (int i=0;i<derivs.size();i++){
	 psi_v[i]=cps->GetParam_real(i);
       }
       //       cerr<<"GRR: "<<(Eigen::MatrixXd::Identity(Hp.rows(),Hp.cols())-tau*(S.inverse()*Hp))<<endl;
       Eigen::MatrixXd newVec=(Eigen::MatrixXd::Identity(Hp.rows(),Hp.cols())-tau*(S.inverse()*Hp))*psi_v;
       for (int i=0;i<newVec.size();i++){
	 cps->SetParam_real(i,newVec(i));
       }
       for (int i=0;i<newVec.size();i++)
	 cerr<<"The current parameters are "<<newVec(i)<<endl;
     }
       
       //     cerr<<"My new vector is "<<(Eigen::MatrixXd::Identity(Hp.rows(),Hp.cols())-tau*(S.inverse()*Hp))*psi_v<<endl;
       



       //     cerr<<"To Multiply by "<<Eigen::MatrixXd::Identity(Hp.rows(),Hp.cols())-tau*(S.inverse()*Hp)<<endl;
     //     cerr<<"To Multiply by "<<Eigen::MatrixXd::Identity(Hp.rows(),Hp.cols())-tau*S.inverse()*Hp<<endl;
     //     for (int i=0;i<H.basis.size();i++){
     //       cerr<<H.basis[i]<<endl;
     //     }
     exit(1);


  }



  void RunFiniteT_exact2(ifstream &infile)
  {
    cerr<<"A"<<endl;

    Hamiltonian2 H;
    cerr<<"B"<<endl;
    H.Init();
    cerr<<"C"<<endl;
    H.GenerateBasis();
    cerr<<"D"<<endl;

   std::map<string,string> input;
   ReadInput(infile,input);
   CommunicatorClass myComm; 
     cerr<<"My processor is "<<myComm.MyProc()<<endl;  
     RandomClass Random(myComm);  
     Random.Init();  

     OptimizeBothClass VMC(Random); 
     list<pair<string,SharedWaveFunctionDataClass* > > wf_list;
/*     //    wf_list.push_back(make_pair("RVB",new PairingFunctionAllBin())); */
     wf_list.push_back(make_pair("CPS",new PairingFunctionMany())); 

     VMC.Init(wf_list,input);
     VMC.opt_equilSweeps=atoi(input["EquilSweeps"].c_str()); 
     VMC.opt_SampleSweeps=atoi(input["SampleSweeps"].c_str());  

     int NumWalkers=atoi(input["NumWalkers"].c_str());  

     cerr<<"Number of Walkers: "<<NumWalkers<<endl; 
     cerr<<"Equil Sweeps: "<<VMC.opt_equilSweeps<<endl; 
     cerr<<"Sample Sweeps: "<<VMC.opt_SampleSweeps<<endl; 

     VMC.EvaluateAll(); 
     VMC.VMC(true); 
     
     CPSClass *cps =((CPSClass*)(*(VMC.wf_list.begin()))); 


     vector<double> Psi(H.basis.size());


     Array<complex<double>,1 > derivs;
     derivs.resize(VMC.NumberOfParams);


     vector<vector<double> > PsiDerivs;
     PsiDerivs.resize(derivs.size());
     for (int i=0;i<PsiDerivs.size();i++)
       PsiDerivs[i].resize(H.basis.size());

     for (int i=0;i<H.basis.size();i++){
       for (int j=0;j<VMC.System.x.size();j++)
	 VMC.System.x(j)=H.basis[i][j];
       Psi[i]=(cps->evaluate(VMC.System)).real();
       cps->RealDerivs(VMC.System,derivs);

/*         cerr<<"My derivs is "; */
/*         for (int ii=0;ii<derivs.size();ii++){  */
/*  	 cerr<<derivs(ii)<<" "<<" ";  */
/*         }  */
/*         cerr<<endl;  */


       for (int d=0;d<derivs.size();d++){
	 PsiDerivs[d][i]=derivs(d).real();
       }
     }
     Eigen::MatrixXd S;
     Eigen::MatrixXd Hp;
     S.resize(derivs.size(),derivs.size());
     for (int i=0;i<S.cols();i++){
       for (int j=0;j<S.rows();j++){
	 S(i,j)=0.0;
	 for (int k=0;k<PsiDerivs[i].size();k++){
	   S(i,j)+=PsiDerivs[i][k]*PsiDerivs[j][k];
	 }
       }
     }
     cerr<<"S is "<<endl<<S<<endl;

     cerr<<"S inverse is "<<endl<<S.inverse()<<endl;

     vector<double> scratch(H.basis.size()); 
     Hp.resize(derivs.size(),derivs.size());     
     for (int i=0;i<S.cols();i++){
       H.mult(PsiDerivs[i],scratch);
       for (int j=0;j<S.rows();j++){
	 Hp(i,j)=0.0;
	 for (int k=0;k<PsiDerivs[i].size();k++){
	   Hp(i,j)+=scratch[k]*PsiDerivs[j][k];
	 }
       }
     }

     cerr<<"Hp is "<<endl;
     cerr<<Hp<<endl;

     double tau=0.01;
     Eigen::VectorXd psi_v;
     psi_v.resize(derivs.size());
     for (int i=0;i<derivs.size();i++){
       psi_v[i]=cps->GetParam_real(i);
     }
     cerr<<"My new vector is "<<(Eigen::MatrixXd::Identity(Hp.rows(),Hp.cols())-tau*(S.inverse()*Hp))*psi_v<<endl;
       //     cerr<<"To Multiply by "<<Eigen::MatrixXd::Identity(Hp.rows(),Hp.cols())-tau*(S.inverse()*Hp)<<endl;
     //     cerr<<"To Multiply by "<<Eigen::MatrixXd::Identity(Hp.rows(),Hp.cols())-tau*S.inverse()*Hp<<endl;
     //     for (int i=0;i<H.basis.size();i++){
     //       cerr<<H.basis[i]<<endl;
     //     }
     exit(1);
     
       /*      ofstream energyFile; */
/*      energyFile.open("Energy.dat"); */

/*      ///INITIALIZATION UP TO HERE SHOULD BE SAME FOR FINITE TEMPERATURE! */
/*      for (int markovSteps=0;markovSteps<1000;markovSteps++){ */
/*      int numSteps=10; */
/*      for (int step=0;step<numSteps;step++){  */
/*        cerr<<"Step number: "<<step<<endl;  */
/* #pragma omp parallel for  */
/*        for (int i=0;i<NumWalkers;i++) {  */
/* 	 VMC_vec[i]->Optimize();  */
/*        }  */
/*        VMC_combine.Combine(VMC_vec); */
/*        VMC_combine.VarDeriv.ParallelCombine(myComm);  */
/*        if (myComm.MyProc()==0){  */
/* 	 VMC_combine.TakeStep(myComm);  */
/*        }  */
/*        VMC_combine.BroadcastParams(myComm);  */

/*        if (myComm.MyProc()==0){ */
/* 	 string fileName="params.dat."; */
/* 	 ostringstream ss; */
/* 	 ss<<fileName<<(step+1); */
/* 	 VMC_combine.SaveParams(ss.str());  */
/* 	 complex<double> theEnergy=VMC_combine.VarDeriv.ComputeEnergy(); */
/* 	 complex<double> theVariance=VMC_combine.VarDeriv.ComputeVariance(); */
/* 	 energyFile<<step<<" "<<theEnergy.real()<<" "<<theEnergy.imag()<<" "<<theVariance.real()<<" "<<theVariance.imag()<<endl; */
/*        } */
/*        VMC_combine.VarDeriv.Clear();  */
/*        for (int i=0;i<NumWalkers;i++)   */
/* 	 VMC_vec[i]->VarDeriv.Clear();  */
/* #pragma omp parallel for  */
/*        for (int i=0;i<NumWalkers;i++)   */
/* 	 VMC_vec[i]->EvaluateAll();  */
/*        //I don't think VMC_Combine needs to evaluate all but I'm a bit worried? */
/*      } */

/*      VMC_vec[0]->VMC(true);  */
/*      cerr<<"The particles are at "; */
/*      for (int i=0;i<VMC_vec[0]->System.x.size();i++) */
/*        cerr<<VMC_vec[0]->System.x(i)<<" "; */
/*      cerr<<endl; */
/*        //       cerr<<"The "<<i<<" particle is at "<<VMC_vec[0]->System.x(i)<<endl; */
/*      cerr<<"My staggered magnetization is "<<MeasureStaggered(*VMC_vec[0])<<endl; */
/*      CPSClass *cps =((CPSClass*)(*(VMC_vec[0]->wf_list.begin()))); */
/*      for (int i=0;i<cps->PF.f0.size();i++){ */
/*        int theBin=cps->PF.corr2Bin(VMC_vec[0]->System.x,i); */
/*        for (int j=0;j<cps->PF.f0[i].size();j++){ */
/* 	 if (theBin==j || i>=4) */
/* 	   cps->PF.f0[i][j]=1+0.001*i+0.002*j; */
/* 	 else */
/*        	   cps->PF.f0[i][j]=0.001+0.001*i+0.002*j; */
/*        } */
/*        //       if (cps->PF.f0[i].size()==2){ */
/*        //	 int spin=VMC_vec[0]->System.x(i); */
/* 	 //	 cps->PF.f0[i][1]=1*(spin)+0.001*(i+1)+0.002*0; */
/* 	 //	 cps->PF.f0[i][0]= -1*(spin-1)+0.001*(i+1)+0.002*1; */
/* /\* 	 if (spin==0){ *\/ */
/* /\* 	   cps->PF.f0[i][1]=1+0.001*(i+1)+0.002*0; *\/ */
/* /\* 	   cps->PF.f0[i][0]=1+0.001*(i+1)+0.002*1; *\/ */
/* /\* 	 } *\/ */
/* /\* 	 else { *\/ */
/* /\* 	   cps->PF.f0[i][1]=1+0.001*(i+1)+0.002*0; *\/ */
/* /\* 	   cps->PF.f0[i][0]=-1+0.001*(i+1)+0.002*1; *\/ */
/* /\* 	 } *\/ */
/* //	 cps->PF.f0[i][1]= 1*(spin)+0.001*(i+1)+0.002*0+0.01; */
/* //	 cps->PF.f0[i][0]= -1*(spin-1)+0.001*(i+1)+0.002*1+0.01; */
/* //       } */
/* //       else { */
/* //	 for (int j=0;j<cps->PF.f0[i].size();j++){ */
/* //	   cps->PF.f0[i][j]=1+0.001*i+0.002*j; */
/* //	 } */
/* //       } */
/*      } */
/*      for (int i=0;i<cps->PF.f0.size();i++){ */
/*        for (int j=0;j<cps->PF.f0[i].size();j++){ */
/* 	 cerr<<i<<" "<<j<<" "<<cps->PF.f0[i][j]<<endl; */
/*        } */
/*      } */
     
/*      VMC_vec[0]->BroadcastParams(myComm);  */
/*      ///WONt BROADCAST RIGHT!!!! BUG: */
/*      for (int i=0;i<NumWalkers;i++) { */
/*        for (int j=0;j<VMC_vec[0]->System.x.size();j++) */
/*        	 VMC_vec[i]->System.x(j)=VMC_vec[0]->System.x(j); */
/*        VMC_vec[i]->EvaluateAll();  */
/*        VMC_vec[i]->VMC(true);  */
/*      } */
/*      //I don't think VMC_Combine needs to evaluate all but I'm a bit worried? */
/*      } */
     
/*      energyFile.close();  */
/*      exit(1); */


  }

  
  //  void RunMultipleOpt(ifstream &infile)
  void RunMultipleOpt(InputClass &myInput)
  {
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
    string waveFunction=myInput.GetVariable("WaveFunction");
    if (waveFunction=="CPS")
      wf_list.push_back(make_pair("CPS",new PairingFunctionMany()));
    else if (waveFunction=="RVB")
      wf_list.push_back(make_pair("RVB",new PairingFunctionAllBin()));
    else 
      assert(1==2);

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
       cerr<<"Outside "<<endl;
       if (myComm.MyProc()==0){ 
	 VMC_combine.TakeStep(myComm); 
       } 
       TakeStepTimer.Stop();
       BroadcastStepTimer.Start();
       VMC_combine.BroadcastParams(myComm); 
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
