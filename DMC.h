#ifndef DMC_H
#define DMC_H


#include <iostream>
#include "Communication/Communication.h"
#include "Random/Random.h"
#include "OptimizeBothClass.h"
#include "QMCTools.h"
#include <omp.h>
#include <bitset>

class DMCDriverClass
{
 public:
  void RunMultipleDMC(double& J2)
{

/*   CommunicatorClass myComm; */
/*   cerr<<"My processor is "<<myComm.MyProc()<<endl; */
/*   RandomClass Random(myComm); */
/*   Random.Init(); */

/*   PairingFunctionMany pf; */
/*   PairingFunctionAllBin pfUp; */
/*   PairingFunctionAllBin pfDown; */

/*   vector<OptimizeBothClass*> VMC_vec; */
/*   vector<OptimizeBothClass*> VMC_vec_branch; */

/*   int NumWalkers=20; */
/*   VMC_vec.resize(NumWalkers); */
/*   VMC_vec_branch.resize(NumWalkers); */
/*   vector<list<pair<SpinSwap,double> > > vals(NumWalkers); */
/*   vector<double> weight(NumWalkers); */
/*   vector<double> FullWeight(NumWalkers); */
/*   for (int i=0;i<FullWeight.size();i++) */
/*     FullWeight[i]=1.0; */
  

/*   for (int i=0;i<NumWalkers;i++){ */
/*     VMC_vec[i]=new OptimizeBothClass(Random); */
/*     VMC_vec[i]->SetJ2(J2); */
/*     VMC_vec[i]->Init(pf,pfUp,pfDown); */
/*     VMC_vec[i]->opt_equilSweeps=1000; */
/*     VMC_vec[i]->opt_SampleSweeps=10000; */
/*   } */

/*   for (int i=0;i<NumWalkers;i++){ */
/*     VMC_vec_branch[i]=new OptimizeBothClass(Random); */
/*     VMC_vec_branch[i]->SetJ2(J2); */
/*     VMC_vec_branch[i]->Init(pf,pfUp,pfDown); */
/*     VMC_vec_branch[i]->opt_equilSweeps=1000; */
/*     VMC_vec_branch[i]->opt_SampleSweeps=10000; */
/*   } */

/*   //  VMC_combine.GetParams("params.dat"); */
/*   //  VMC_combine.EvaluateAll(); */
  

/* #pragma omp parallel for */
/*   for (int i=0;i<NumWalkers;i++){ */
/*     cerr<<"On walker "<<i<<endl; */
/*     VMC_vec[i]->EvaluateAll(); */
/*     VMC_vec[i]->VMC(true); */
/*   } */


/* #pragma omp parallel for */
/*   for (int i=0;i<NumWalkers;i++){ */
/*     cerr<<"On walker "<<i<<endl; */
/*     VMC_vec_branch[i]->EvaluateAll(); */
/*     VMC_vec_branch[i]->VMC(true); */
/*   } */


/*   for (int DMCSteps=0;DMCSteps<1000;DMCSteps++){ */
/*     for (int i=0;i<NumWalkers;i++){ */
/*       for (list<HamiltonianClass*>::iterator ham_iter=VMC_vec[i]->Ham.begin(); */
/* 	   ham_iter!=VMC_vec[i]->Ham.end();ham_iter++){ */
/* 	(*ham_iter)->AllConnected(vals[i],VMC_vec[i]->System,VMC_vec[i]->wf_list); */
/*       } */
/*       int mySample=QMCTools::Sample(vals[i],Random,weight[i]); */
/*       cerr<<"The sample that I'm applying is "<<i<<" "<<mySample<<endl; */
/*       FullWeight[i]*=weight[i]; */
/*       VMC_vec[i]->Apply(vals[i],VMC_vec[i]->System,VMC_vec[i]->wf_list,mySample); */
      
/*     } */
    
/*     QMCTools::Branch(VMC_vec, */
/* 		     VMC_vec_branch, */
/* 		     FullWeight, */
/* 		     Random); */
    
/*     for (int i=0;i<FullWeight.size();i++) */
/*       FullWeight[i]=1.0; */
/*   } */


}


void RunMultipleDMC_CT(double& J2)
{

/*   CommunicatorClass myComm; */
/*   cerr<<"My processor is "<<myComm.MyProc()<<endl; */
/*   RandomClass Random(myComm); */
/*   Random.Init(); */

/*   PairingFunctionMany pf; */
/*   PairingFunctionAllBin pfUp; */
/*   PairingFunctionAllBin pfDown; */

/*   vector<OptimizeBothClass*> VMC_vec; */
/*   vector<OptimizeBothClass*> VMC_vec_branch; */

/*   int NumWalkers=80; */
/*   VMC_vec.resize(NumWalkers); */
/*   VMC_vec_branch.resize(NumWalkers); */
/*   vector<list<pair<SpinSwap,double> > > vals(NumWalkers); */
/*   vector<double> weight(NumWalkers); */
/*   vector<double> FullWeight(NumWalkers); */
/*   for (int i=0;i<FullWeight.size();i++) */
/*     FullWeight[i]=1.0; */
  

/*    for (int i=0;i<NumWalkers;i++){ */
/*      VMC_vec[i]=new OptimizeBothClass(Random); */
/*      VMC_vec[i]->SetJ2(J2); */
/*      VMC_vec[i]->Init(pf,pfUp,pfDown); */
/*      VMC_vec[i]->opt_equilSweeps=100; //0; */
/*      VMC_vec[i]->opt_SampleSweeps=1000; //0; */
/*    } */

/*   for (int i=0;i<NumWalkers;i++){ */
/*     VMC_vec_branch[i]=new OptimizeBothClass(Random); */
/*     VMC_vec_branch[i]->SetJ2(J2); */
/*     VMC_vec_branch[i]->Init(pf,pfUp,pfDown); */
/*     VMC_vec_branch[i]->opt_equilSweeps=100; //0; */
/*     VMC_vec_branch[i]->opt_SampleSweeps=1000; //0; */
/*   } */

/*   //  VMC_combine.GetParams("params.dat"); */
/*   //  VMC_combine.EvaluateAll(); */
  

/* #pragma omp parallel for */
/*   for (int i=0;i<NumWalkers;i++){ */
/*     cerr<<"On walker "<<i<<endl; */
/*     VMC_vec[i]->EvaluateAll(); */
/*     VMC_vec[i]->VMC(true); */
/*   } */


/* #pragma omp parallel for */
/*   for (int i=0;i<NumWalkers;i++){ */
/*     cerr<<"On walker "<<i<<endl; */
/*     VMC_vec_branch[i]->EvaluateAll(); */
/*     VMC_vec_branch[i]->VMC(true); */
/*   } */
/*   cerr<<"STARTING DMC"<<endl; */
  
/*   for (int DMCStep=0;DMCStep<1000000;DMCStep++){ */
/*     cerr<<"On DMC Step "<<DMCStep<<endl; */
/*     QMCTools::CalcEnergy(VMC_vec); */
/*     for (int i=0;i<NumWalkers;i++){ */
/*       //      cerr<<"On walker number "<<i<<endl; */
/*       double T_branch=0.1; */
/*       double T_left=T_branch; */

/*       ///I think I'm not clearing vals[i] when I shoudl be!! */
/*       while (T_left >0){     */
/* 	vals[i].clear(); */
/* 	//	cerr<<"T_left is "<<T_left<<endl; */
/* 	for (list<HamiltonianClass*>::iterator ham_iter=VMC_vec[i]->Ham.begin(); */
/* 	     ham_iter!=VMC_vec[i]->Ham.end();ham_iter++){ */
/* 	  (*ham_iter)->AllConnected(vals[i],VMC_vec[i]->System,VMC_vec[i]->wf_list); */
/* 	} */
/* 	double Z_xxp=0.0; */
/* 	for (list<pair<SpinSwap,double> >::iterator iter= ++(vals[i].begin());iter!=vals[i].end();iter++){ */
/* 	  Z_xxp+=(*iter).second/(VMC_vec[i]->System.tau); */
/* 	} */
/* 	//	cerr<<"Z_xxp is "<<Z_xxp<<endl; */
/* 	//	double El=((*(vals[i].begin())).second-1)/(VMC_vec[i]->System.tau)-Z_xxp; */
/* 	double El=((*(vals[i].begin())).second-1)/(-VMC_vec[i]->System.tau)-Z_xxp; */
	
/* 	double El_check=0; */
	
/* 	//	for (int k=0;i<vals[k].size();k++) */
/* 	//	  El_check+=vals[k]/VMC_vec[k]->System.tau; */
     
/* 	for (list<pair<SpinSwap,double> >::iterator iter= vals[i].begin(); */
/* 	     iter!=vals[i].end();iter++){ */
/* 	  if (iter==vals[i].begin()){ */
/* 	    El_check+=(iter->second-1.0)/(VMC_vec[i]->System.tau); */
/* 	  } */
/* 	  else{ */
/* 	    El_check+=iter->second/(VMC_vec[i]->System.tau); */
/* 	  } */
/* 	} */

/* 	//	cerr<<"El is "<<El<<endl; */
/* 	//	cerr<<"El_check is "<<El_check<<endl; */
/* 	//	for (list<HamiltonianClass*>::iterator ham_iter=VMC_vec[i]->Ham.begin(); */
/* 	//	     ham_iter!=VMC_vec[i]->Ham.end();ham_iter++){ */
/* 	  //	  cerr<<"ENERGY: "<<(*ham_iter)->Energy(VMC_vec[i]->System,VMC_vec[i]->wf_list)<<endl; */
/* 	//	} */
	
	
/* 	double T_try=min(-log(Random.ranf())/Z_xxp,T_left); */
/* 	FullWeight[i]*=exp(-T_try*El); */
/* 	T_left=T_left-T_try; */
	
/* 	if (T_left>0){ */
/* 	  vals[i].pop_front(); */
/* 	  int mySample=QMCTools::Sample(vals[i],Random,weight[i]); */
/* 	  VMC_vec[i]->Apply(vals[i],VMC_vec[i]->System,VMC_vec[i]->wf_list,mySample); */
/* 	} */
/*       } */
/*     } */

/* //     //HACK! */
/* //     for (int i=0;i<FullWeight.size();i++) */
/* //       FullWeight[i]=1.0; */
/* //     //HACK! */

/*     QMCTools::Branch(VMC_vec, */
/* 		     VMC_vec_branch, */
/* 		     FullWeight, */
/* 		     Random); */

/* //     for (int i=0;i<NumWalkers;i++){ */
/* //       VMC_vec[i]->EvaluateAll(); */
/* //       VMC_vec_branch[i]->EvaluateAll(); */
/* //     } */

/*     for (int i=0;i<FullWeight.size();i++) */
/*       FullWeight[i]=1.0; */
/*   } */
}



  

};


#endif
