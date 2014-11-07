#ifndef DERIV_CLASS_H
#define DERIV_CLASS_H

class DerivClass
{						       
public:


  vector<complex<double> > Psi_alpha_over_Psip;
  vector<complex<double> > El_times_Psi_alpha_over_Psip;
  double NumTimes;
  complex<double> E_avgp;
  complex<double> E_avg2p;

  Eigen::MatrixXd S;
  Eigen::MatrixXd S_inverse;

  Eigen::MatrixXd Sp;

  //  Array<double,2> S;
  //  Array<double,2> S_inverse;

  bool doS;
  void Init(int NumParams)
  {
    assert(1==2);
  }
  void Init(int NumParams,bool t_doS)
  {
    doS=t_doS;
    Psi_alpha_over_Psip.resize(NumParams);
    El_times_Psi_alpha_over_Psip.resize(NumParams);
    if (doS){
      S.resize(NumParams,NumParams);

      S_inverse.resize(NumParams,NumParams);
      Sp.resize(NumParams,NumParams);
    }

   
  }
  void Clear()
  {
    E_avgp=0.0;
    E_avg2p=0.0;
    NumTimes=0;
    for (int i=0;i<Psi_alpha_over_Psip.size();i++){
      El_times_Psi_alpha_over_Psip[i]=0.0;
      Psi_alpha_over_Psip[i]=0.0;
      if (doS){
	for (int j=0;j<Psi_alpha_over_Psip.size();j++){
	  S(i,j)=0;
	}
      }
    }
    
  }

  void Combine(vector<DerivClass*> &Derivs)
  {
    cerr<<"THE COMBINE SIZE IS "<<Derivs.size()<<endl;
    for (int i=0;i<Derivs.size();i++){
      NumTimes+=Derivs[i]->NumTimes;
      E_avgp+=Derivs[i]->E_avgp;
      E_avg2p+=Derivs[i]->E_avg2p;
      if (doS){
	for (int ii=0;ii<Psi_alpha_over_Psip.size();ii++)
	  for (int j=0;j<Psi_alpha_over_Psip.size();j++){
	    S(ii,j)+=Derivs[i]->S(ii,j);
	  }
      }
      for (int j=0;j<Psi_alpha_over_Psip.size();j++){
	Psi_alpha_over_Psip[j]+=Derivs[i]->Psi_alpha_over_Psip[j];
	El_times_Psi_alpha_over_Psip[j]+=Derivs[i]->El_times_Psi_alpha_over_Psip[j];
      }
    }
  }


  void ParallelCombine(CommunicatorClass &myComm)
  {
    NumTimes=myComm.Sum(NumTimes); 
    E_avgp=myComm.Sum(E_avgp);
    E_avg2p=myComm.Sum(E_avg2p);
    vector<complex<double>  > Psi_alpha_over_Psi2(Psi_alpha_over_Psip.size(),0);
    vector<complex<double>  > El_times_Psi_alpha_over_Psi2(El_times_Psi_alpha_over_Psip.size(),0);
    myComm.Sum(Psi_alpha_over_Psip,Psi_alpha_over_Psi2);
    myComm.Sum(El_times_Psi_alpha_over_Psip,El_times_Psi_alpha_over_Psi2);
    for (int i=0;i<Psi_alpha_over_Psip.size();i++){
    Psi_alpha_over_Psip[i]=Psi_alpha_over_Psi2[i];
    }
    for (int i=0;i<El_times_Psi_alpha_over_Psip.size();i++){
      El_times_Psi_alpha_over_Psip[i]=El_times_Psi_alpha_over_Psi2[i];
    }
    if (doS){
      Sp=S;
      myComm.Sum(Sp,S);
    }
    
    
/*       NumTimes+=Derivs[i]->NumTimes; */
/*       E_avgp+=Derivs[i]->E_avgp; */
/*       E_avg2p+=Derivs[i]->E_avg2p; */
/*       if (doS){ */
/* 	for (int ii=0;ii<Psi_alpha_over_Psip.size();ii++) */
/* 	  for (int j=0;j<Psi_alpha_over_Psip.size();j++){ */
/* 	    S(ii,j)+=Derivs[i]->S(ii,j); */
/* 	  } */
/*       } */
/*       for (int j=0;j<Psi_alpha_over_Psip.size();j++){ */
/* 	Psi_alpha_over_Psip[j]+=Derivs[i]->Psi_alpha_over_Psip[j]; */
/* 	El_times_Psi_alpha_over_Psip[j]+=Derivs[i]->El_times_Psi_alpha_over_Psip[j]; */
/*       } */
//    }

      }


  
  void GetSInverse()
  {
    for (int i=0;i<S.cols();i++){
      for (int j=0;j<S.rows();j++){
        S(i,j)=S(i,j)/NumTimes-(Psi_alpha_over_Psip[i].real()/NumTimes)*
          (Psi_alpha_over_Psip[j].real()/NumTimes);
      }
    }
    for (int i=0;i<S.cols();i++)
      S(i,i)=S(i,i)+1e-8;
	//    S_inverse=MatrixOps::Inverse(S);
    S_inverse=S.inverse();
    //    cerr<<"Diagonal "<<endl;
    //    cerr<<S.diagonal()<<endl;
    //    cerr<<endl;
    //    cerr<<S_inverse<<endl;
    //    cerr<<endl;
    //    cerr<<endl;
    //    exit(1);
  }
  
  void Add(double &val_El, 
	   Array<complex<double>,1> &derivs)
  {
    NumTimes=NumTimes+1;
    E_avgp+=val_El;
    E_avg2p+=val_El*val_El;
    assert(derivs.size()==Psi_alpha_over_Psip.size());
    for (int i=0;i<derivs.size();i++){

      Psi_alpha_over_Psip[i]+=derivs(i);
      El_times_Psi_alpha_over_Psip[i]+=derivs(i)*val_El;
      if (doS){
	for (int j=0;j<derivs.size();j++){
	  S(i,j)+=derivs(i).real()*derivs(j).real();
	}
      }
    }
    
  }
  complex<double> ComputeDerivSR(int param)
    {
      //switched i and param     
      double myDeriv=0;
      for (int i=0;i<S.cols();i++){
        myDeriv=myDeriv+(S_inverse(i,param)*
                         2.0 *(El_times_Psi_alpha_over_Psip[i].real()/NumTimes -
			       (Psi_alpha_over_Psip[i].real()/NumTimes)*
                               (E_avgp.real()/NumTimes)));
      }
      return myDeriv;
    }
  complex<double> ComputeEnergy()
    {
      return (E_avgp.real()/NumTimes);
    }
  complex<double> ComputeVariance()
    {
      return (E_avgp.real()/NumTimes)*(E_avgp.real()/NumTimes)-E_avg2p.real()/NumTimes;
    }
  complex<double> ComputeDerivp(int param)
    {
      return 2.0 * (El_times_Psi_alpha_over_Psip[param].real()/NumTimes - (Psi_alpha_over_Psip[param].real()/NumTimes)*(E_avgp.real()/NumTimes));

    }
  complex<double> ComputeDerivp_imag(int param)
    {
      return 2.0 * (El_times_Psi_alpha_over_Psip[param].imag()/NumTimes - (Psi_alpha_over_Psip[param].imag()/NumTimes)*(E_avgp.real()/NumTimes));

    }


};

#endif
