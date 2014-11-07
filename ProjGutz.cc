#include "ProjGutz.h"
#include "SystemClass.h"
#include "MatrixOps.h"


using namespace std;
complex<double> 
ProjectedGutzweilerClass::evaluateRatio_honeycomb(SystemClass &system,TinyVector<int,6>  &honeycomb_locs,
						  TinyVector<int,6> &honeycomb_backup,int amt)
{
    vector<int> upList;
    vector<int> downList;
    for (int i=0;i<6;i++){
      if (system.x(honeycomb_locs[i])==0)
	upList.push_back(honeycomb_locs[i]);
      else
	downList.push_back(honeycomb_locs[i]);
    }
    
    Array<double,2> upDet(upList.size(),upList.size());
    upDet=0.0;
    for (int i=0;i<upDet.extent(0);i++)
      for (int j=0;j<upDet.extent(1);j++){
	for (int ki=0;ki<system.kList.size();ki++){
	  int pos2=DetPos(upList[j]);
	  upDet(i,j)+=eigs(ki,upList[i])*MInverse[0](ki,pos2);
	}
      }

    double rat1;
    if (upList.size()!=0)
      rat1=MatrixOps::Determinant(upDet);
    else
      rat1=1.0;

    Array<double,2> downDet(downList.size(),downList.size());
    downDet=0.0;
    for (int i=0;i<downDet.extent(0);i++)
      for (int j=0;j<downDet.extent(1);j++){
	for (int ki=0;ki<system.kList.size();ki++){
	  int pos2=DetPos(downList[j]);
	  downDet(i,j)+=eigs(ki,downList[i])*MInverse[1](ki,pos2);
	}
      }

    double rat2;
    if (downList.size()!=0)
      rat2=MatrixOps::Determinant(downDet);
    else 
      rat2=1.0;

    return rat1*rat2*(-1);
}



complex<double> 
ProjectedGutzweilerClass::evaluateRatio(SystemClass &system,TinyVector<int,4>  &rhombus_locs,
					TinyVector<int,4> &rhombus_backup)
{
  double oldSign=RealSign;
    vector<int> upList;
    vector<int> downList;
    for (int i=0;i<4;i++){
      if (system.x(rhombus_locs[i])==0)
	upList.push_back(rhombus_locs[i]);
      else
	downList.push_back(rhombus_locs[i]);
    }
    
    Array<double,2> upDet(upList.size(),upList.size());
    upDet=0.0;
    for (int i=0;i<upDet.extent(0);i++)
      for (int j=0;j<upDet.extent(1);j++){
	for (int ki=0;ki<system.kList.size();ki++){
	  int pos2=DetPos(upList[j]);
	  upDet(i,j)+=eigs(ki,upList[i])*MInverse[0](ki,pos2);
	}
      }

    double rat1;
    if (upList.size()!=0)
      rat1=MatrixOps::Determinant(upDet);
    else
      rat1=1.0;

    

    Array<double,2> downDet(downList.size(),downList.size());
    downDet=0.0;
    for (int i=0;i<downDet.extent(0);i++)
      for (int j=0;j<downDet.extent(1);j++){
	for (int ki=0;ki<system.kList.size();ki++){
	  int pos2=DetPos(downList[j]);
	  downDet(i,j)+=eigs(ki,downList[i])*MInverse[1](ki,pos2);
	}
      }

    double rat2;
    if (downList.size()!=0)
      rat2=MatrixOps::Determinant(downDet);
    else 
      rat2=1.0;


    Swap(rhombus_locs[0],rhombus_locs[1]);
    Swap(rhombus_locs[1],rhombus_locs[2]);
    Swap(rhombus_locs[2],rhombus_locs[3]);
    double newSign=RealSign;
    Swap(rhombus_locs[2],rhombus_locs[3]);
    Swap(rhombus_locs[1],rhombus_locs[2]);
    Swap(rhombus_locs[0],rhombus_locs[1]);
    

    
    return rat1*rat2*oldSign*newSign;
}



void 
ProjectedGutzweilerClass::RotateHoneycomb(TinyVector<int,6> &honeycomb,
					  
					  int amt)
{

   for (int i=0;i<honeycomb.length();i++)
     honeycomb_backup(i)=DetPos(honeycomb(i));
   for (int i=0;i<honeycomb.length();i++)
     DetPos(honeycomb((i+amt+6) %6))=honeycomb_backup(i);
}


void 
ProjectedGutzweilerClass::Init(SystemClass &system)
{
  
  NumSpinUp=system.x.size()/2;
  for (int i=0;i<2;i++){
    Dets[i].resize(NumSpinUp,NumSpinUp);
    MInverse[i].resize(NumSpinUp,NumSpinUp);
  }
  u[0].resize(NumSpinUp);
  u[1].resize(NumSpinUp);
  MInverseu.resize(NumSpinUp);
  MInverse_k.resize(NumSpinUp);
  DetPos.resize(system.x.size());
  //  eigs.resize(NumSpinUp,NumSpinUp);
  GetEigs(system);
  FillDet(system,0);
  FillDet(system,1);
  ParticleOrder.resize(system.x.size());
  for (int i=0;i<ParticleOrder.size();i++){
    ParticleOrder(i)=i;
  }
  RealSign=1;
    
}

void
ProjectedGutzweilerClass::Swap(int i,int j)
{
  swap(DetPos(i),DetPos(j));
  swap(ParticleOrder(i),ParticleOrder(j));
  if (abs(i-j) % 2==1)
    RealSign*=-1;
}

void 
ProjectedGutzweilerClass::FillDet(SystemClass &system,int spin)
{
  int numUp=0;
  int numDown=0;
  for (int i=0;i<system.x.size();i++){
    if (system.x(i)==0)
      numUp++;
    else 
      numDown++;
  }
  //  cerr<<"Num up and down is "<<numUp<<" "<<numDown<<endl;
  
  Array<double,2> &Det(Dets[spin]);
  Det=0.0;
  int detLoc=0;
  for (int ri=0;ri<system.x.size();ri++)
    if (system.x(ri)==spin){
      for (int ki=0;ki<system.kList.size();ki++){
	assert(detLoc<Det.extent(0));
	assert(ki<Det.extent(1));
	Det(detLoc,ki)=eigs(ki,ri);
	DetPos(ri)=detLoc;
      }
      detLoc++;
    }
  MInverse[spin]=MatrixOps::Inverse(Det);  
}


//BUG: I THINK THE SIGN IS OFF HERE BECAUSE IT
//DEPENDS ON HOW THE UP AND DOWN SPINS ARE INTERLACED
complex<double> 
ProjectedGutzweilerClass::evaluate(SystemClass &system)
{
  FillDet(system,0);
  double spinUp=MatrixOps::Determinant(Dets[0]);
  FillDet(system,1);
  double spinDown=MatrixOps::Determinant(Dets[1]);
  cerr<<"Spin up and spin down is "<<spinUp<<" "<<spinDown<<endl;
  return spinUp*spinDown;
  //  return log(abs(spinUp))+log(abs(spinDown));
}

//BUG: I THINK THE SIGN IS OFF HERE BECAUSE IT
//DEPENDS ON HOW THE UP AND DOWN SPINS ARE INTERLACED
complex<double> 
ProjectedGutzweilerClass::logevaluate(SystemClass &system,int &sign)
{
  FillDet(system,0);
  double spinUp=MatrixOps::Determinant(Dets[0]);
  FillDet(system,1);
  double spinDown=MatrixOps::Determinant(Dets[1]);
  sign = int(abs(spinUp)/spinUp * abs(spinDown)/spinDown);
  //  CheckMInverse();
  return log(abs(spinUp))+log(abs(spinDown));  

}



//BUG: I THINK THE SIGN IS OFF HERE BECAUSE IT
//DEPENDS ON HOW THE UP AND DOWN SPINS ARE INTERLACED
complex<double> 
ProjectedGutzweilerClass::evaluateRatio(SystemClass &system,int pos1, int pos2)
{
  double ratio1=0.0;
  for  (int ki=0;ki<system.kList.size();ki++){
    ratio1+=eigs(ki,pos1)*MInverse[system.x(pos1)](ki,DetPos(pos1));
  }
  double ratio2=0.0;
  for  (int ki=0;ki<system.kList.size();ki++){
    ratio2+=eigs(ki,pos2)*MInverse[system.x(pos2)](ki,DetPos(pos2));
  }
//   ///CHECK!
//   cerr<<"CHECKING: "<<pos1<<" "<<pos2<<" "<<evaluateRatio_check(system,pos1,pos2)<<" "<<ratio2*ratio1*-1<<endl;
//   system.Swap(pos1,pos2);
//   Swap(pos1,pos2);
//   evaluateRatio_check(system,pos1,pos2);
//   system.Swap(pos1,pos2);
//   Swap(pos1,pos2);
// //   //CHECK!
//   ///The -1 because it only matters during the swap
//   cerr<<"Ratios: "<<ratio1<<" "<<ratio2<<endl;
//  cerr<<"Real Sign is "<<pos1<<" "<<pos2<<" "<<RealSign<<endl;
//  return RealSign*ratio2*ratio1;
  return (-1)*ratio2*ratio1;
//  return (-1)*ratio2*ratio1;
}

 
complex<double> 
ProjectedGutzweilerClass::evaluateRatio_check(SystemClass &system, int swap1, int swap2)
{
   system.Swap(swap1,swap2);
   Swap(swap1,swap2);
   int signa;
   complex<double> before=logevaluate(system,signa);
   cerr<<"Before: "<<exp(before.real())<<" "<<before.real()<<endl;
   system.Swap(swap1,swap2);
   Swap(swap1,swap2);
   int signb;
   complex<double> after=logevaluate(system,signb);
   cerr<<"After: "<<exp(after.real())<<" "<<after.real()<<endl;
   return exp(after.real()-before.real())*signa*signb;
}





void 
ProjectedGutzweilerClass::CheckMInverse()
{
  Array<double,2> MInverse_check(MInverse[0].extent(0),MInverse[0].extent(1));
  bool correct=true;
  MInverse_check=MatrixOps::Inverse(Dets[0]);
  for (int i=0;i<MInverse_check.extent(0);i++)
    for (int j=0;j<MInverse_check.extent(1);j++)
      if (fabs(MInverse_check(i,j)-MInverse[0](i,j))>1e-3){
	correct=false;
	cerr<<"CORRECT  INCORRECT"<<endl;
	cerr<<i<<" "<<j<<" "<<MInverse_check(i,j)<<" "<<MInverse[0](i,j)<<endl;
      }
  if (!correct)
    cerr<<"The first determinant is not correct"<<endl;
  else 
    cerr<<"First determinant: correct"<<endl;
  correct=true;
  
  MInverse_check=MatrixOps::Inverse(Dets[1]);
  for (int i=0;i<MInverse_check.extent(0);i++)
    for (int j=0;j<MInverse_check.extent(0);j++)
      if (fabs(MInverse_check(i,j)-MInverse[1](i,j))>1e-3){
	correct=false;
	cerr<<"CORRECT  INCORRECT"<<endl;
	cerr<<i<<" "<<j<<" "<<MInverse_check(i,j)<<" "<<MInverse[1](i,j)<<" "<<MInverse_check(i,j)-MInverse[1](i,j)<<endl;
      }
  if (!correct)
    cerr<<"The second determinant is not correct"<<endl;
  else 
    cerr<<"Second determinant: correct"<<endl;
}

void 
ProjectedGutzweilerClass::UpdateDets(SystemClass &system,int swap1, int swap2)
{
  Dets[system.x(swap1)](DetPos(swap1),Range::all())=eigs(Range::all(),swap1);
  Dets[system.x(swap2)](DetPos(swap2),Range::all())=eigs(Range::all(),swap2);
  firstIndex i;  secondIndex j;  thirdIndex k;
  int myDet=system.x(swap1);
  u[myDet](Range::all())=eigs(Range::all(),swap1)-eigs(Range::all(),swap2);
  MInverseu=sum(MInverse[myDet](j,i)*u[myDet](j),j);
  MInverse_k(Range::all()) = MInverse[myDet](Range::all(),DetPos(swap1));
  double s=sum(u[myDet](i)*MInverse_k(i))+1.0;
  //    cerr<<"s1 is "<<s<<endl;
  MInverse[myDet]=MInverse[myDet]-1.0/s*(MInverseu(j)*MInverse_k(i));

  myDet=system.x(swap2);
  u[myDet](Range::all())=eigs(Range::all(),swap2)-eigs(Range::all(),swap1);
  MInverseu=sum(MInverse[myDet](j,i)*u[myDet](j),j);
  MInverse_k(Range::all()) = MInverse[myDet](Range::all(),DetPos(swap2));
  s=sum(u[myDet](i)*MInverse_k(i))+1.0;
  //    cerr<<"s2 is "<<s<<endl;
  MInverse[myDet]=MInverse[myDet]-1.0/s*(MInverseu(j)*MInverse_k(i));
  //  CheckMInverse();

}



void 
ProjectedGutzweilerClass::GetEigs(SystemClass &system)
{
  eigs.resize(system.kList.size(),system.rList.size());
  ifstream infile;
  infile.open("eigs.txt");
  for (int i=0;i<system.kList.size();i++)
    for (int j=0;j<system.rList.size();j++)
      infile>>eigs(i,j);
  infile.close();
}


