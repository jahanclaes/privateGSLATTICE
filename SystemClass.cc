#include "SystemClass.h"
#include <fstream>
#include "input.h"

double SystemClass::minDist(dVec r1, dVec r2)
{
  assert(1==2);
  double minDist2=999999;
  for (int i=-1;i<=1;i++)
    for (int j=-1;j<=1;j++){
      dVec diff(r1[0]+a1[0]*N1*i+a2[0]*j*N2-r2[0],r1[1]+a1[1]*i*N1+a2[1]*j*N2-r2[1]);
      double dist2=dot(diff,diff);
      minDist2=min(dist2,minDist2);
    }
  return sqrt(minDist2);
}


dVec SystemClass::minDiff(dVec r1, dVec r2)
{
  assert(1==2);
  dVec minD;
  double minDist2=999999;
  for (int i=-1;i<=1;i++)
    for (int j=-1;j<=1;j++){
      dVec diff(r1[0]+a1[0]*N1*i+a2[0]*j*N2-r2[0],r1[1]+a1[1]*i*N1+a2[1]*j*N2-r2[1]);
      double dist2=dot(diff,diff);
      if (dist2<minDist2){
	minDist2=dist2;
	minD=diff;
      }					       

	
    }
  return minD;
}

dVec SystemClass::minDiffStable(dVec r1, dVec r2)
{
  assert(1==2);
  dVec minD;
  double minDist2=999999;
  for (int i=-1;i<=1;i++)
    for (int j=-1;j<=1;j++){
      dVec diff(r1[0]+a1[0]*N1*i+a2[0]*j*N2-r2[0],r1[1]+a1[1]*i*N1+a2[1]*j*N2-r2[1]);
      double dist2=dot(diff,diff);
      if (dist2<minDist2){
	minDist2=dist2;
	minD=diff;
      }					       
      else if (dist2==minDist2 && (diff[0]<minD[0] || (diff[0]==minD[0] && diff[1]<minD[1]))){
	minDist2=dist2;
	minD=diff;
	
      }
    }
  return minD;
}

void SystemClass::Move(int site, int end_site, int spin)
{
  int not_spin = (spin==1) ? -1 : 1;
  if (x(site)==2)
    x(site)=not_spin;
  else 
    x(site)=0;
  if (x(end_site)==not_spin)
    x(end_site)=2;
  else 
    x(end_site)=spin;

  return;
}

void SystemClass::GenerateKList()
{
   ifstream infile;
  infile.open("kPoints.txt");
  dVec k;
  kList.clear();
  while (!infile.eof()){
    //    cerr<<"Trying to read"<<endl;
    infile>>k[0];
    infile>>k[1];
#if NDIM==3
    infile>>k[2];
#endif
    //    cerr<<k[0]<<" "<<k[1]<<endl;
    if (!infile.eof())
      kList.push_back(k);
    //    cerr<<"Pushed back"<<endl;
  }
  infile.close();
  cerr<<"SIZE OF K LIST IS "<<kList.size()<<endl;

}


void SystemClass::GenerateRList()
{
  cerr<<"reading r points"<<endl;
  ifstream infile;
  infile.open("rPoints.txt");
  
  infile>>a1[0];
  infile>>a1[1];
#if NDIM==3
  infile>>a1[2];
#endif
  infile>>a2[0];
  infile>>a2[1];
#if NDIM==3
  infile>>a2[2];
#endif
  infile>>N1;
  infile>>N2;
  cerr<<a1<<" "<<a2<<" "<<N1<<" "<<N2<<endl;
  
  
  dVec r;
  rList.clear();
  while (!infile.eof()){
    infile>>r[0];
    infile>>r[1];
#if NDIM==3
    infile>>r[2];
#endif
    if (!infile.eof())
      rList.push_back(r);
  //cerr << r << endl;  
  }
  infile.close();
  cerr<<"done reading r points"<<endl;
  cerr<<"list size "<<rList.size()<<endl;
  
}



void SystemClass::ReadNeighbors()
{
  
  // ifstream infile;
  // infile.open("neighbors.txt");
  // int numNeighbors;
  // infile>>numNeighbors;
  // //    cerr<<"numNeighbors is "<<numNeighbors<<" "<<rList.size()<<endl;
  // neighbors.resize(rList.size(),numNeighbors);
  // int latticeSite;
  // int neighbor;
  // while (!infile.eof()){
  //   infile>>latticeSite;
  //   if (!infile.eof())
  //     for (int i=0;i<numNeighbors;i++)
  // 	infile>>neighbors(latticeSite,i);
  // }
}


void SystemClass::SetupABSites()
{
  assert(1==2);
  // ABSites.resize(rList.size());
  // ABSites(Range::all())=-1;
   
  // //Figure out what the AB Sites are 
  // //for the bipartite lattice!
  // //Should be eventually moved into system!
  // ABSites(0)=0;
  // bool notDone=true;
  // while (notDone){
  //   notDone=false;
  //   for (int i=0;i<neighbors.extent(0);i++)
  //     for (int j=0;j<neighbors.extent(1);j++){
  // 	if (ABSites(i)==0)
  // 	  ABSites(neighbors(i,j))=1;
  // 	else if (ABSites(i)==1)
  // 	  ABSites(neighbors(i,j))=0;
  // 	else if (ABSites(i)==-1)
  // 	  notDone=true;
  //     }
  // }
  // for (int i=0;i<ABSites.size();i++)
  //   cerr<<"ABSITES: "<<i<<" "<<ABSites(i)<<endl;
}


void SystemClass::Init(InputClass &input)
{
  GenerateRList();
  x.resize(rList.size());
  cerr<<"Reading neighbors"<<endl;
  ReadNeighbors();
  ////  cerr<<"Reading k list"<<endl;
  ////  GenerateKList();
  double doping=input.toDouble(input.GetVariable("doping"));
  
  //  SetupABSites();
  tau=0.1;
  cerr<<"Staggering"<<endl;
  Stagger(doping);
  cerr<<"I am done staggering"<<endl;
  int numUp=CountElectrons(-1,x.size(),1);
  cerr<<"my numUp is "<<numUp<<endl;
  for (int i=0;i<x.size();i++)
    cerr<<x(i)<<" ";
  cerr<<endl;
  for (int i=0;i<numUp;i++){
    dVec garbage;
    kList.push_back(garbage);
  }
  cerr<<"SIZE OF K LIST IS "<<kList.size()<<endl;
  cerr<<"SystemClass init done"<<endl;
}

//Assumes your neighbor is on a different
//bipartite lattice
void SystemClass::Stagger_kondo(double &doping)
{
  for (int i=0;i<x.size();i++){
    x(i)=0;
  }
  int numUp=0.5*doping*(x.size()/2);
  cerr<<"I want to set the num Up to be "<<numUp<<endl;
  {
    int pos=-2;
    for (int ii=0;ii<numUp;ii++){
      pos+=2;
      if (pos>=x.size()/2)
	pos=1;
      x(pos)=1;
    }
  }
  {
    int pos=-1;
    for (int ii=0;ii<numUp;ii++){
      pos+=2;
      if (pos>=x.size()/2)
	pos=0;
      x(pos) = ( (x(pos)==1) ? 2 : -1);
    }
  }
  x(x.size()/2)=1;
  for (int ii=x.size()/2+1;ii<x.size();ii++){
    x(ii)=x(ii-1)*-1;
  }

  cerr<<"The current system is "<<endl;
  for (int i=0;i<x.size();i++){
    cerr<<x(i)<<endl;
  }
 return;
    
  //  assert(1==2);
  //  RandomClass Random;
  x(0)=1;
  //  for (int i=1;i<x.size();i++){
  //    x(i)=x(i-1)*-1;
  //  }
  x(0)=2;
  x(1)=2;
  x(2)=1;
  x(3)=-1;

  //  x(x.size()-6)=0;
  //  x(x.size()-5)=0;
  //  x(x.size()-4)=0;
  //  x(x.size()-3)=0;
  //  x(x.size()-2)=0;
  //  x(x.size()-1)=0;
  //  for (int i=1;i<x.size();i++){
  //    x(i)=x(i-1)*-1+1;
  //  }
  //  swap(x(0),x(1));
  //  for (int i=0;i<x.size();i++)
  //    swap(x(i),x(Random.randInt(x.size())));
  //  swap(x(2),x(3));
  //     int N=x.size();
  //   x(Range::all())=1;
  //    for (int i=0;i<N/2;i++)
  //     x(i)=0;
  
    
}



//Assumes your neighbor is on a different
//bipartite lattice
void SystemClass::Stagger(double &doping)
{
  for (int i=0;i<x.size();i++){
    x(i)=0;
  }
  int numUp=0.5*doping*(x.size());
  cerr<<"I want to set the num Up to be "<<numUp<<endl;
  {
    int pos=-2;
    for (int ii=0;ii<numUp;ii++){
      pos+=2;
      if (pos>=x.size())
	pos=1;
      x(pos)=1;
    }
  }
  {
    int pos=-1;
    for (int ii=0;ii<numUp;ii++){
      pos+=2;
      if (pos>=x.size())
	pos=0;
      x(pos) = ( (x(pos)==1) ? 2 : -1);
    }
  }

  cerr<<"The current system is "<<endl;
  for (int i=0;i<x.size();i++){
    cerr<<x(i)<<endl;
  }
 return;
    
}


// bool SystemClass::inList(int site,TinyVector<int,6> &myList)
// {
//   assert(1==2);
//   for (int i=0;i<6;i++){
//     if (site==myList(i))
//       return true;
//   }
//   return false;
// }



// int SystemClass::calcABSign()
// {
//   assert(1==2);
//   Array<int,1> temp_x(x.size());
//   temp_x=x;
//   bool ok=false;
//   int needSwap=0;
//   cerr<<"INIT: ";
//   for (int i=0;i<temp_x.size();i++)
//     cerr<<temp_x(i)<<" ";
//   cerr<<endl;
//   while (!ok){
//     ok=true;
//     for (int i=0;i<temp_x.size()-1;i++){
//       if (temp_x(i)==1 && i<temp_x.size()/2){
// 	ok=false;
//       }
//       if (temp_x(i)==1 && temp_x(i+1)==0){
// 	swap(temp_x(i),temp_x(i+1));
// 	needSwap++;
//       }
//     }
//   }
//   return (needSwap %2 ? -1: 1);
// }

// bool SystemClass::notNeighbor(int i,int neighbor)
// {
//   assert(1==2);
//   for (int j=0;j<neighbors.extent(1);j++){
//     if (neighbors(i,j)==neighbor)
//       return false;
//   }
//   return true;
// }

void SystemClass::Swap(int i,int j)
{
  swap(x(i),x(j));

}

//Han-Yi Chou 11/12
int SystemClass::CountElectrons(int i,int j,int spin)
{
  int numElectrons=0;
  int myMin=min(i,j);
  int myMax=max(i,j);
  for (int count=myMin+1;count<myMax;count++){
    numElectrons += ( (x(count)==spin  || x(count)==2) ? 1: 0 );
    //    numElectrons+=abs(x(count));
  }
  return numElectrons;

}
