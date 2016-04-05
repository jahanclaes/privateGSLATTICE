#ifndef EXACT_H
#define EXACT_H
#define NSITES 4

#include <bitset>
#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
using namespace std;

bool bitsetLess( std::bitset<NSITES> const& left, 
		 std::bitset<NSITES> const& right )
{
  return left.to_ulong() < right.to_ulong();
}


class Hamiltonian2
{
 public:
  
  vector<pair<int,int> > bondList;
  double J;
  int NumSites;
  int NumUp;
  vector<bitset<NSITES> > basis;
  void mult(vector<double> const & vin, vector<double> & vout)   const
  {
    for (int i=0;i<basis.size();i++){
      vout[i]=0.0;
      MultForBit(basis[i],vin,vout);

    }
    
  }

  void Init()
  {
    Init("J1.txt");
  }
  void Init(string fileName)
  {
    NumSites=NSITES;
    NumUp=NumSites/2;
    GenerateBasis();


    ifstream infile;
    infile.open(fileName.c_str());

    if (!infile){
      cerr<<"Could not open heisenberg hamiltonian "<<fileName<<endl;
      cerr<<"Aborting!"<<endl;
      exit(1);
    }
    
    while (!infile.eof()){
      int i;
      int j;
      infile>>i;
      infile>>j;
      if (!infile.eof())
	bondList.push_back(make_pair(i,j));
    }
    infile.close();
    J=1.0;
  }

  double FindIndex(bitset<NSITES> const &myBit)    const
  {
    vector<bitset<NSITES > >::const_iterator iter=std::lower_bound (basis.begin(), basis.end(),myBit,bitsetLess);
    int index=iter-basis.begin();
    return index;
  }

  void MultForBit(bitset<NSITES> myBit,vector<double> const & vin,vector<double> &vout)  const
   { 

     double Jz=1.0;
     double Jx=1.0; 
     int ii=FindIndex(myBit); 
     complex<double> sum=0.0;  
     double diagonal=0.0;  
     for (int counter=0;counter<bondList.size();counter++){  
       int i=bondList[counter].first;  
       int j=bondList[counter].second;  
       if (myBit[i]==myBit[j]){  
	 diagonal+=0.25*Jz;  
       }  
       else {  
     	 diagonal-=0.25*Jz;  
	 myBit[i].flip();  
	 myBit[j].flip();  
	 int jj=FindIndex(myBit);  

	 myBit[i].flip();  
	 myBit[j].flip();  
	 vout[ii]+=vin[jj]*0.5*Jx; 
       }  
       
       
     } 
     vout[ii]+=vin[ii]*diagonal; 
   }

  bool NextBit(bitset<NSITES> &myBit)
  {
    int i=0;
    while (i!=NSITES && myBit[i]==1){
      myBit.set(i,0);
      i++;
    }
    if (i==NSITES)
      return false;
    myBit.set(i,1);
    return true;
  }

  void GenerateBasis()
  {
    cerr<<"Going to generate basis "<<endl;
    //    exit(1);
    bool continueMe;
    bitset<NSITES> basisElement;
    basisElement.reset();
    for (int i=0;i<NumUp;i++)
      basisElement.set(i,1);
    basis.push_back(basisElement);
    cerr<<basisElement<<endl;
    while (continueMe){
      continueMe=NextBit(basisElement);
      while (continueMe && basisElement.count()!=NumUp)
	continueMe=NextBit(basisElement);
      if (continueMe)
	basis.push_back(basisElement);
	//      cerr<<basisElement<<endl;
    }
    std::sort( basis.begin(), basis.end(), bitsetLess );
    cerr<<"My basis size is "<<basis.size()<<endl;
  }
  

  
};


#endif
