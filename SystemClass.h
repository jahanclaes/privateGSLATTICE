#ifndef SYSTEM_CLASS_H
#define SYSTEM_CLASS_H
#include <vector>
#include "Blitz.h"
class InputClass;

class SystemClass
{
 public:
  Array<int,1> x;
  //  Array<int,1> ABSites;
  void Init(InputClass &input);
  double tau;
  //  Array<int,2> neighbors;
  vector<dVec> rList;
  vector<dVec> kList;
  int N1; int N2;
  dVec a1; dVec a2;
  double minDist(dVec r1, dVec r2);
  dVec minDiff(dVec r1, dVec r2);
  dVec minDiffStable(dVec r1, dVec r2);
  void GenerateRList();
  void GenerateKList();
  void ReadNeighbors();
  void SetupABSites();
  void Stagger(double &doping);
  void Swap(int i, int j);
  void Move(int site, int end_site, int spin);
  int CountElectrons(int i,int j,int spin);
  //  bool notNeighbor(int i,int neighbor);
  //  bool inList(int site,TinyVector<int,6> &myList);
  //  int calcABSign();
};

#endif
