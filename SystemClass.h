#ifndef SYSTEM_CLASS_H
#define SYSTEM_CLASS_H
#include <vector>
#include "Blitz.h"


class SystemClass
{
 public:
  Array<int,1> x;
  Array<int,1> ABSites;
  void Init();
  double tau;
  vector<vector<int > > neighbors;
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
  void ReadNeighbors_old();
  void SetupABSites();
  void Stagger();
  void GetHoneycomb(int i,TinyVector<int,6> &locs);
  void GetRhombus(int i,TinyVector<int,4> &locs);
  void GetRhombusB(int i,TinyVector<int,4> &locs);
  bool notNeighbor(int i,int neighbor);
  void Swap(int i, int j);
  void Flip(int i);
  void RotateHoneycomb(TinyVector<int,6> &honeycomb,
		       TinyVector<int,6> &honeycomb_backup,
		       int amt);
  void RotateRhombus(TinyVector<int,4> &honeycomb,
		     TinyVector<int,4> &honeycomb_backup,
		     int amt);

  bool inList(int site,TinyVector<int,6> &myList);
  int calcABSign();
};

#endif
