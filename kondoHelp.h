#ifndef KONDO_HELP
#define KONDO_HELP
#include <vector>
using namespace std;

class KondoHelp
{
 public:
  vector<int> layer1Sites;
  vector<int> layer2Sites;
  
  vector<pair<int,int> > heisenbergBonds;
  void Read(string &layer1File, string &layer2File, string &bondFile)
  {
    ifstream infile;
    infile.open(layer1File.c_str());
    assert(infile);
    while (!infile.eof()){
      int site;
      infile>>site;
      if (!infile.eof())
	layer1Sites.push_back(site);
    }
    infile.close();

    infile.open(layer2File.c_str());
    assert(infile);
    while (!infile.eof()){
      int site;
      infile>>site;
      if (!infile.eof())
	layer2Sites.push_back(site);
    }
    infile.close();


    infile.open(bondFile.c_str());
    assert(infile);
    while (!infile.eof()){
      int site1;
      int site2;
      infile>>site1;
      infile>>site2;
      if (!infile.eof())
	heisenbergBonds.push_back(make_pair(site1,site2));
    }
    infile.close();
  }

};

#endif
