#ifndef input_h
#define input_h
#include <vector>
#include <string>

using namespace std;

class InputTree
{
 public:
  string name;
  vector<InputTree*> sections;
  map<string,string> input;
  InputTree* parent;

};


class InputClass
{
 public:
  int toInteger(string myString)
  {
    return atoi(myString.c_str());
  }
  double toDouble(string myString)
  {
    return atof(myString.c_str());
  }
  bool IsVariable(string myVar)
    {
      return tree.input.count(myVar)==1;
    }
  string GetVariable(string myVar)
  {
    assert(IsVariable(myVar));
    return tree.input[myVar];
  }
  
  InputTree tree;
  void RemoveTrailingWhiteSpace(string &s)
  {
    int lastPos=s.find_last_not_of(" \t\n\v\f\r");
    if (lastPos!=string::npos){
      s=s.substr(0,lastPos+1);
    }
  }

  void RemoveLeadingWhiteSpace(string &s)
  {
    int firstPos=s.find_first_not_of(" \t\n\v\f\r");
    if (firstPos!=string::npos)
      s=s.substr(firstPos,s.size());
  }
  void RemoveWhiteSpace(string &s){
    RemoveTrailingWhiteSpace(s);
    RemoveLeadingWhiteSpace(s);
  }
  void Read(ifstream &infile)
  {
    tree.name="root";
    InputTree *currentTree=&tree;

    while (!infile.eof()){
      string line;
      infile>>line;
      RemoveWhiteSpace(line);
      if (line[line.size()-1]==':'){
	line.substr(0,line.size()-1);
	InputTree *t=new InputTree();
	t->name=line;
	t->parent=currentTree;
	currentTree->sections.push_back(t);
	currentTree=t;
      }
      else if (line=="END"){
	currentTree=currentTree->parent;
      }
      else if (line.size()==0){
      }
      else {
	int eqLoc=line.find('=');
	assert(eqLoc!=string::npos);
	string token=line.substr(0,eqLoc);
	RemoveWhiteSpace(token);
	string outVals=line.substr(eqLoc+1);
	RemoveWhiteSpace(outVals);
	currentTree->input[token]=outVals;
      }
      

    }

  };
  
  


};


#endif
