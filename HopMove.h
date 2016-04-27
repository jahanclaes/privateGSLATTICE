#ifndef HOP_MOVE_H
#define HOP_MOVE_H

class HopMove
{
 public:
  int start;
  int end;
  int spin;
  HopMove()
  {
    start=-1;
    end=-1;
    spin=-1;

  }

  HopMove(int &t_start, int &t_end, int &t_spin){
    start=t_start;
    end=t_end;
    spin=t_spin;
  }
  void ReverseHop()
  {
    swap(start,end);

  }
};


class ExchangeMove
{
 public:
  int site1;
  int spin1;
  int site2;
  int spin2;
  ExchangeMove(int &t_site1,int &t_spin1, int &t_site2, int &t_spin2){
    site1=t_site1; site2=t_site2; spin1=t_spin1; spin2=t_spin2;
  }
  ExchangeMove(){
    site1=-1; site2=-1; spin1=-1; spin2=-1;
  }
    
};

#endif
