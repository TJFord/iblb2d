#ifndef IBM_2D_H
#define IBM_2D_H

#include "lb.h"
#include "cell.h"

class IBM{
public:
  LB *plb;
  Cell *pCell;
  int nc;//# of cells 
  double w1[4],w2[4];
  
  IBM(LB*,Cell*,int nc_=1);
  ~IBM();

  void interpret();
  void spread();
  void ph1(double r);
  void ph2(double r);
  void ph3(double r);
  double phCos(double r);
  void output(const std::string filename);
  void printInfor();
  void periodic();
  void moveSolidTo(double x, double y);
};
#endif
