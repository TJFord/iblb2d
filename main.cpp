#include <cstdlib>
#include <iostream>

#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>

using namespace std;

#include "lb.h"
#include "boundary.h"
#include "units.h"

int main(int argc, char *argv[])
{
  //fstream in("input.txt",ios::in);
  string in="input.txt";
  //string in="input_simpleSquare.txt";
  string out="rst.txt";

  LB a;
  a.readInput(in);
  a.printInfor();
  /*Units para(0.1, 1.0, 1.0, 1e-6,1.0);
  para.calculateLBPara();
  double Re2 = para.getRe();
  */
 // cout<<"Re2"<<Re2<<endl;
  int nSave = 100;
  a.init();
  
  clock_t begin = clock();
  for (int i=0;i<1200;i++)
  {
    //  cout<<"loop i="<<i<<endl;
    //a.applyBC();
    //a.stream();
    //a.collide();
    a.collideSwap();
    a.streamSwap();
    a.applyBC();
    if (i%nSave ==0 )
      a.output(out);
  }
  clock_t end = clock();
  double elapsedSecs = double(end-begin)/CLOCKS_PER_SEC;
  cout<<"time elapsed "<<elapsedSecs<<endl;
  return 0;
}
