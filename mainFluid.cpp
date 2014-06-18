#include <cstdlib>
#include <iostream>

#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>
//#include <climits>

using namespace std;

#include "lb.h"
#include "boundary.h"
#include "units.h"
#include "cell.h"
#include "ibm.h"

int main(int argc, char *argv[])
{
  //fstream in("input.txt",ios::in);
  //string in="input.txt";
  //string out="rst.txt";
  string cin="cellInput.txt";
  //string in="inputForce.txt";
  //string fin="inputChannel.txt";
  //string fin="shear.txt";
  string fin="vonkarman.txt";
  string fgeom="fgeom.txt";
  string cellout="cellRst.txt";
  string fluidout="fluidRst.txt";

  LB channel;
  channel.readInput(fin);
  channel.init();
  channel.printInfor();
  channel.writeGeometry(fgeom);
  
  channel.writeGeometry(fgeom);
/*  Cell rbc;
  rbc.readInput(cin);
  rbc.init();
  rbc.nondimension(*channel.pUnits);
  //rbc.output(out);

  //LB a;
  //a.readInput(in);
  IBM cellInChanl(&channel,&rbc);
  */
  int nSave = 500;
  int nts =100000;
  //a.init();
  
  //a.printInfor();// this one should come after init();
  //a.output(out);
  
  clock_t begin = clock();
 /* for (int i=0;i<1200;i++)
  {
    channel.computeVelocity();
    
    //a.collideSwap();
    //a.streamSwap();
    //a.applyBC();
    //if (i%nSave ==0 )
    //  a.output(out);

    //channel.collideSwap();
    //channel.streamSwap();
    
    channel.collide();
    channel.stream();
    channel.applyBC();

  }*/
  //clock_t begin = clock();
  //cellInChanl.output(cellout); 
  
  for (int i=0;i<nts;i++){
    //channel.computeVelocity();
    //cellInChanl.interpret();
    //rbc.updateHalf();
    //rbc.computeForce();
    //cellInChanl.spread();

    channel.applyForce();
    channel.collide();
    channel.stream();
    channel.applyBC();
    //rbc.update();
    
    if (i%nSave ==0 ){
      channel.writeVelocity(fluidout);
      //rbc.output(cellout);
      cout<<"time step "<<i<<" finsished"<<endl;
    }
  }
  clock_t end = clock();
  double elapsedSecs = double(end-begin)/CLOCKS_PER_SEC;
  cout<<"time elapsed "<<elapsedSecs<<endl;

  return 0;
}
