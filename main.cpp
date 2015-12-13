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
#include "chain.h"
#include "ibm.h"

int main(int argc, char *argv[])
{
  //string cin="MultiCells.txt";
  //string win="particles.txt";
  string cin="sphere.txt";
  //string fin="Velchannel.txt";
  string fin="channel.txt";
  //string fin="vonkarman.txt";
  string fgeom="fgeom.txt";
  string cellout="cellRst.txt";
  string cellForce="cellForce.txt";
  string cellVelocity="cellVelocity.txt";
  //string wormout="chainRst.txt";
  //string wormForce="chainForce.txt";
  //string wormVelocity="chainVelocity.txt";
  string fluidout="fluidRst.txt";
  string fluidForce="fluidForce.txt";
  string log="Log.txt";

  LB channel;
  channel.readInput(fin);
  channel.init();
  channel.printInfor();
  channel.writeGeometry(fgeom);
  channel.writeLog(log);
  
  Cell rbc;
  rbc.readInput(cin);
  rbc.init();
  rbc.nondimension(*channel.pUnits);
  //rbc.output(out);
  rbc.writeLog(log);
  
  IBM cellInChanl(&channel,&rbc);
 
  int nSave =channel.ntsOut;//250000;//100000;//5000;
  int nts =channel.nts;//12500001;//12500001;//5000001;//250001;//5000001;//250001;//100000;
  
/*  
  for (int i=0;i<12000;i++){
    //channel.collideSwap();
    //channel.streamSwap();
    
    channel.collide();
    channel.stream();
    channel.applyBC();

  }*/
 
  clock_t begin = clock();
  
  for (int i=0;i<nts;i++){
    //---compute fluid velocity and interpret velocity---//
    cellInChanl.interpret();
    //---update temporary position at half time step---//
    rbc.updateHalf();
    //---compute solid force based on temporary position---// 
    rbc.computeForce();
    //---spread and apply to fluid---// 
    cellInChanl.spread();
    channel.applyForce();//spread will overwrite plb->force=0

    //---LB fluid solver---// 
    channel.collide();
    channel.stream();
    channel.applyBC();
    
    //---compute fluid velocity and interpret velocity after force spreading---//
    cellInChanl.interpret();
    //---update position at a full time step---//
    rbc.update();

    if (i%nSave ==0 ){
      channel.writeVelocity(fluidout);
      channel.writeForce(fluidForce);//only writes the last spread force
      rbc.writeGeometry(cellout);
      rbc.writeForce(cellForce);
      rbc.writeVelocity(cellVelocity);
      cout<<"time step "<<i<<" finsished"<<endl;
      cellInChanl.writeLog(log,i);
    }
  }
  clock_t end = clock();
  double elapsedSecs = double(end-begin)/CLOCKS_PER_SEC;
  cout<<"time elapsed "<<elapsedSecs<<" seconds"<<endl;
  cellInChanl.writeLog(log,int(elapsedSecs));

  return 0;
}
