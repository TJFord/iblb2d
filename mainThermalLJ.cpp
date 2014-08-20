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
  //fstream in("input.txt",ios::in);
  //string in="input.txt";
  //string out="rst.txt";
  //string cin="cellInput.txt";
  string cin="MultiCells.txt";
  string win="MultiWorms.txt";
  //string cin="circle.txt";
  //string cin="sphere.txt";
  //string cin="chain.txt";
  //string in="inputForce.txt";
  //string fin="inputChannel.txt";
  //string fin="shear.txt";
  //string fin="Velchannel.txt";
  string fin="noFlowChannel.txt";
  //string fin="channel.txt";
  //string fin="vonkarman.txt";
  string fgeom="fgeom.txt";
  string cellout="cellRst.txt";
  string cellForce="cellForce.txt";
  string cellVelocity="cellVelocity.txt";
  string wormout="chainRst.txt";
  string wormForce="chainForce.txt";
  string fluidout="fluidRst.txt";
  string fluidForce="fluidForce.txt";
  string log="Log.txt";

  LB channel;
  channel.readInput(fin);
  channel.init();
  channel.printInfor();
  channel.writeGeometry(fgeom);
  channel.writeLog(log);
 /* 
  Cell rbc;
  rbc.readInput(cin);
  rbc.init();
  rbc.nondimension(*channel.pUnits);
  //rbc.output(out);
  rbc.writeLog(log);*/
  
  Chain worm;
  worm.readInput(win);
  worm.init();
  worm.nondimension(*channel.pUnits);
  //rbc.output(out);
  worm.writeLog(log);
  

  //IBM cellInChanl(&channel,&rbc);
  IBM wormInChanl(&channel,&worm);

  worm.initLJ();
  
  int nSave =5000;
  int nts =250001;//100000;
  //a.init();
  
  //a.printInfor();// this one should come after init();
  //a.output(out);
  
  //clock_t begin = clock();
  /*
  for (int i=0;i<12000;i++){
    //channel.collideSwap();
    //channel.streamSwap();
    
    channel.collide();
    channel.stream();
    channel.applyBC();

  }*/
     // channel.writeVelocity(fluidout);
  clock_t begin = clock();
  //cellInChanl.output(cellout); 
  
  for (int i=0;i<nts;i++){
    //---compute fluid velocity and interpret velocity---//
    //channel.computeVelocity();
    //cellInChanl.interpret();
    wormInChanl.interpret();
    //---update temporary position at half time step---//
    //rbc.updateHalf();
    worm.updateHalf();
    //---compute solid force based on temporary position---// 
    //rbc.computeForce();
    worm.computeForce();
    //rbc.computeReference();
    //rbc.computeRigidForce();
    
    //---spread force to fluid---// 
    //cellInChanl.spread();
    wormInChanl.spread();


    //---LB fluid solver---// 
    channel.applyForce();
    channel.collide();
    channel.stream();
    channel.applyBC();
    
    //---compute fluid velocity and interpret velocity after force spreading---//
    //channel.computeVelocity();
    //cellInChanl.interpret();
    wormInChanl.interpret();
    //---update position at a full time step---//
    //rbc.update();
    worm.update();
    //worm.thermalFluctuation();
    
    if (i%nSave ==0 ){
      channel.writeVelocity(fluidout);
      channel.writeForce(fluidForce);
      //rbc.writeGeometry(cellout);
      //rbc.writeForce(cellForce);
      //rbc.writeVelocity(cellVelocity);
      worm.writeGeometry(wormout);
      worm.writeForce(wormForce);
      cout<<"time step "<<i<<" finsished"<<endl;
      //cout<<"area "<<rbc.computeArea()/rbc.A0<<endl;
      //cellInChanl.writeLog(log,i);
      wormInChanl.writeLog(log,i);
    }
  }
  clock_t end = clock();
  double elapsedSecs = double(end-begin)/CLOCKS_PER_SEC;
  cout<<"time elapsed "<<elapsedSecs<<endl;

  return 0;
}
