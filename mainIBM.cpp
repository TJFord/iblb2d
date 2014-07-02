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
  //string cin="circle.txt";
  //string cin="sphere.txt";
  //string cin="chain.txt";
  //string in="inputForce.txt";
  //string fin="inputChannel.txt";
  string fin="shear.txt";
  //string fin="channel.txt";
  //string fin="vonkarman.txt";
  string fgeom="fgeom.txt";
  string cellout="cellRst.txt";
  string cellForce="cellForce.txt";
  string cellVelocity="cellVelocity.txt";
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
  
  int nSave =50000;
  int nts =2500001;//100000;
  //a.init();
  
  //a.printInfor();// this one should come after init();
  //a.output(out);
  
  //clock_t begin = clock();
  for (int i=0;i<12000;i++){
    //channel.collideSwap();
    //channel.streamSwap();
    
    channel.collide();
    channel.stream();
    channel.applyBC();

  }
     // channel.writeVelocity(fluidout);
  clock_t begin = clock();
  //cellInChanl.output(cellout); 
  
  for (int i=0;i<nts;i++){
    //---compute fluid velocity and interpret velocity---//
    //channel.computeVelocity();
    cellInChanl.interpret();
    //---update temporary position at half time step---//
    rbc.updateHalf();
    //---compute solid force based on temporary position---// 
    rbc.computeForce();
    //rbc.computeReference();
    //rbc.computeRigidForce();
    
    //---spread force to fluid---// 
    cellInChanl.spread();

    //---LB fluid solver---// 
    channel.applyForce();
    channel.collide();
    channel.stream();
    channel.applyBC();
    
    //---compute fluid velocity and interpret velocity after force spreading---//
    //channel.computeVelocity();
    cellInChanl.interpret();
    //---update position at a full time step---//
    rbc.update();
    
    if (i%nSave ==0 ){
      channel.writeVelocity(fluidout);
      channel.writeForce(fluidForce);
      rbc.writeGeometry(cellout);
      rbc.writeForce(cellForce);
      rbc.writeVelocity(cellVelocity);
      cout<<"time step "<<i<<" finsished"<<endl;
      cout<<"area "<<rbc.computeArea()/rbc.A0<<endl;
      cellInChanl.writeLog(log,i);
    }
  }
  clock_t end = clock();
  double elapsedSecs = double(end-begin)/CLOCKS_PER_SEC;
  cout<<"time elapsed "<<elapsedSecs<<endl;

  return 0;
}
