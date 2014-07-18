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
  //string cin="cellInput.txt";
  //string cin="circle.txt";
  string cin="sphere.txt";
  //string cin="chain.txt";
  //string in="inputForce.txt";
  //string fin="inputChannel.txt";
  //string fin="shear.txt";
  string fin="channel.txt";
  //string fin="Velchannel.txt";
  //string fin="vonkarman.txt";
  string fgeom="fgeom.txt";
  string cellout="cellRst.txt";
  string cellRef="cellRef.txt";
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
  
  int nSave =10000;//10000;//2000;
  int nts =500001;//500000;

  //cellInChanl.moveSolidTo(395,50);
  //a.init();
  //a.printInfor();// this one should come after init();
  //a.output(out);
  
  //clock_t begin = clock();
  /*
  for (int i=0;i<1200;i++)
  {
    //channel.computeVelocity();
    
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
  clock_t begin = clock();
  //cellInChanl.output(cellout); 
  
  for (int i=0;i<nts;i++){
    cellInChanl.interpret();
    rbc.updateHalf();
    //rbc.computeForce();
    rbc.computeReference();
    rbc.computeRigidForce();

    cellInChanl.spread();

    channel.applyForce();
    channel.collide();
    channel.stream();
    channel.applyBC();
    //---compute fluid velocity and interpret velocity after force spreading---//
    cellInChanl.interpret();
    //---update position at a full time step---//
    rbc.update();
    
    cellInChanl.periodic();

    if (i%nSave ==0 ){
      channel.writeVelocity(fluidout);
      channel.writeForce(fluidForce);
      rbc.writeGeometry(cellout);
      rbc.writeForce(cellForce);
      rbc.writeReferenceGeometry(cellRef);
      rbc.writeVelocity(cellVelocity);
      cout<<"time step "<<i<<" finished"<<endl;
      //cout<<"area "<<rbc.computeArea()/rbc.A0<<endl;
      cellInChanl.writeLog(log,i);
    }
  }
  clock_t end = clock();
  double elapsedSecs = double(end-begin)/CLOCKS_PER_SEC;
  cout<<"time elapsed "<<elapsedSecs<<endl;

  return 0;
}
