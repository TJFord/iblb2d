#ifndef SOLID_2D_H
#define SOLID_2D_H

#include "units.h"

class Solid{
  public:
    int nn;//totoal number of solid nodes;
    int ns;//# of solid strucutes, e.g., cells
    int lx,ly; // fluid domain size
    double *x;
    double *v;
    //double xc[2];
    double *xc;
    double *force;
  public:
    Solid();
    virtual ~Solid();
    /*virtual void readInput(const std::string filename)=0;
    virtual void computeForce()=0;
    virtual void init()=0;
    virtual void nondimension(const Units&)=0;
    virtual void updateHalf()=0;
    virtual void update()=0;
*/

};
#endif

