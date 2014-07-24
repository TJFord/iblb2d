#ifndef SOLID_2D_H
#define SOLID_2D_H

#include "units.h"

struct Bond{
  int size;
  int *first;
  int *next;
  double *L0;
public:
  Bond(){size=0;
    first=NULL;next=NULL;L0=NULL;
  }
  Bond(int size_){
    size = size_;
    first = new int[size];
    next = new int[size];
    L0 = new double[size];
  }
  Bond& operator=(const Bond& rhs){
    if (this !=&rhs){
      delete [] first;
      delete [] next;
      delete [] L0;
      size = rhs.size;
      first = new int[size];
      next = new int[size];
      L0 = new double[size];
      for (int i=0;i<size;i++){
        first[i]=rhs.first[i];
        next[i]=rhs.next[i];
        L0[i]=rhs.L0[i];
      }
    }
    return *this;
  }
  ~Bond(){
    delete [] first;
    delete [] next;
    delete [] L0;
  }
};
struct Angle{
  int size;
  int *left;
  int *middle;
  int *right;
  double *ang0;
public:
  Angle(){size=0;
    left=NULL;middle=NULL;right=NULL;ang0=NULL;
  }
  Angle(int size_){
    size = size_;
    left = new int[size];
    middle = new int[size];
    right = new int[size];
    ang0 = new double[size];
  }
  Angle& operator=(const Angle& rhs){
    if (this !=&rhs){
      delete [] left;
      delete [] middle;
      delete [] right;
      delete [] ang0;
      size = rhs.size;
      left = new int[size];
      middle = new int[size];
      right = new int[size];
      ang0 = new double[size];
      for (int i=0;i<size;i++){
        left[i]=rhs.left[i];
        middle[i]=rhs.middle[i];
        right[i]=rhs.right[i];
        ang0[i]=rhs.ang0[i];
      }
    }
    return *this;
  }
  ~Angle(){
    delete [] left;
    delete [] right;
    delete [] middle;
    delete [] ang0;
  }
};

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

