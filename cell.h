#ifndef CELL_2D_H
#define CELL_2D_H

# include "units.h"

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
class Cell{
public:
  double *x;
  double *xtmp;
  double *xR;//reference position for rigids
  double *v;
  //double *a;
  double *x0;
  double *force;
  
  Bond *pBond;
  Angle *pAngle;

  double rho;
  double ks;
  double kb;
  double kp;

  double m;
  double g;

  double radius;//reserved for rigid spheres
  double angRef;
  double xc0[2];
  double xc[2];

  int nn, nb, na;
public:
  // creator
  Cell();
  ~Cell();
  Cell& operator=(const Cell& rhs);
  void readInput(const std::string filename);
  void init();
  // manipulator
  void update();
  void updateHalf();
  void nondimension(const Units&);
  void computeForce();
  void bondHarmonicForce();
  void angleBendForce();
  void velocityVerletIntegration();
  void computeReference();
  void computeRigidForce();
  // acessor
  void writeGeometry(const std::string filename);
  void writeReferenceGeometry(const std::string filename);
  void writeForce(const std::string filename);
};
#endif
