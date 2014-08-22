#ifndef CHAIN_2D_H
#define CHAIN_2D_H
//# include <random>
# include "random_mars.h"
# include "units.h"
# include "solid.h"
# include "cell.h"

#define EMPTY -1
/*
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
};*/

class Chain: public Solid{
public:
  //Solid *pCell;
  Cell *pCell;
  //double *x;
  double *xtmp;
  //double *v;
  //double *force;
  
  Bond *pBond;
  Angle *pAngle;

  double rho;
  double diffusionCoef;
  double diffDist;
  double ks;
  double kb;
  double kp;
  double kBT; //thermal energy
  //double damp;// Brownian dynamics dampping,in time units, 100 
  double m;
  double g;
  double diameter;

  //double xc[2];

  //int nn, nb, na;
  int nb, na;
  int ns;//# of cells
  int nForOne;// # of nodes for each solid object
  
  int periodicX, periodicY;
  //double halfX, halfY;
  //random number generator:
  // this one needs c++11 suport, include headfile <random>
  //std::default_random_engine generator;
  //std::normal_distribution<double> ndist;
  //portable random file from lammps
  RanMars *random;
  int seed;

  //linked list cells for pairwise potential
  int *lscl;//list of cell link
  int *head; //head for each cell
  double rCut,rrCut;//cut of distance and its square
  double sig2,sig6;
  int cx,cy;//cell number in x,y direction
  double epsilon, sigma;
  double zeta;//friction coefficient

public:
  // creator
  Chain();
  ~Chain();
  Chain& operator=(const Chain& rhs);
  void readInput(const std::string filename);
  void reReadPosition(const std::string filename);
  void init();
  void setCells(Cell * pCell_);
  // manipulator
  void update();
  void updateHalf();
  void nondimension(const Units&);
  void computeEquilibrium();
  void computeForce();
  void bondHarmonicForce();
  void angleBendForce();
  void velocityVerletIntegration();
  void moveTo(double x,double y);
  void randomForce();
  void randomDisplacement();
  void penetrationRemoval();
  void initLJ();
  void buildLinkList();
  void pairWiseInteraction();
  void LJForce(int i,int j);
  int particleInsideCell(double xp, double yp, int idx);
  void moveOutside(double xp, double yp, int idc, int idp);
  // acessor
  void writeGeometry(const std::string filename);
  void writeForce(const std::string filename);
  void writeVelocity(const std::string filename);
  void writeLog(const std::string filename);
};
#endif
