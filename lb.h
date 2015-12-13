#ifndef LB_2D_H
#define LB_2D_H

# include <iostream>
# include <fstream>
# include <iomanip>
//# include <string>

# include "units.h"
# include "boundary.h"

//using namespace std;

class LB;

typedef void (*funType)(double* f, void* selfData);
typedef void (LB::*collType)(int id, void *selfData);
typedef struct{
  funType dynf;
  void* selfData;
} Dynamics;

class LB
{
public:
  int lx,ly;//bounding box
  int *xy2idx;//semi direct address index
  int *nbList;//neighbor list
  int nf; // total fluid nodes
  int nt; // all nodes including buffer
  double *f;//density distribution
  double *ft;//copy of f
  int *coor;
  double *v;//velocity
  double g[2];//gravity
  std::string collisionScheme;

  Units *pUnits;

  double omega; // relaxation = 1/tau;
  static const int d;//dimension
  static const int q;//discrete velocity component
  static const double cs2; // square of sound speed
  static const int c[9][2]; //velocity vector
  static const double w[9];// weight factor
  // reserved for regularized scheme
  double Qxx[9], Qxy[9], Qyy[9];//Q tensor for regularized method
  double feq[9], fNeq[9];
  double neqPixx, neqPixy, neqPiyy;

  //double F[9]; // force term for one lattice
  int nts;//total time steps
  int ntsOut;//output timestep

  int nbc;//number of boundary conditions

  Dynamics *pBC;

  collType collisionFun;

  int nv,no,nb;//number of velocity, open boundaries;
  velData *pVel; 
  bbData *pBB;
  openData *pOpen;
public:
  double *force;
  
public:
  //creator
  LB();
  ~LB();
  LB& operator=(const LB& rhs);
  void readInput(const std::string filename); 
  void reReadVelocity(const std::string filename); //restart.
  // manipulator
  void init();//initialize density with 0 velocity
  void stream();
  void collide();
  void streamSwap();
  void collideSwap();
  void applyBC();
  void applyForce();
  //void applyForce(int id, double fx, double fy);
  void bgk(int id, void* selfData);
  void regularized(int id, void* selfData);
  void stokes(int id, void* selfData);
  
  double computeEqStokes(int idx, double rho,double ux, double uy);
  double computeEquilibrium(int idx, double rho,double ux, double uy,double uSqr);
  // accessor  
  void computeMacros(int id, double *rho, double * ux, double *uy);
  void computeVelocity();
  void writeVelocity(const std::string out);
  void writeGeometry(const std::string out);
  void writeForce(const std::string out);
  void printInfor();
  void writeLog(const std::string out);
};

#endif
