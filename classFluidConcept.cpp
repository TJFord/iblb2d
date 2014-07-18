class LB
{
  int lx,ly;//bounding box
  int *IJidx;//semi direct address index
  int *nbList;//neighbor list
  int nf; // total fluid nodes
  double *f;//density distribution
  int *coor;

  Units *pUnits;

  double omega; // relaxation = 1/tau;
  static const int d;//dimension
  static const int q;//discrete velocity component
  static const double cs2; // square of sound speed

  static const int c[9][2]; //velocity vector
  static const double w[9];// weight factor

  int nbc;//number of boundary conditions

  Dynamics *pBC;

  collType collisionFun;

  int nv,no,nb;//number of velocity, open boundaries;
  velData *pVel; 
  bbData *pBB;
  openData *pOpen;
public:
  //creator
  LB();
  ~LB();
  LB& operator=(const LB& rhs);
  void readInput(const std::string filename); //only read nn, nlist, etc.
  // manipulator
  void init();//initialize density with 0 velocity
  void stream();
  void collide();
  void streamSwap();
  void collideSwap();
  void applyBC();
  //void applyForce(int id, double fx, double fy);
  void bgk(int id, void* selfData);
  double computeEquilibrium(int idx, double rho,double ux, double uy,double uSqr);
  // accessor  
  void computeMacros(int id, double *rho, double * ux, double *uy);
  void output(const std::string out);
  void printInfor();
};

class Solid
{
  double *x;
  double *v;
  double *a;
  double *x0;
  double *force;
  double xc[3];//center position for periodic BC
  int nn, nb, na;
public:
  virtual void readInput(const std::string filename)=0;
  virtual void computeForce()=0;
  virtual void init()=0;
  virtual void nondimension(const Units&);
  virtual void updateHalf()=0;
  virtual void update()=0;
}

class IBM
{
  Fluid *pFluid;
  Solid *pSolid;
public:
  void interpret(pFluid, pSolid, *velocity){
    for i in pSolid.node
      x1 = x(i)-1:x(i)+1;
      y1 = y(i)-1:y(i)+1;
      computeMacros(IJidx(x1,y1),rho, ux,uy);
      pSolid->ux[i] += ux*weight;
      pSolid->uy[i] += uy*weight;
  }
  void spread(pFluid, pSolid){
    for i in pSolid.node
      x1 = x(i)-1:x(i)+1;
      y1 = y(i)-1:y(i)+1;
      pFluid->force_x[IJidx(x1,y1)] += pSolid->force[x]*weight;
      pFluid->force_y[IJidx(x1,y1)] += pSolid->force[y]*weight; 
  }
  void ph1();
  void ph2();
  void ph3();
}
main()
{
  int ts;
  IBM ibm;
  ibm.pFluid->readInput();
  ibm.pSolid->readInput();

  for (int i=0;i<ts;i++)
  {
    ibm.interpret()  
    if (i%tsave == 0)
      saveResults();
  }	
}
