# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>
# include "boundary.h"

static const int oppositeOf[9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };

void bounceBack(double* f, void* selfData)
{
  bbData *data = (bbData*)selfData;
  double fTmp[9];
  int id,st,size;
  int iPop;
  size = data->size;
  for(int i=0;i<size;i++){
    id = data->node[i];
    st = id*9;
    for (iPop=0; iPop<9; ++iPop){
      fTmp[iPop] = f[st + oppositeOf[iPop]];
    }
    for (iPop=0; iPop<9; ++iPop){
      f[st + iPop] = fTmp[iPop];
    }
  }
}

inline static double leftRho(double* f, int id, double ux) {
  int st = id*9;
  return 1./(1.-ux) * (
    f[st] + f[st+2] + f[st+4] + 2*(f[st+3]+f[st+6]+f[st+7])
  );
}

inline static double rightRho(double* f, int id, double ux) {
  int st = id*9;
  return 1./(1.+ux) * (
    f[st] + f[st+2] + f[st+4] + 2*(f[st+1]+f[st+5]+f[st+8])
  );
}

inline static double upperRho(double* f, int id, double uy) {
  int st = id*9;
  return 1./(1.+uy) * (
    f[st] + f[st+3] + f[st+1] + 2*(f[st+6]+f[st+2]+f[st+5])
  );
}

  /* Compute density on wall from bulk information on
     lower boundary. */

inline static double lowerRho(double* f, int id, double uy) {
  int st = id*9;
  return 1./(1.-uy) * (
    f[st] + f[st+3] + f[st+1] + 2*(f[st+8]+f[st+4]+f[st+7])
  );
}

inline static void completeLeft(double* f,int id,
                                double ux, double uy, double rho)
{
  int st = id*9;
  f[st+5] = f[st+7] + 0.5 *(f[st+4]-f[st+2])
                    + rho*ux/6. + rho*uy/2.;
  f[st+8] = f[st+6] + 0.5 *(f[st+2]-f[st+4])
                    + rho*ux/6. - rho*uy/2.;
  f[st+1] = f[st+3] + 2./3.*rho*ux;
}

inline static void completeRight(double* f,int id,
                                double ux, double uy, double rho)
{
  int st = id*9;
  f[st+6] = f[st+8] + 0.5 *(f[st+4]-f[st+2]) 
                    - rho*ux/6 + rho*uy/2.;
  f[st+7] = f[st+5] + 0.5 *(f[st+2]-f[st+4]) 
                    - rho*ux/6 - rho*uy/2.;
  f[st+3] = f[st+1] - 2./3.*rho*ux;
}

inline static void completeUpper(double* f,int id,
                                 double ux, double uy, double rho)
{
  int st = id*9;
  f[st+7] = f[st+5] + 0.5 * (f[st+1]-f[st+3])
                    - rho*uy/6 - rho*ux/2.;
  f[st+8] = f[st+6] + 0.5 *(f[st+3]-f[st+1]) 
                    - rho*uy/6 + rho*ux/2.;
  f[st+4] = f[st+2] - 2./3.*rho*uy;
}

inline static void completeLower(double* f, int id,
                                 double ux, double uy, double rho)
{
  int st = id*9;
  f[st+6] = f[st+8] + 0.5 *(f[st+1]-f[st+3]) 
                    + rho*uy/6 - rho*ux/2.;
  f[st+5] = f[st+7] + 0.5 *(f[st+3]-f[st+1]) 
                    + rho*uy/6 + rho*ux/2.;
  f[st+2] = f[st+4] + 2./3.*rho*uy;
}


  // ZouHe velocity boundaries on upper, lower, left and right
  //   boundaries
void leftZouHe(double* f, void* selfData) {
  velData* data = (velData*) selfData;
  double rho;
  int id,size;
  double ux,uy;
  size = data->size;
  for (int i=0;i<size;i++){
    id = data->node[i];
    ux = data->ux[i];
    uy = data->uy[i];
    rho = leftRho(f,id,ux);
    completeLeft(f,id,ux,uy,rho);
  }
}

// right zou he velocity boundary 
void rightZouHe(double* f, void* selfData) {
  velData* data = (velData*) selfData;
  double rho;
  int id,size;
  double ux,uy;
  size = data->size;
  for (int i=0;i<size;i++){
    id = data->node[i];
    ux = data->ux[i];
    uy = data->uy[i];
    rho = rightRho(f,id,ux);
    completeRight(f,id,ux,uy,rho);
  }
}

void upperZouHe(double* f, void* selfData) {
  velData* data = (velData*) selfData;
  double rho;
  int id,size;
  double ux,uy;
  size = data->size;
  for (int i=0;i<size;i++){
    id = data->node[i];
    ux = data->ux[i];
    uy = data->uy[i];
    rho = upperRho(f,id,uy);
    completeUpper(f,id,ux,uy,rho);
  }
}

void lowerZouHe(double* f, void* selfData) {
  velData* data = (velData*) selfData;
  double rho;
  int id,size;
  double ux,uy;
  size = data->size;
  for (int i=0;i<size;i++){
    id = data->node[i];
    ux = data->ux[i];
    uy = data->uy[i];
    rho = lowerRho(f,id,uy);
    completeLower(f,id,ux,uy,rho);
  }
}

void rightOpen (double* f, void* selfData)
{
  openData *data = (openData*)selfData;
  int id,st,size,next,tgt;
  size = data->size;
  for (int i=0;i<size;i++){
    id = data->node[i];
    next = data->nodeNext[i];
    st = id*9;
    tgt = next*9;
    f[st+3] = f[tgt+3];
    f[st+6] = f[tgt+6];
    f[st+7] = f[tgt+7];
  }
/*  int st = id*9;
  int *leftID =(int*)selfData;
  //std::cout<<*leftID<<std::endl;
  int tgt = *leftID;
  f[st+3] = f[tgt*9+3];
  f[st+6] = f[tgt*9+6];
  f[st+7] = f[tgt*9+7];*/
/*	
  f[st+1] = f[tgt*9+1];
  f[st+2] = f[tgt*9+2];
  f[st+4] = f[tgt*9+4];
  f[st+5] = f[tgt*9+5];
  f[st+8] = f[tgt*9+8];
  f[st] = f[tgt*9];*/
}

void leftOpen(double *f, void* selfData){}
void upperOpen(double *f, void* selfData){}
void lowerOpen(double *f, void* selfData){}


/*
void velocityZouHe(double* f, int id, int *n, void* selfData)
//void velocityZouHe(double *f, int id, void*selfData)
{
  velData* data = (velData*) selfData;
  double ux = data->ux;
  double uy = data->uy;
  int inplane[3]={0,0,0};
  int leaving[3]={0,0,0};
  int unknown[3]={0,0,0};
  int tmp(0),tmp2(0),tmp3(0);
  for (int i=0;i<q;i++)
  {
      int dotpt=c[i][0]*n[0]+c[i][1]*n[1];
      if (dotpt==0)
      {
              inplane[tmp]=i;
              tmp++;
      }
      else if (dotpt > 0)
      {
              leaving[tmp2]=i;
              tmp2++;
      }
      else 
      {
              unknown[tmp3]=i;
              tmp3++;
      }
    }
    int st = id*q;
    if (n[0]==0)//along y direction
    double rho = 1./(1-uy)*(f[st+inplane[0]]+f[st+inplane[1]]+f[st+inplane[2]]
                            +2*(f[st+leaving[0]]+f[st+leaving[1]]+f[st+leaving[2]]));
// double rho = rightRho(f,id, data->ux);
// completeRight(f,id, data->ux, data->uy, rho);
                            //not correct from below
    f[st+unknown[0]] = f[st+oppositeOf[unknown[0]]] + 0.5 *(f[st+4]-f[st+2]) 
                  - rho*ux/6 + rho*uy/2.;
f[st+7] = f[st+5] + 0.5 *(f[st+2]-f[st+4]) 
                  - rho*ux/6 - rho*uy/2.;
f[st+3] = f[st+1] - 2./3.*rho*ux;

    else if (n[1]==0)
    else
            cerr<<"can not handle nonorthogonal boundary"<<endl;
            

}
*/
