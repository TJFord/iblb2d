# include <cmath>

# include "ibm.h"
# include "lb.h"
# include "cell.h"

# define PI 3.14159265

IBM::IBM(LB* plb_, Cell* pCell_, int nc_){
  plb = plb_;
  pCell = pCell_;
  nc = nc_; 
  for (int i=0;i<4;i++){
    w1[i]=0.;
    w2[i]=0.;
  }
}

IBM::~IBM(){
}

void IBM::ph1(double r){
  double q;
  q = sqrt(1.+4.*r*(1.-r));
  w1[0]=(3.-2.*r-q)/8;
  w1[1]=(3.-2.*r+q)/8;
  w1[2]=(1.+2.*r+q)/8;
  w1[3]=(1.+2.*r-q)/8;
  /*w1[0]=phCos(r-1);
  w1[1]=phCos(r);
  w1[2]=phCos(r+1);
  w1[3]=phCos(r+2);*/
}

void IBM::ph2(double r){
  double q;
  q = sqrt(1.+4.*r*(1.-r));
  w2[0]=(3.-2.*r-q)/8;
  w2[1]=(3.-2.*r+q)/8;
  w2[2]=(1.+2.*r+q)/8;
  w2[3]=(1.+2.*r-q)/8;
  /*w2[0]=phCos(r-1);
  w2[1]=phCos(r);
  w2[2]=phCos(r+1);
  w2[3]=phCos(r+2);*/
}

double IBM::phCos(double r){
  if (r > 2.)
    return 0.;
  else if(r<-2.)
    return 0.;
  else
    return 0.25*(1+cos(PI*r*0.5));
}

void IBM::printInfor(){
  for (int i=0;i<4;i++)
    std::cout<<w1[i]<<" "<<std::endl;
}
  
void IBM::output(const std::string filename){
  using namespace std; 
  std::ofstream out(filename.c_str(),std::ios::out | std::ios::app);
    if (out.is_open()){
      out<<setw(10)<<"ks"<<setw(10)<<pCell->ks<<endl;
      out<<setw(10)<<"kb"<<setw(10)<<pCell->kb<<endl;
      out<<setw(10)<<"kp"<<setw(10)<<pCell->kp<<endl;
      out<<setw(10)<<"node"<<endl;
      out<<setw(10)<<pCell->nn<<endl;
      for (int i=0;i<pCell->nn;i++)
        out<<setw(10)<<pCell->x[2*i]<<setw(10)<<pCell->x[2*i+1]<<endl;
      
      out<<setw(10)<<"bond"<<endl;
      out<<setw(10)<<pCell->nb<<endl;
      for (int i=0;i<pCell->nb;i++)
        out<<setw(10)<<pCell->pBond->first[i]<<setw(10)<<pCell->pBond->next[i]<<setw(10)<<pCell->pBond->L0[i]<<endl;
      
      out<<setw(10)<<"angle"<<endl;
      out<<setw(10)<<pCell->na<<endl;
      for (int i=0;i<pCell->na;i++)
        out<<setw(10)<<pCell->pAngle->left[i]<<setw(10)<<pCell->pAngle->middle[i]<<setw(10)<<pCell->pAngle->right[i]<<setw(10)<<pCell->pAngle->ang0[i]<<endl;

      
    //out.close();
    }else{
      std::cout<<"cannot open cell output file"<<std::endl;
    }
   
}
/*
void IBM::interpret(){
  for (int i=0;i<nc;i++){
    int nn = pCell[i].nn;
    int nfID,lx,ly;
    double s1,s2,r1,r2;
    int i1,i2;
    double rho, ux, uy;
    lx = plb->lx;
    ly = plb->ly;
    for (int j=0;j<nn;j++){
      s1 = pCell[i].x[2*j];
      s2 = pCell[i].x[2*j+1];
      i1 =(int)floor(s1);
      i2 =(int)floor(s2);
      r1 = s1-i1;
      r2 = s2-i2;
      ph1(r1);
      ph2(r2);
      for (int k=0;k<4;k++){
        for (int m=0;m<4;m++){
          nfID = plb->IJidx[(i1-1+k)*lx+i2-1+m];
          plb->computeMacros(nfID,&rho, &ux,&uy);
          pCell[i].v[2*j] += w1[k]*w2[m]*ux;
          pCell[i].v[2*j+1] += w1[k]*w2[m]*uy;
        }
      }
    } 
  }
}*/

void IBM::interpret(){
    int nn = pCell->nn;
    int nfID,lx,ly;
    double s1,s2,r1,r2;
    int i1,i2;
    lx = plb->lx;
    ly = plb->ly;
    for(int j=0;j<nn;j++){
      pCell->v[2*j] =0.;
      pCell->v[2*j+1] =0.;
    }
    for (int j=0;j<nn;j++){
      s1 = pCell->x[2*j];
      s2 = pCell->x[2*j+1];
      i1 =(int)floor(s1);
      i2 =(int)floor(s2);
      r1 = s1-i1;
      r2 = s2-i2;
      ph1(r1);
      ph2(r2);
      for (int k=0;k<4;k++){
        for (int m=0;m<4;m++){
          //nfID = plb->IJidx[(i2-1+k)*lx+i1-1+m];
          nfID = plb->IJidx[(i2-2+k)*lx+i1-2+m];
          //plb->computeMacros(nfID,&rho, &ux,&uy);
          pCell->v[2*j] += w1[k]*w2[m]*plb->v[2*nfID];
          pCell->v[2*j+1] += w1[k]*w2[m]*plb->v[2*nfID+1];
        }
      }
    } 
  }


/*
void IBM::spread(){
  for (int i=0;i<nc;i++){
    int nn = pCell[i].nn;
    int nfID,lx,ly;
    double s1,s2,r1,r2;
    int i1,i2;
    double rho, ux, uy;
    lx = plb->lx;
    ly = plb->ly;
    for (int j=0;j<nn;j++){
      s1 = pCell[i].xhalf[2*j];
      s2 = pCell[i].xhalf[2*j+1];
      i1 =(int)floor(s1);
      i2 =(int)floor(s2);
      r1 = s1-i1;
      r2 = s2-i2;
      ph1(r1);
      ph2(r2);
      for (int k=0;i<4;k++){
        for (int m=0;m<4;m++){
          nfID = plb->IJidx[(i1-1+k)*lx+i2-1+m];
          plb->computeMacros(nfID,&rho, &ux,&uy);
          pCell[i].computeForce();
          plb->force[2*nfID] +=w1[k]*w2[m]*pCell[i].force[2*j];
          plb->force[2*nfID+1] +=w1[k]*w2[m]*pCell[i].force[2*j+1];
        }
      }
    }
  }
}*/

void IBM::spread(){
    int nn = pCell->nn;
    int nfID,lx,ly;
    double s1,s2,r1,r2;
    int i1,i2;
   // double rho, ux, uy;
    lx = plb->lx;
    ly = plb->ly;
    int nf = plb->nf;
    for (int i=0;i<nf;i++){
      plb->force[2*i]=0.;
      plb->force[2*i+1]=0.;
    }
    for (int j=0;j<nn;j++){
      s1 = pCell->x[2*j];
      s2 = pCell->x[2*j+1];
      i1 =(int)floor(s1);
      i2 =(int)floor(s2);
      r1 = s1-i1;
      r2 = s2-i2;
     // std::cout<<"x "<<s1<<" i1"<<i1<<" left"<<i1-1<<" rgt"<<i1+2<<std::endl;
      ph1(r1);
      ph2(r2);
      for (int k=0;k<4;k++){
        for (int m=0;m<4;m++){
          //nfID = plb->IJidx[(i2-1+k)*lx+i1-1+m];
          nfID = plb->IJidx[(i2-2+k)*lx+i1-2+m];
          if (nfID>plb->nt){ 
            std::cout<<"error! spread to node out of fluid domain"<<std::endl;
            std::cout<<"node id "<<nfID<<"at x="<<i1-2+m<<" y="<<i2-2+k<<std::endl;
          }
         // plb->computeMacros(nfID,&rho, &ux,&uy);
          //pCell->computeForce();
          plb->force[2*nfID] +=w1[k]*w2[m]*pCell->force[2*j];
          plb->force[2*nfID+1] +=w1[k]*w2[m]*pCell->force[2*j+1];
        }
      }
    }
  }

void IBM::periodic(){
  double radius;
  radius =1.5*pCell->radius;

  if(pCell->xc[0] > plb->lx - radius){
    for (int i=0;i<pCell->nn;i++){
      pCell->x[2*i]= pCell->x[2*i]-(plb->lx-2*radius);
    }
  }
  
}

void IBM::moveSolidTo(double x, double y){
  pCell->computeReference();
  double dx = x - pCell->xc[0];
  double dy = y - pCell->xc[1];
  for (int i=0;i<pCell->nn;i++){
    pCell->x[2*i] += dx;
    pCell->x[2*i+1] += dy;
  }
}
