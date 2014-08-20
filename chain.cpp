# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>
# include <cmath>

# include "chain.h"

Chain::Chain(){
  x=NULL;
  xc=NULL;
  xtmp=NULL;
  v=NULL;
  force=NULL;
  pBond=NULL;
  pAngle=NULL;
  random=NULL;

  lscl=NULL;
  head=NULL;

  rho=0.;
  m=0.;
  ks=0.;
  kb=0.;
  kp=0.;
  g=9.8;
  diffusionCoef=0;
  diffDist=0.;

  periodicX=0;
  periodicY=0;
}

Chain& Chain::operator=(const Chain& rhs){
  if(this!=&rhs){
    delete [] x;
    delete [] xc;
    delete [] xtmp;
    delete [] v;
    delete [] force;
    delete pBond;
    delete pAngle;
    delete random;

    delete [] lscl;
    delete [] head;
    
    nn=rhs.nn; nb=rhs.nb; na=rhs.na;
    ns=rhs.ns;
    rho=rhs.rho; m=rhs.m;
    ks=rhs.ks; kb=rhs.kb;
    kp=rhs.kp;
    g=rhs.g;
    lx=rhs.lx; ly=rhs.ly;
    diffusionCoef=rhs.diffusionCoef;
    diffDist=rhs.diffDist;
    
    rCut=rhs.rCut;
    cx=rhs.cx;
    cy=rhs.cy;

    periodicX = rhs.periodicX;
    periodicY = rhs.periodicY;
    
    lscl=new int[nn]; //have to rebuild each time
    head=new int[cx*cy];

    x=new double[nn*2];
    xtmp=new double[nn*2];
    v=new double[nn*2];
    force=new double[nn*2];
    for (int i=0;i<nn;i++){
      x[2*i]=rhs.x[2*i];
      x[2*i+1]=rhs.x[2*i+1];
      xtmp[2*i]=rhs.xtmp[2*i];
      xtmp[2*i+1]=rhs.xtmp[2*i+1];
      v[2*i]=rhs.v[2*i];
      v[2*i+1]=rhs.v[2*i+1];
      force[2*i]=rhs.force[2*i];
      force[2*i+1]=rhs.force[2*i+1];
    }
    xc=new double[ns*2];
    for (int i=0;i<ns;i++){
      xc[2*i]=rhs.xc[2*i];
      xc[2*i+1]=rhs.xc[2*i+1];
    }
    int size;
    pBond = new Bond;
    size = rhs.pBond->size;
    *pBond = Bond(size);
    for (int i=0;i<pBond->size;i++){
      pBond->first[i]=rhs.pBond->first[i];
      pBond->next[i]=rhs.pBond->next[i];
      pBond->L0[i]=rhs.pBond->L0[i];
    }

    pAngle = new Angle;
    size = rhs.pAngle->size;
    *pAngle = Angle(size);
    for (int i=0;i<pAngle->size;i++){
      pAngle->left[i]=rhs.pAngle->left[i];
      pAngle->middle[i]=rhs.pAngle->middle[i];
      pAngle->right[i]=rhs.pAngle->right[i];
      pAngle->ang0[i]=rhs.pAngle->ang0[i];
    }
  }
  return *this;
}

void Chain::readInput(const std::string filename){
  std::ifstream in(filename.c_str(),std::ios::in);
  if (in.is_open()){
    int tmp=0;
    std::string str;
    while(in>>str){
      //in>>str;
      if(str.compare("m")==0){
        in>>m;
      }else if(str.compare("diffusion")==0){
        in>>diffusionCoef;
      }else if(str.compare("ks")==0){
        in>>ks;
      }else if(str.compare("kb")==0){
        in>>kb;
      }else if(str.compare("kp")==0){
        in>>kp;
      }else if(str.compare("ljcut")==0){
        in>>epsilon>>sigma>>rCut;
      }else if(str.compare("kBT")==0){
        in>>kBT;
      }else if(str.compare("periodic")==0){
        in>>periodicX>>periodicY;
      }else if(str.compare("ns")==0){
        in>>ns;
        xc=new double[ns*2];
      }else if(str.compare("node")==0){
        in>>nn;
        x=new double[nn*2];
        xtmp=new double[nn*2];
        v=new double[nn*2];
        force=new double[nn*2];

        for (int i=0;i<nn;i++)
          in>>x[2*i]>>x[2*i+1];
      }else if(str.compare("bond")==0){
        in>>nb;
        pBond = new Bond;
        *pBond = Bond(nb);
        for (int i=0;i<nb;i++)
          in>>pBond->first[i]>>pBond->next[i]>>pBond->L0[i];
      }else if(str.compare("angle")==0){
        in>>na;
        tmp++;
        pAngle = new Angle;
        *pAngle = Angle(na);
        for (int i=0;i<na;i++){
          in>>pAngle->left[i]>>pAngle->middle[i]>>pAngle->right[i]>>pAngle->ang0[i];
        }
      }else{
        std::cout<<"unknown keywords in cell input"<<std::endl;
        return;
      }
    }
    std::cout<<"successfully read cell input"<<std::endl;
  }else{
    std::cout<<"cannot open input file"<<std::endl;
  }
}

Chain::~Chain(){
  delete [] x;
  delete [] xc;
  delete [] xtmp;
  delete [] v;
  delete [] force;
  delete pBond;
  delete pAngle;
  delete random;
  delete [] lscl;
  delete [] head;
}

void Chain::reReadPosition(const std::string filename){
  std::ifstream in(filename.c_str(),std::ios::in);
  if (in.is_open()){
    for (int i=0;i<nn;i++)
      in>>x[2*i]>>x[2*i+1];
  }else{
    std::cout<<"cannot reread cell position"<<std::endl;
  }
}

void Chain::init(){
  for (int i=0;i<nn;i++){
    v[2*i]=0.;
    v[2*i+1]=0.;
    //std::cout<<"xhalf "<<xhalf[2*i]<<" "<<xhalf[2*i+1]<<std::endl;
  }
  nForOne=nn/ns;
  damp=100;
  /*int tmp;
  for (int i=0;i<ns;i++){
    xc[2*i]=0.;
    xc[2*i+1]=0.;
    for (int j=0;j<nForOne;j++){
      tmp = i*nForOne + j;
      xc[2*i] += x[2*tmp];
      xc[2*i+1] += x[2*tmp+1];
    }
    xc[2*i] /= nForOne;
    xc[2*i+1] /= nForOne;*/
  
/*
  A0=computeArea();
  //x[0] *= 1.5;
  //x[4] -=0.5;
  //x[1] -= 0.5;
  //x[2*4] *= 1.02;
  xc0[0]/=nn;
  xc0[1]/=nn;
  dx = x[0]-xc0[0];
  dy = x[1]-xc0[1];
  radius = sqrt(dx*dx+dy*dy);
  angRef[0] = acos(dx/radius);
  for(int i=1;i<nn;i++){
    dx=x[2*i]-xc0[0];
    dy=x[2*i+1]-xc0[1];
    angRef[i]=acos(dx/radius);
    //if (i<20)
    //  x[2*i+1] += 1;
  }*/
  //computeEquilibrium();
  random = new RanMars(12379);
}

void Chain::update(){
  for (int i=0;i<nn;i++){
    x[2*i] = xtmp[2*i]+v[2*i];
    x[2*i+1] =xtmp[2*i+1]+v[2*i+1];

    if (periodicX){
      if (x[2*i] < 0.) x[2*i] += lx;
      if (x[2*i] >= lx) x[2*i] -= lx;
    }
    if (periodicY){
      if (x[2*i+1] < 0.) x[2*i+1] += ly;
      if (x[2*i+1] >= ly) x[2*i+1] -= ly;
    }
    //bounceback at top bottom wall
    if (x[2*i+1]>ly) x[2*i+1] = 2*ly-x[2*i+1];
    if (x[2*i+1]<0) x[2*i+1] = -x[2*i+1];
  }
  /*
  int tmp;
  for (int i=0;i<ns;i++){
    xc[2*i]=0.;
    xc[2*i+1]=0.;
    for (int j=0;j<nForOne;j++){
      tmp = i*nForOne+j;
      x[2*tmp] = xtmp[2*tmp] + v[2*tmp];
      x[2*tmp+1] = xtmp[2*tmp+1] + v[2*tmp+1];
      xc[2*i] += x[2*tmp];
      xc[2*i+1] += x[2*tmp+1];
    }
    xc[2*i] /= nForOne;
    xc[2*i+1] /= nForOne;
  }*/
}

void Chain::updateHalf(){
  for (int i=0;i<nn;i++){
    xtmp[2*i] =x[2*i];
    xtmp[2*i+1]=x[2*i+1];
    x[2*i] += 0.5*v[2*i];
    x[2*i+1] += 0.5*v[2*i+1];
  }
}

void Chain::nondimension(const Units& unt){
  // ks: energy Es = 0.5*ks*dr^2
  ks *= unt.dt*unt.dt/unt.dm;
  // kb: energy units, Eb = 0.5*kb*dtheta^2
  kb *= unt.dt*unt.dt/unt.dm/unt.dx/unt.dx;
  // kp: energy Ep = 0.5*kp*dA^2,for 2d, this 
  // is not derived rigorously based on derivative
  // thus, the force is simplified as F= kp*dA;
  kp *= unt.dt*unt.dt*unt.dx/unt.dm;
  kBT *= unt.dt*unt.dt/unt.dm/unt.dx/unt.dx;
  m /= unt.dm;
  g *= unt.dt*unt.dt/unt.dx;
  diffusionCoef *= unt.dt/unt.dx/unt.dx;
  diffDist = sqrt(2*diffusionCoef);
  rCut /= unt.dx;
  sigma /= unt.dx;
  epsilon *= unt.dt*unt.dt/unt.dm/unt.dx/unt.dx;
  std::cout<<"rCut,sigma "<<rCut<<" "<<sigma<<std::endl;
  std::cout<<"diffDist"<<diffDist<<std::endl;
  std::cout<<"kBT "<<kBT<<" m "<<m<<std::endl;
  std::cout<<"ks, kb, kp "<<ks<<" "<<kb<<" "<<kp<<std::endl;
}

void Chain::computeEquilibrium(){
  double dx,dy,rsq;
    int i1,i2,i3;
    for (int i=0;i<nb;i++){
      i1=pBond->first[i];
      i2=pBond->next[i];
      dx=x[2*i1]-x[2*i2];
      dy=x[2*i1+1]-x[2*i2+1];
      rsq = dx*dx+dy*dy;
      pBond->L0[i]=sqrt(rsq);
    }
    double dx1, dy1, dx2,dy2, rsq1, rsq2, r1,r2;
    double c, s;//cosine, sine
    for (int i=0;i<na;i++){
      i1 = pAngle->left[i];
      i2 = pAngle->middle[i];
      i3 = pAngle->right[i];
      dx1 = x[2*i1]-x[2*i2];
      dy1 = x[2*i1+1]-x[2*i2+1];
      dx2 = x[2*i3]-x[2*i2];
      dy2 = x[2*i3+1]-x[2*i2+1];
      rsq1 = dx1*dx1+dy1*dy1;
      rsq2 = dx2*dx2+dy2*dy2;
      r1 = sqrt(rsq1);
      r2 = sqrt(rsq2);
      c = dx1*dx2 + dy1*dy2;
      c /= r1*r2;

      if (c>1.0) c=1.0;
      if (c<-1.0) c=-1.0;
      s = sqrt(1.0 - c*c);
      if (s<0.001) s=0.001;
      s = 1.0/s;
      pAngle->ang0[i]=acos(c);
    }
}

void Chain::bondHarmonicForce(){
  //--calculate bond force with F = ks*(r-r0)---//
  double dx,dy,rsq,r,fv,dr;
  int i1,i2;
  double halfX, halfY;
  halfX=0.5*lx;
  halfY=0.5*ly;
  for (int i=0;i<nb;i++){
    i1=pBond->first[i];
    i2=pBond->next[i];
    dx=x[2*i1]-x[2*i2];
    dy=x[2*i1+1]-x[2*i2+1];
    if (periodicX){
      if (std::abs(dx)>halfX){
        if (dx > 0.)   
          dx=x[2*i1]-x[2*i2]-lx;
        else
          dx=x[2*i1]-x[2*i2]+lx;
      }
    }
    if (periodicY){
      if (std::abs(dy)>halfY){
        if (dy > 0.)   
          dy=x[2*i1+1]-x[2*i2+1]-ly;
        else
          dy=x[2*i1+1]-x[2*i2+1]+ly;
      }
    }
    rsq = dx*dx+dy*dy;
  /*  if (i==nb-1){
      std::cout<<"x "<<i1<<" "<<xhalf[2*i1]<<" "<<xhalf[2*i1+1]<<std::endl;
      std::cout<<"x "<<xhalf[2*i2]<<" "<<xhalf[2*i2+1]<<std::endl;
      std::cout<<"rsq"<<rsq<<std::endl;
    }*/
    //std::cout<<"rsq"<<rsq<<" i1 "<<i1<<" "<<i2<<std::endl;
    r=sqrt(rsq);
    //std::cout<<"i "<<i<<" r "<<r<<" L0 "<<pBond->L0[i]<<std::endl;
    dr = r-(pBond->L0[i]);
    
    if (r>0.0)
      fv = -ks*dr/r;
    else
      fv = 0.0;

    force[2*i1] += fv*dx;
    force[2*i1+1] += fv*dy;
    force[2*i2] -= fv*dx;
    force[2*i2+1] -= fv*dy;
  }
}

void Chain::angleBendForce(){
  //--angle bending force Energy = 0.5*kb(dtheta)^2*--//
  int i1, i2, i3;
  double dx1, dy1, dx2,dy2, rsq1, rsq2, r1,r2;
  double c, s;//cosine, sine
  double dtheta,a,a11,a12,a22,f1x,f1y,f3x,f3y;
  double halfX = 0.5*lx;
  double halfY = 0.5*ly;
  for (int i=0;i<na;i++){
    i1 = pAngle->left[i];
    i2 = pAngle->middle[i];
    i3 = pAngle->right[i];
    dx1 = x[2*i1]-x[2*i2];
    dy1 = x[2*i1+1]-x[2*i2+1];
    dx2 = x[2*i3]-x[2*i2];
    dy2 = x[2*i3+1]-x[2*i2+1];
    if (periodicX){
      if (std::abs(dx1)>halfX){
        if (dx1 > 0.)   
          dx1=x[2*i1]-x[2*i2]-lx;
        else
          dx1=x[2*i1]-x[2*i2]+lx;
      }
      if (std::abs(dx2)>halfX){
        if (dx2 > 0.)   
          dx2=x[2*i3]-x[2*i2]-lx;
        else
          dx2=x[2*i3]-x[2*i2]+lx;
      }
    }
    if (periodicY){
      if (std::abs(dy1)>halfY){
        if (dy1 > 0.)   
          dy1=x[2*i1+1]-x[2*i2+1]-ly;
        else
          dy1=x[2*i1+1]-x[2*i2+1]+ly;
      }
      if (std::abs(dy2)>halfY){
        if (dy2 > 0.)   
          dy2=x[2*i3+1]-x[2*i2+1]-ly;
        else
          dy2=x[2*i3+1]-x[2*i2+1]+ly;
      }
    }
    rsq1 = dx1*dx1+dy1*dy1;
    rsq2 = dx2*dx2+dy2*dy2;
    r1 = sqrt(rsq1);
    r2 = sqrt(rsq2);
    c = dx1*dx2 + dy1*dy2;
    c /= r1*r2;

    if (c>1.0) c=1.0;
    if (c<-1.0) c=-1.0;
    s = sqrt(1.0 - c*c);
    if (s<0.001) s=0.001;
    s = 1.0/s;
    
    dtheta = acos(c)-pAngle->ang0[i];
    a = -kb*dtheta*s;
    a11 = a*c/rsq1;
    a12 = -a/(r1*r2);
    a22 = a*c/rsq2;
    f1x = a11*dx1 + a12*dx2;
    f1y = a11*dy1 + a12*dy2;
    f3x = a22*dx2 + a12*dx1;
    f3y = a22*dy2 + a12*dy1;
    //if (i==0){
    //  std::cout<<"f1,f3 "<<f1x<<" "<<f1y<<" "<<f3x<<" "<<f3y<<std::endl; 
    //}
    force[2*i1] += f1x;
    force[2*i1+1] += f1y;
    force[2*i2] -= f1x+f3x;
    force[2*i2+1] -= f1y+f3y;
    force[2*i3] += f3x;
    force[2*i3+1] += f3y;

  }
}


void Chain::computeForce(){
  for (int i=0;i<nn;i++){
    force[2*i]=0.;
    force[2*i+1]=0.;
  }
  if (nb){
    bondHarmonicForce();
    //std::cout<<"chain bond force"<<std::endl;
  }
  if (na){
    angleBendForce();
    //std::cout<<"chain angle force"<<std::endl;
  }
  if (rCut){
    buildLinkList();
    pairWiseInteraction();
    std::cout<<"chain LJ force"<<std::endl;
  }
  //std::cout<<"ks"<<ks<<std::endl;
  //for (int i=0;i<nn;i++)
  //  std::cout<<i<<" force "<<force[2*i]<<" "<<force[2*i+1]<<std::endl; 
}

void Chain::velocityVerletIntegration(){
  computeForce();
  for (int i=0;i<nn;i++){
    v[2*i] += 0.5*force[2*i]/m;
    v[2*i+1] += 0.5*force[2*i+1]/m;
    x[2*i] += v[2*i];
    x[2*i+1] += v[2*i+1];
  }
  computeForce();
  for(int i=0;i<nn;i++){
    v[2*i] += 0.5*force[2*i]/m;
    v[2*i+1] += 0.5*force[2*i+1]/m;
  }
}
/*
void Chain::output(const std::string filename){
  using namespace std;*/ 
/*
void Chain::output(const std::string filename){
  using namespace std; 
  std::ofstream out(filename.c_str(),std::ios::out | std::ios::app);
    if (out.is_open()){
      out<<setw(10)<<"ks"<<setw(10)<<ks<<endl;
      out<<setw(10)<<"kb"<<setw(10)<<kb<<endl;
      out<<setw(10)<<"kp"<<setw(10)<<kp<<endl;
      out<<setw(10)<<"node"<<endl;
      out<<setw(10)<<nn<<endl;
      for (int i=0;i<nn;i++)
        out<<setw(10)<<x[2*i]<<setw(10)<<x[2*i+1]<<endl;
      
      out<<setw(10)<<"bond"<<endl;
      out<<setw(10)<<nb<<endl;
      for (int i=0;i<nb;i++)
        out<<setw(10)<<pBond->first[i]<<setw(10)<<pBond->next[i]<<setw(10)<<pBond->L0[i]<<endl;
      
      out<<setw(10)<<"angle"<<endl;
      out<<setw(10)<<na<<endl;
      for (int i=0;i<na;i++)
        out<<setw(10)<<pAngle->left[i]<<setw(10)<<pAngle->middle[i]<<setw(10)<<pAngle->right[i]<<setw(10)<<pAngle->ang0[i]<<endl;

      
    out.close();
    }else{
      std::cout<<"cannot open cell output file"<<std::endl;
    }
    
}
*/

void Chain::setCells(Cell* pCell_){
  pCell = pCell_;
}

int Chain::particleInsideCell(double xp, double yp, int idc){
  int i,j,c=0;
  int nForOne = pCell->nn/pCell->ns;
  double xci, yci, xcj,ycj;
  if (pCell->edgeFlag[idc]){
    //std::cout<<"cell "<<idc<<" edge "<<pCell->edgeFlag[idc]<<std::endl;
    return 0;
  }else{
    for (i=0,j=nForOne-1;i<nForOne; j=i++){
      xci=pCell->x[2*(idc*nForOne+i)];
      yci=pCell->x[2*(idc*nForOne+i)+1];
      xcj=pCell->x[2*(idc*nForOne+j)];
      ycj=pCell->x[2*(idc*nForOne+j)+1];
      if (((yci>yp)!=(ycj>yp)) && (xp<(xcj-xci)*(yp-yci)/(ycj-yci)+xci))
        c = !c;
    }
    return c;
  }
}

void Chain::moveTo(double x, double y){
/*  computeReference();
  double dx = x - xc[0];
  double dy = y - xc[1];
  for (int i=0;i<nn;i++){
    x[2*i] += dx;
    x[2*i+1] += dy;
  }*/
}

void Chain::moveOutside(double xp, double yp, int idc, int idp){
  double dx,dy,r;
  double dxMin,dyMin,rMin;//vector
  rMin = 1000000.;//just a big number
  dxMin=0.;
  dyMin=0.;
  int nForOne = pCell->nn/pCell->ns;
  for (int i=0;i<nForOne;i++){
    dx = pCell->x[2*(idc*nForOne+i)]-xp;
    dy = pCell->x[2*(idc*nForOne+i)+1]-yp;
    r = dx*dx + dy*dy;
    if (r < rMin){
      rMin = r;
      dxMin = dx;
      dyMin = dy;
    }
  }
  x[2*idp] += dxMin ;
  x[2*idp+1] += dyMin;
}

void Chain::thermalFluctuation(){
  double dx, dy;
  //double xorg,yorg;
  int inside=0;
  for (int i=0;i<nn;i++){
    dx = diffDist*random->gaussian();
    dy = diffDist*random->gaussian();
    //xorg = x[2*i];
    //yorg = x[2*i+1];
    x[2*i] += dx;
    x[2*i+1] += dy;
    for (int idc=0;idc<pCell->ns;idc++){
      inside = particleInsideCell(x[2*i],x[2*i+1],idc);
      if (inside){
        //std::cout<<"inside cell "<<idc<<" particle "<<i<<std::endl;
        x[2*i] -= 2*dx;
        x[2*i+1] -= 2*dy;
        //x[2*i] = xorg;
        //x[2*i+1] = yorg;
        break;
      }
    }
  }
  /* 
  for (int i=0;i<nn;i++){
    for (int idc=0;idc<pCell->ns;idc++){
      inside = particleInsideCell(x[2*i],x[2*i+1],idc);
      if (inside){
        //std::cout<<"cell "<<idc<<" particle  "<<i<<std::endl;
        moveOutside(x[2*i],x[2*i+1],idc,i);
        break;
      }
    }
  }*/
  /*double fx, fy, factor;
  factor = sqrt(kBT*m/damp);//dt=1
  for (int i=0;i<nn;i++){
    fx = factor*(random->uniform()-0.5);
    fy = factor*(random->uniform()-0.5);
    force[2*i] += fx;
    force[2*i+1] += fy;
  }*/
}

void Chain::initLJ(){
  if (lx==0 || ly==0)
    std::cout<<"particle doesn't know fluid geometry size"<<std::endl;
  cx = floor(lx/rCut);
  cy = floor(ly/rCut);
  std::cout<<"cx,cy "<<cx<<" "<<cy<<std::endl;
  lscl=new int[nn];
  head=new int[cx*cy];
  //std::cout<<"lscl,head size"<<nn<<" "<<cx*cy<<std::endl;
  rrCut=rCut*rCut;
  sig2=sigma*sigma;
  sig6=sig2*sig2*sig2;
}

void Chain::buildLinkList(){
  int cellSize=cx*cy;
  int idx,idy,idc;
  double dx,dy;
  //std::cout<<"lx,ly,cx,cy "<<lx<<" "<<ly<<" "<<cx<<" "<<cy<<std::endl;
  dx = double(lx/cx);
  dy = double(ly/cy);
  //if (dx<1 || dy < 1) std::cout<<"error, the cut off distance is too small"<<std::endl;
  //std::cout<<"lx,ly,cx,cy "<<lx<<" "<<ly<<" "<<cx<<" "<<cy<<std::endl;
  for (int i=0;i<cellSize;i++)
    head[i]=EMPTY;
  
  for (int i=0;i<nn;i++){
    idx = floor(x[2*i]/dx);
    idy = floor(x[2*i+1]/dy);
    idc = idy*cx + idx;
    //std::cout<<"dx,dy, i,idc "<<dx<<" "<<dy<<" "<<i<<" "<<idc<<std::endl;
    lscl[i]=head[idc];
    head[idc]=i;
    if (idc >cx*cy) std::cout<<"idc "<<idc<<" "<<cx*cy<<std::endl;
  }
}

void Chain::pairWiseInteraction(){
  int idc,idc_nb;
  int i,j;
  //double dx,dy,rsq;//r
  //double sig2,sig6,ri2,ri6,fcVal;

  //sig2=sigma*sigma;
  //sig6=sig2*sig2*sig2;
  // scan all the cells through x, and y direction
  for (int idcy=0;idcy<cy;idcy++){
    for (int idcx=0;idcx<cx;idcx++){
      idc = idcy*cx + idcx;
      if (head[idc]==EMPTY) continue;
      //scan neighbor cells
      for (int nby=idcy-1;nby<=idcy+1;nby++){
        for (int nbx=idcx-1;nbx<=idcx+1;nbx++){
          //when cell is on the edge, -1 could lead to out of boundary
          
          if (periodicX && periodicY){
            idc_nb = ((nby+cy)%cy)*cx+((nbx+cx)%cx);
          }else if (periodicX){
            if (nby== -1 || nby == cy) 
              continue;
            else 
              idc_nb = nby*cx + ((nbx+cx)%cx);
          }else if (periodicY){
            if (nbx == -1 || nbx == cx)
              continue;
            else
              idc_nb = ((nby+cy)%cy)*cx + nbx;
          }else{
            if (nby== -1 || nby == cy) 
              continue;
            if (nbx == -1 || nbx == cx)
              continue;
             
            idc_nb = nby*cx + nbx;
          }


          if (head[idc_nb]==EMPTY) continue;
          i=head[idc];
          while(i!=EMPTY){
            j=head[idc_nb];
            while(j!=EMPTY){
              if (i<j) LJForce(i,j);
              j=lscl[j];
            }//end of loop j
            i=lscl[i];
          }//eof loop i
        }//eof nbx
      }//eof nby
    }//eof idcx
  }//eof idcy
  
}

void Chain::LJForce(int i,int j){
  double dx,dy,rsq,ri2,ri6,fcVal;
  double halfX, halfY;
  halfX = 0.5*lx;
  halfY = 0.5*ly;
  dx = x[2*i]-x[2*j];
  dy = x[2*i+1]-x[2*j+1];
  if (periodicX){
    if (std::abs(dx)>halfX){
      if (dx > 0.)   
        dx=x[2*i]-x[2*j]-lx;
      else
        dx=x[2*i]-x[2*j]+lx;
    }
  }
  if (periodicY){
    if (std::abs(dy)>halfY){
      if (dy > 0.)   
        dy=x[2*i+1]-x[2*j+1]-ly;
      else
        dy=x[2*i+1]-x[2*j+1]+ly;
    }
  }
  rsq = dx*dx + dy*dy;
  //std::cout<<"i,j"<<i<<" "<<j<<" "<<rsq<<std::endl;
  // periodic conditions not considered yet 
  if (rsq < rrCut){
    ri2=1.0/rsq;ri6=ri2*ri2*ri2;
    //r=sqrt(rsq);
    fcVal=48.0*epsilon*ri2*ri6*sig6*(sig6*ri6-0.5);
    force[2*i] += fcVal*dx;
    force[2*i+1] += fcVal*dy;
    force[2*j] -= fcVal*dx;
    force[2*j+1] -= fcVal*dy;
  }
}

void Chain::writeGeometry(const std::string filename){
  using namespace std; 
  std::ofstream out(filename.c_str(),std::ios::out | std::ios::app);
    if (out.is_open()){
      for (int i=0;i<nn;i++)
        out<<setw(10)<<x[2*i]<<" "<<setw(10)<<x[2*i+1]<<endl; 
    }else{
      cout<<"cannot open output file"<<endl;
    }
}


void Chain::writeForce(const std::string filename){
  using namespace std; 
  ofstream out(filename.c_str(),ios::out | ios::app);
    if (out.is_open()){
      for (int i=0;i<nn;i++)
        out<<setw(10)<<force[2*i]<<" "<<setw(10)<<force[2*i+1]<<endl; 
      }else{
        cout<<"cannot open output file"<<endl;
      }
   
}

void Chain::writeVelocity(const std::string filename){
  using namespace std; 
  ofstream out(filename.c_str(),ios::out | ios::app);
    if (out.is_open()){
      for (int i=0;i<nn;i++)
        out<<setw(10)<<v[2*i]<<" "<<setw(10)<<v[2*i+1]<<endl; 
      }else{
        cout<<"cannot open output file"<<endl;
      }
}

void Chain::writeLog(const std::string filename){
  using namespace std; 
  ofstream out(filename.c_str(),ios::out | ios::app);
    if (out.is_open()){
      out<<"-----------Solid Chain-----------"<<endl;
      out<<"ks "<<ks<<" kb "<<kb<<" kp "<<kp<<endl;
      out<<"node: "<<nn<<" bond: "<<nb<<" angle: "<<na<<endl;
    }else{
      cout<<"cannot open log file"<<endl;
    }
}
