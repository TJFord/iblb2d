# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>
# include <cmath>

# include "cell.h"

Cell::Cell(){
  x=NULL;
  xtmp=NULL;
  v=NULL;
  x0=NULL;
  xR=NULL;
  force=NULL;
  pBond=NULL;
  pAngle=NULL;
  
  rho=0.;
  m=0.;
  ks=0.;
  kb=0.;
  kp=0.;
  g=9.8;
}

Cell& Cell::operator=(const Cell& rhs){
  if(this!=&rhs){
    delete [] x;
    delete [] xtmp;
    delete [] v;
    delete [] x0;
    delete [] xR;
    delete [] force;
    delete pBond;
    delete pAngle;
   
    nn=rhs.nn; nb=rhs.nb; na=rhs.na;
    rho=rhs.rho; m=rhs.m;
    ks=rhs.ks; kb=rhs.kb;
    kp=rhs.kp;

    x=new double[nn*2];
    xtmp=new double[nn*2];
    x0=new double[nn*2];
    xR=new double[nn*2];
    v=new double[nn*2];
    force=new double[nn*2];
    for (int i=0;i<nn;i++){
      x[2*i]=rhs.x[2*i];
      x[2*i+1]=rhs.x[2*i+1];
      xtmp[2*i]=rhs.xtmp[2*i];
      xtmp[2*i+1]=rhs.xtmp[2*i+1];
      x0[2*i]=rhs.x0[2*i];
      x0[2*i+1]=rhs.x0[2*i+1];
      xR[2*i]=rhs.xR[2*i];
      xR[2*i+1]=rhs.xR[2*i+1];
      v[2*i]=rhs.v[2*i];
      v[2*i+1]=rhs.v[2*i+1];
      force[2*i]=rhs.force[2*i];
      force[2*i+1]=rhs.force[2*i+1];
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

void Cell::readInput(const std::string filename){
  std::ifstream in(filename.c_str(),std::ios::in);
  if (in.is_open()){
    int tmp=0;
    std::string str;
    while(in>>str){
      //in>>str;
      if(str.compare("m")==0){
        in>>m;
      }else if(str.compare("ks")==0){
        in>>ks;
      }else if(str.compare("kb")==0){
        in>>kb;
      }else if(str.compare("kp")==0){
        in>>kp;
      }else if(str.compare("node")==0){
        in>>nn;
        x=new double[nn*2];
        xtmp=new double[nn*2];
        x0=new double[nn*2];
        xR=new double[nn*2];
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

Cell::~Cell(){
  delete [] x;
  delete [] xtmp;
  delete [] v;
  delete [] x0;
  delete [] xR;
  delete [] force;
  delete pBond;
  delete pAngle;
}

void Cell::init(){
  double dx,dy;
  xc0[0]=0.;
  xc0[1]=0.;
  for (int i=0;i<nn;i++){
    xR[2*i] =x[2*i];
    xR[2*i+1]=x[2*i+1];
    x0[2*i] =x[2*i];
    x0[2*i+1]=x[2*i+1];
    xc0[0] += x[2*i];
    xc0[1] += x[2*i+1];
    //std::cout<<"xhalf "<<xhalf[2*i]<<" "<<xhalf[2*i+1]<<std::endl;
  }
  xc0[0]/=nn;
  xc0[1]/=nn;
  dx = x[0]-xc0[0];
  dy = x[1]-xc0[1];
  radius = sqrt(dx*dx+dy*dy);
  angRef = acos(dx/radius);
}

void Cell::update(){
  for (int i=0;i<nn;i++){
    x[2*i] = xtmp[2*i]+v[2*i];
    x[2*i+1] =xtmp[2*i+1]+v[2*i+1];
    //x[2*i] = xtmp[2*i]+100*v[2*i];
    //x[2*i+1] =xtmp[2*i+1]+100*v[2*i+1];
  }
}

void Cell::updateHalf(){
  for (int i=0;i<nn;i++){
    xtmp[2*i] =x[2*i];
    xtmp[2*i+1]=x[2*i+1];
    x[2*i] += 0.5*v[2*i];
    x[2*i+1] += 0.5*v[2*i+1];
    //x[2*i] += 100*v[2*i];
    //x[2*i+1] += 100*v[2*i+1];
  }
}

void Cell::nondimension(const Units& unt){
  ks *= unt.dt*unt.dt/unt.dm;
  kb *= unt.dt*unt.dt/unt.dm/unt.dx/unt.dx;
  kp *= unt.dt*unt.dt/unt.dm/unt.dx;
  m /= unt.dm;
  g *= unt.dt*unt.dt/unt.dx;
}

void Cell::bondHarmonicForce(){
  double dx,dy,rsq,r,fv,dr;
  int i1,i2;
  for (int i=0;i<nb;i++){
    i1=pBond->first[i];
    i2=pBond->next[i];
    dx=x[2*i1]-x[2*i2];
    dy=x[2*i1+1]-x[2*i2+1];
    rsq = dx*dx+dy*dy;
  /*  if (i==nb-1){
      std::cout<<"x "<<i1<<" "<<xhalf[2*i1]<<" "<<xhalf[2*i1+1]<<std::endl;
      std::cout<<"x "<<xhalf[2*i2]<<" "<<xhalf[2*i2+1]<<std::endl;
      std::cout<<"rsq"<<rsq<<std::endl;
    }*/
    //std::cout<<"rsq"<<rsq<<" i1 "<<i1<<" "<<i2<<std::endl;
    r=sqrt(rsq);
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

void Cell::angleBendForce(){
  int i1, i2, i3;
  double dx1, dy1, dx2,dy2, rsq1, rsq2, r1,r2;
  double c, s;//cosine, sine
  double dtheta,a,a11,a12,a22,f1x,f1y,f3x,f3y;
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
    
    dtheta = acos(c)-pAngle->ang0[i];
    a = -kb*dtheta*s;
    a11 = a*c/rsq1;
    a12 = -a/(r1*r2);
    a22 = a*c/rsq2;
    f1x = a11*dx1 + a12*dx2;
    f1y = a11*dy1 + a12*dy2;
    f3x = a22*dx2 + a12*dx1;
    f3y = a22*dy2 + a12*dy1;
    
    force[2*i1] += f1x;
    force[2*i1+1] += f1y;
    force[2*i2] -= f1x+f3x;
    force[2*i2+1] -= f1y+f3y;
    force[2*i3] += f3x;
    force[2*i3+1] += f3y;

  }
}

void Cell::computeForce(){
  for (int i=0;i<nn;i++){
    force[2*i]=0.;
    force[2*i+1]=0.;
  }
  bondHarmonicForce();
  angleBendForce();
  //for (int i=0;i<nn;i++)
  //  force[2*i] += g;
}

void Cell::velocityVerletIntegration(){
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
void Cell::output(const std::string filename){
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

void Cell::computeReference(){
  double dx,dy;
  double angle,dtheta;
  double c,s;
  xc[0]=0.;
  xc[1]=0.;
  for(int i=0;i<nn;i++){
    xc[0] += x[2*i];
    xc[1] += x[2*i+1];
  }
  xc[0] /= nn;
  xc[1] /= nn;
  //std::cout<<"center "<<center[0]<<" "<<center[1]<<std::endl;
  //std::cout<<"xc "<<xc[0]<<" "<<xc[1]<<std::endl;
  dx = x[0]-xc[0];
  dy = x[1]-xc[1];
  if (dx > radius)
    dx = radius;
  else if (dx < -radius)
    dx = -radius;
  angle = acos(dx/radius);
  dtheta = angle - angRef;
 // std::cout<<"dx "<<dx<<" "<<dy<<" "<<angle<<std::endl;
  c = cos(dtheta);
  s = sin(dtheta);
  for (int i=0;i<nn;i++){
    dx = x0[2*i]-xc0[0];
    dy = x0[2*i+1]-xc0[1];
    xR[2*i] = xc[0]+ c*dx - s*dy;
    xR[2*i+1] = xc[1]+ s*dx + c*dy;
    //std::cout<<"X0 "<<x0[2*i]<<" "<<x0[2*i+1]<<" "<<c<<" "<<s<<std::endl;
  }
}

void Cell::computeRigidForce(){
 double dx,dy;
  for (int i=0;i<nn;i++){
    dx=x[2*i]-xR[2*i];
    dy=x[2*i+1]-xR[2*i+1];
    //rsq = dx*dx+dy*dy;
  /*  if (i==nb-1){
      std::cout<<"x "<<i1<<" "<<xhalf[2*i1]<<" "<<xhalf[2*i1+1]<<std::endl;
      std::cout<<"x "<<xhalf[2*i2]<<" "<<xhalf[2*i2+1]<<std::endl;
      std::cout<<"rsq"<<rsq<<std::endl;
    }*/
    //std::cout<<"rsq"<<rsq<<" i1 "<<i1<<" "<<i2<<std::endl;

    force[2*i] =-ks*dx;
    force[2*i+1] = -ks*dy;
    force[2*i] += 1e-2*g;
    //std::cout<<"restore force"<<force[2*i]<<" "<<force[2*i+1]<<" g"<<g<<std::endl;
  } 
}

void Cell::writeGeometry(const std::string filename){
  using namespace std; 
  std::ofstream out(filename.c_str(),std::ios::out | std::ios::app);
    if (out.is_open()){
      for (int i=0;i<nn;i++)
        out<<setw(10)<<x[2*i]<<" "<<setw(10)<<x[2*i+1]<<endl; 
    }else{
      cout<<"cannot open output file"<<endl;
    }
}


void Cell::writeReferenceGeometry(const std::string filename){
  using namespace std; 
  ofstream out(filename.c_str(),ios::out | ios::app);
    if (out.is_open()){
      for (int i=0;i<nn;i++)
        out<<setw(10)<<xR[2*i]<<" "<<setw(10)<<xR[2*i+1]<<endl; 
    }else{
      cout<<"cannot open output file"<<endl;
    } 
}

void Cell::writeForce(const std::string filename){
  using namespace std; 
  ofstream out(filename.c_str(),ios::out | ios::app);
      if (out.is_open()){
        for (int i=0;i<nn;i++)
          out<<setw(10)<<force[2*i]<<" "<<setw(10)<<force[2*i+1]<<endl; 
      }else{
        cout<<"cannot open output file"<<endl;
      }
   
}
