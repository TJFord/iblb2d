# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>
# include <cmath>

# include "rigid.h"

Rigid::Rigid(){
  x=NULL;
  xc=NULL;
  xc0=NULL;
  xtmp=NULL;
  v=NULL;
  x0=NULL;
  xR=NULL;
  force=NULL;
 
  radius = NULL;
  angRef = NULL;

  rho=0.;
  m=0.;
  ks=0.;
  kb=0.;
  kp=0.;
  g=9.8;
}

Rigid& Rigid::operator=(const Rigid& rhs){
  if(this!=&rhs){
    delete [] x;
    delete [] xc;
    delete [] xc0;
    delete [] xtmp;
    delete [] v;
    delete [] x0;
    delete [] xR;
    delete [] force;
    
    delete [] angRef;
    delete [] radius;

    nn=rhs.nn; nb=rhs.nb; na=rhs.na;
    ns=rhs.ns;
    rho=rhs.rho; m=rhs.m;
    ks=rhs.ks; kb=rhs.kb;
    kp=rhs.kp;
    g=rhs.g;

    x=new double[nn*2];
    xtmp=new double[nn*2];
    x0=new double[nn*2];
    xR=new double[nn*2];
    v=new double[nn*2];
    force=new double[nn*2];
    angRef = new double[nn];
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
      angRef[i]=rhs.angRef[i];
    }
    xc=new double[ns*2];
    xc0=new double[ns*2];
    radius=new double[ns];
    for (int i=0;i<ns;i++){
      xc[2*i]=rhs.xc[2*i];
      xc[2*i+1]=rhs.xc[2*i+1];
      xc0[2*i]=rhs.xc0[2*i];
      xc0[2*i+1]=rhs.xc0[2*i+1];
      radius[i]=rhs.radius[i];
    }
  }
  return *this;
}

void Rigid::readInput(const std::string filename){
  std::ifstream in(filename.c_str(),std::ios::in);
  if (in.is_open()){
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
      }else if(str.compare("ns")==0){
        in>>ns;
        xc=new double[ns*2];
        xc0=new double[ns*2];
        radius=new double[ns];
      }else if(str.compare("node")==0){
        in>>nn;
        x=new double[nn*2];
        xtmp=new double[nn*2];
        x0=new double[nn*2];
        xR=new double[nn*2];
        v=new double[nn*2];
        force=new double[nn*2];
        angRef=new double[nn];

        for (int i=0;i<nn;i++)
          in>>x[2*i]>>x[2*i+1];
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

Rigid::~Rigid(){
  delete [] x;
  delete [] xc;
  delete [] xc0;
  delete [] xtmp;
  delete [] v;
  delete [] x0;
  delete [] xR;
  delete [] force;
  delete [] angRef;
  delete [] radius;
}

void Rigid::init(){
  double dx,dy;
  for (int i=0;i<nn;i++){
    xR[2*i] =x[2*i];
    xR[2*i+1]=x[2*i+1];
    x0[2*i] =x[2*i];
    x0[2*i+1]=x[2*i+1];
    //xc0[0] += x[2*i];
    //xc0[1] += x[2*i+1];
    v[2*i]=0.;
    v[2*i+1]=0.;
    //std::cout<<"xhalf "<<xhalf[2*i]<<" "<<xhalf[2*i+1]<<std::endl;
  }
  nForOne = nn/ns;
  int tmp;
  for (int i=0;i<ns;i++){
    xc0[2*i]=0.;
    xc0[2*i+1]=0.;
    for (int j=0;j<nForOne;j++){
      tmp = i*nForOne + j;
      xc0[2*i] += x[2*tmp];
      xc0[2*i+1] += x[2*tmp+1];
    }
    xc0[2*i] /= nForOne;
    xc0[2*i+1] /= nForOne;
    tmp = i*nForOne;
    dx = x[2*tmp]-xc0[2*i];
    dy = x[2*tmp+1]-xc0[2*i+1];
    radius[i] = sqrt(dx*dx+dy*dy);// same radius rigid
    angRef[i*nForOne]=acos(dx/radius[i]);
    for (int j=1;j<nForOne;j++){
      tmp = i*nForOne+j;
      dx = x[2*tmp]-xc0[2*i];
      dy = x[2*tmp+1]-xc0[2*i+1];
      angRef[tmp]=acos(dx/radius[i]);
    }
  }

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
}

void Rigid::update(){
  for (int i=0;i<nn;i++){
    x[2*i] = xtmp[2*i]+v[2*i];
    x[2*i+1] =xtmp[2*i+1]+v[2*i+1];
    //x[2*i] = xtmp[2*i]+100*v[2*i];
    //x[2*i+1] =xtmp[2*i+1]+100*v[2*i+1];
  }
}

void Rigid::updateHalf(){
  for (int i=0;i<nn;i++){
    xtmp[2*i] =x[2*i];
    xtmp[2*i+1]=x[2*i+1];
    x[2*i] += 0.5*v[2*i];
    x[2*i+1] += 0.5*v[2*i+1];
  }
}

void Rigid::nondimension(const Units& unt){
  // ks: energy Es = 0.5*ks*dr^2
  ks *= unt.dt*unt.dt/unt.dm;
  // kb: energy units, Eb = 0.5*kb*dtheta^2
  kb *= unt.dt*unt.dt/unt.dm/unt.dx/unt.dx;
  // kp: energy Ep = 0.5*kp*dA^2,for 2d, this 
  // is not derived rigorously based on derivative
  // thus, the force is simplified as F= kp*dA;
  kp *= unt.dt*unt.dt*unt.dx/unt.dm;
  m /= unt.dm;
  g *= unt.dt*unt.dt/unt.dx;
  std::cout<<"ks, kb, kp "<<ks<<" "<<kb<<" "<<kp<<std::endl;
}

void Rigid::velocityVerletIntegration(){
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
void Rigid::output(const std::string filename){
  using namespace std;*/ 
/*
void Rigid::output(const std::string filename){
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

void Rigid::computeReference(){
  double dx,dy;
  double angle,dtheta;
  double c,s;
  int half;
  int tmp;
  if (nForOne%2)
    half = (nForOne+1)/2;
  else
    half = nForOne/2;
  for (int j=0;j<ns;j++){
    xc[2*j]=0.;
    xc[2*j+1]=0.;
    for(int i=0;i<nForOne;i++){
      tmp = j*nForOne + i;
      xc[2*i] += x[2*tmp];
      xc[2*i+1] += x[2*tmp+1];
    }
    xc[2*j] /= nForOne;
    xc[2*j+1] /= nForOne;

    dtheta=0.;
    tmp = j*nForOne;
    dx = x[2*tmp]-xc[2*j];
    dy = x[2*tmp+1]-xc[2*j+1];
    if (dx > radius[j])
      dx = radius[j];
    else if (dx < -radius[j])
      dx = -radius[j];
    angle = acos(dx/radius[j]);

    dtheta += angle - angRef[tmp];
    for (int i=1;i<half;i++){
      tmp = j*nForOne + i;
      dx=x[2*tmp]-xc[j];
      if (dx > radius[j])
        dx = radius[j];
      else if (dx < -radius[j])
        dx = -radius[j];
      angle = acos(dx/radius[j]);
      dtheta += angle - angRef[tmp];
    }
    dtheta /= half;
     // std::cout<<"dx "<<dx<<" "<<dy<<" "<<angle<<std::endl;
    c = cos(dtheta);
    s = sin(dtheta);
    for (int i=0;i<nForOne;i++){
      tmp = j*nForOne + i;
      dx = x0[2*tmp]-xc0[2*j];
      dy = x0[2*tmp+1]-xc0[2*j+1];
      xR[2*tmp] = xc[2*j]+ c*dx - s*dy;
      xR[2*tmp+1] = xc[2*j+1]+ s*dx + c*dy;
        //std::cout<<"X0 "<<x0[2*i]<<" "<<x0[2*i+1]<<" "<<c<<" "<<s<<std::endl;
    }
  }
}

void Rigid::computeForce(){
  double dx,dy;
  /*for (int i=0;i<nn;i++){
    force[2*i]=0.;
    force[2*i+1]=0.;
  }*/
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
    force[2*i] += 1e-3*g;//1e-4*g;
    //std::cout<<"restore force"<<force[2*i]<<" "<<force[2*i+1]<<" g"<<g<<std::endl;
  } 
}

void Rigid::moveTo(double x, double y){
/*  computeReference();
  double dx = x - xc[0];
  double dy = y - xc[1];
  for (int i=0;i<nn;i++){
    x[2*i] += dx;
    x[2*i+1] += dy;
  }*/
}

void Rigid::writeGeometry(const std::string filename){
  using namespace std; 
  std::ofstream out(filename.c_str(),std::ios::out | std::ios::app);
    if (out.is_open()){
      for (int i=0;i<nn;i++)
        out<<setw(10)<<x[2*i]<<" "<<setw(10)<<x[2*i+1]<<endl; 
    }else{
      cout<<"cannot open output file"<<endl;
    }
}


void Rigid::writeReferenceGeometry(const std::string filename){
  using namespace std; 
  ofstream out(filename.c_str(),ios::out | ios::app);
    if (out.is_open()){
      for (int i=0;i<nn;i++)
        out<<setw(10)<<xR[2*i]<<" "<<setw(10)<<xR[2*i+1]<<endl; 
    }else{
      cout<<"cannot open output file"<<endl;
    } 
}

void Rigid::writeForce(const std::string filename){
  using namespace std; 
  ofstream out(filename.c_str(),ios::out | ios::app);
    if (out.is_open()){
      for (int i=0;i<nn;i++)
        out<<setw(10)<<force[2*i]<<" "<<setw(10)<<force[2*i+1]<<endl; 
      }else{
        cout<<"cannot open output file"<<endl;
      }
   
}

void Rigid::writeVelocity(const std::string filename){
  using namespace std; 
  ofstream out(filename.c_str(),ios::out | ios::app);
    if (out.is_open()){
      for (int i=0;i<nn;i++)
        out<<setw(10)<<v[2*i]<<" "<<setw(10)<<v[2*i+1]<<endl; 
      }else{
        cout<<"cannot open output file"<<endl;
      }
}

void Rigid::writeLog(const std::string filename){
  using namespace std; 
  ofstream out(filename.c_str(),ios::out | ios::app);
    if (out.is_open()){
      out<<"-----------Solid Rigid-----------"<<endl;
      out<<"ks "<<ks<<" kb "<<kb<<" kp "<<kp<<endl;
      out<<"node: "<<nn<<" bond: "<<nb<<" angle: "<<na<<endl;
    }else{
      cout<<"cannot open log file"<<endl;
    }
}
