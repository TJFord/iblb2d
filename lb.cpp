# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>
# include "lb.h"
# include "units.h"
# include "boundary.h"

//using namespace std;

const int LB::d = 2;
const int LB::q = 9;
const double LB::cs2 = 1.0/3.0;

const int LB::c[9][2] = {
    {0,0},
    {1,0}, {0,1}, {-1,0}, {0,-1},
    {1,1}, {-1,1}, {-1,-1}, {1,-1}
};

const double LB::w[9] = { 4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

void swap(double &f1, double &f2){
  double tmp = f1;
  f1 = f2;
  f2 = tmp;
}

LB::LB(){
  lx=0;ly=0;nf=0;nt=0;
  xy2idx=NULL;
  nbList=NULL;
  f=NULL;
  ft=NULL;
  coor=NULL;
  v=NULL;
  pUnits=NULL;
  omega=1.0;
  nbc=0;
  pBC=NULL;
  collisionFun=NULL;
  nv=0;no=0;nb=0;
  pVel=NULL;
  pOpen=NULL;
  pBB=NULL;
  force = NULL;
}

LB& LB::operator=(const LB& rhs){
  if (this!=&rhs){
    delete [] xy2idx;
    delete [] nbList;
    delete [] f;
    delete [] ft;
    delete [] coor;
    delete [] v;
    delete [] pBC;
    delete [] pVel;
    delete [] pOpen;
    delete [] pBB;
    delete [] force;
    delete pUnits;

    lx=rhs.lx; ly=rhs.ly; nf=rhs.nf; nt=rhs.nt;
    omega = rhs.omega;
    nbc=rhs.nbc; nv=rhs.nv; no=rhs.no; nb=rhs.nb;
    collisionScheme = rhs.collisionScheme;
    collisionFun = rhs.collisionFun;
   
    int tmp;
    xy2idx = new int[lx*ly];
    tmp = lx*ly;
    for (int i=0;i<tmp;i++)
      xy2idx[i]=rhs.xy2idx[i];
    
    nbList = new int[nf*(q-1)];
    tmp = nf*(q-1);
    for (int i=0;i<tmp;i++)
      nbList[i]=rhs.nbList[i];

    f = new double[nt*q];
    tmp = nt*q;
    for (int i=0;i<tmp;i++)
      f[i]=rhs.f[i];

    ft = new double[nt*q];
    tmp = nt*q;
    for (int i=0;i<tmp;i++)
      ft[i]=rhs.ft[i];
    
    coor = new int[nf*d];
    tmp = nf*d;
    for (int i=0;i<tmp;i++)
     coor[i]=rhs.coor[i];
    
    v = new double[nf*d];
    tmp = nf*d;
    for (int i=0;i<tmp;i++)
     v[i]=rhs.v[i];
    
    force = new double[nf*d];
    tmp = nf*d;
    for (int i=0;i<tmp;i++)
     force[i]=rhs.force[i];
  
    /*
    for (int i=0;i<nbc;i++){
      pBC[i].dynf=rhs.pBC[i].dynf;
      *(pBC[i].selfData) = *(rhs.pBC[i].selfData);
    }*/
    //-----velocity boundary-----------//
    if(nv){ 
      pVel = new velData[nv];
      for(int i=0;i<nv;i++){
        pVel[i]=velData(rhs.pVel[i].size);
        //pVel[i].size=rhs.pVel[i].size;
        pVel[i].norm[0]=rhs.pVel[i].norm[0];
        pVel[i].norm[1]=rhs.pVel[i].norm[1];
        for(int j=0;j<pVel[i].size;j++){
          pVel[i].node[j]=rhs.pVel[i].node[j];
          pVel[i].ux[j]=rhs.pVel[i].ux[j];
          pVel[i].uy[j]=rhs.pVel[i].uy[j];
        }
      }
     }
    //-----open boundary-----------//
    if(no){
      pOpen = new openData[no];
      for(int i=0;i<no;i++){
        pOpen[i]=openData(rhs.pOpen[i].size);
        //pOpen[i].size=rhs.pOpen[i].size;
        pOpen[i].norm[0]=rhs.pOpen[i].norm[0];
        pOpen[i].norm[1]=rhs.pOpen[i].norm[1];
        for(int j=0;j<pOpen[i].size;j++){
          pOpen[i].node[j]=rhs.pOpen[i].node[j];
          pOpen[i].nodeNext[j]=rhs.pOpen[i].nodeNext[j];
        }
      }
    }
    //-----bounceback boundary-----------//
    if(nb){
      pBB = new bbData[nb];
      for(int i=0;i<nb;i++){
        pBB[i]=bbData(rhs.pBB[i].size);
        //pBB[i].size=rhs.pBB[i].size;
        for(int j=0;j<pBB[i].size;j++){
          pBB[i].node[j]=rhs.pBB[i].node[j];
        }
      }
    }
    pBC = new Dynamics[nbc];
    tmp=0;
    if(nv){
      for (int i=0;i<nv;i++){
        pBC[tmp].selfData = (void*)&pVel[i];
        pBC[tmp].dynf = rhs.pBC[tmp].dynf;
        tmp++;
      } 
    }
    if(no){
      for (int i=0;i<no;i++){
        pBC[tmp].selfData = (void*)&pOpen[i];
        pBC[tmp].dynf = rhs.pBC[tmp].dynf;
        tmp++;
      } 
    }
    if(nb){
      for (int i=0;i<nb;i++){
        pBC[tmp].selfData = (void*)&pBB[i];
        pBC[tmp].dynf = rhs.pBC[tmp].dynf;
        tmp++;
      } 
    }
    pUnits = new Units;
    *pUnits = *(rhs.pUnits);
  } 
  return *this;
}

void LB::readInput(const std::string filename)
{
  using namespace std;
  ifstream in(filename.c_str(),ios::in);
  if (in.is_open())
    {
     // cout<<"I get here"<<endl;
      string str; 
      int tmp,size(0);
      pUnits = new Units;
      while (in>>str){
        //in>>str;
        if(str.compare("scheme")==0){
          in >>collisionScheme;
        }else if(str.compare("dx")==0){
          in >> pUnits->dx;
        } else if (str.compare("dt")==0){
          in >> pUnits->dt;
        }else if(str.compare("tau")==0){
          in >> pUnits->tau;
        }else if(str.compare("rho")==0){
          in >> pUnits->rho;
        }else if(str.compare("vis")==0){
          in >> pUnits->vis;
        }else if(str.compare("u")==0){
          in >> pUnits->u;
        }else if(str.compare("lx") == 0)
        {in >> lx;}
        else if (str.compare("ly") == 0)
        { in >> ly;}
        else if(str.compare("nf")==0){
        in>>nf;
        }else if(str.compare("nt")==0){
          in>>nt;f=new double[nt*q];
          ft=new double[nt*q];
        }
       // else if (str.compare("nf") == 0)
       // { in >> nf; f = new double[q*(nf+1)];}
        else if (str.compare("nbList") == 0){
          nbList = new int[nf*(q-1)];
          coor = new int[nf*d];
          force = new double[nf*d];
          v = new double[nf*d];
          int nadj = q - 1;
          for (int i = 0; i < nf; i++){
            in >> tmp;//current node number
            for (int j = 0; j < nadj; j++){
              in >> nbList[i*nadj + j];
            }
            in >> coor[i*d] >> coor[i*d + 1]; //x,y
          }
        }else if(str.compare("xy2idx") == 0){
          xy2idx = new int[lx*ly];
            for (int row = 0; row < ly; row++)
              for (int col = 0; col < lx; col++)
                in >> xy2idx[row*lx + col];
        }else if(str.compare("boundaries")==0){
          in >> tmp;
          while(  in >> str){
            if(str.compare("velocity")==0){
              in >> nv;
              pVel = new velData[nv];
              for (int i = 0; i < nv; i++){
                in >> size;
                pVel[i]=velData(size);
                nbc++;
                in >>pVel[i].norm[0] >>pVel[i].norm[1];
                for (int j = 0; j < pVel[i].size; j++){
                  in >> pVel[i].node[j] >> pVel[i].ux[j] >> pVel[i].uy[j];
                }
              }
            }else if(str.compare("open")==0){
              in >> no;
              pOpen = new openData[no];
              for (int i = 0; i < no; i++){
                in >> size;
                pOpen[i]=openData(size);
                nbc++;
                in >>pOpen[i].norm[0] >>pOpen[i].norm[1];
                for (int j = 0; j < pOpen[i].size; j++){
                  in >> pOpen[i].node[j]>>pOpen[i].nodeNext[j];
                }
              }
            }else if(str.compare("bounceback")==0){
              nb = 1;
              pBB = new bbData[1];
              in>>size; pBB[0]=bbData(size);
              nbc++;
              for (int j = 0; j < pBB[0].size; j++){
                in >> pBB[0].node[j];
              }
            }
          }
        }else{
          cout<<"wrong keywords in input! "<<str<<endl;
          return;
          //break;
        }
      }

      cout<<"succesfully read fluid file"<<endl;
      //pUnits->calculateLBPara();

      pBC = new Dynamics[nbc];
      tmp = 0;
      if(nv){
        for(int i=0;i<nv;i++){
          pBC[tmp].selfData = (void*)&pVel[i];
          if (pVel[i].norm[0]==0){
            if (pVel[i].norm[1]==1){
              pBC[tmp].dynf = &upperZouHe;
              tmp++;cout<<"upperZouHe"<<endl;
            }else if(pVel[i].norm[1]==-1){
              pBC[tmp].dynf = &lowerZouHe;
              tmp++;cout<<"lowerZouHe"<<endl;
            }else
              cout<<"wrong norm"<<endl;
          }else if(pVel[i].norm[0]==1){
            pBC[tmp].dynf = &rightZouHe;
            tmp++;cout<<"rightZouHe"<<endl;
          }else if(pVel[i].norm[0]==-1){
            pBC[tmp].dynf = &leftZouHe;
            tmp++;cout<<"leftZouHe"<<endl;
          }
        }
      }
      if(no){
        for (int i=0;i<no;i++){
          pBC[tmp].selfData = (void*)&pOpen[i];
          if (pOpen[i].norm[0]==0){
            if (pOpen[i].norm[1]==1){
              pBC[tmp].dynf = &upperOpen;
              tmp++;
            }else if(pOpen[i].norm[1]==-1){
              pBC[tmp].dynf = &lowerOpen;
              tmp++;
            }else
              cout<<"wrong norm"<<endl;
          }else if(pOpen[i].norm[0]==1){
            pBC[tmp].dynf = &rightOpen;
            tmp++;cout<<"rightOpen"<<endl;
          }else if(pOpen[i].norm[0]==-1){
            pBC[tmp].dynf = &leftOpen;
            tmp++;cout<<"leftOpen"<<endl;
          }
        }
      }
      if(nb){
        pBC[tmp].selfData = (void*)&pBB[0];
        pBC[tmp].dynf = &bounceBack;
      }
      //in.close();
     } else {
      cout<<"cannot open input file"<<endl;
    }
}

void LB::reReadVelocity(const std::string filename)
{
  using namespace std;
  ifstream in(filename.c_str(),ios::in);
  int st;
  double uSqr;// = ux*ux + uy*uy;
  if (in.is_open()){
    for (int i=0;i<nf;i++){
      in>>v[2*i]>>v[2*i+1];//read velocity
      uSqr = v[2*i]*v[2*i]+v[2*i+1]*v[2*i+1];
      st = i*q;
      for (int j=0; j<q; j++){
        f[st+j] = computeEquilibrium(j, 1.0, v[2*i], v[2*i+1], uSqr); //initialize with zero velocity 
        ft[st+j]=f[st+j];
      }
    }
  }else{
    cout<<"reread velocity file file error"<<endl;
  }
}
 
LB::~LB(){
  delete [] xy2idx;
  delete [] nbList;
  delete [] f;
  delete [] ft;
  delete [] coor;
  delete [] v;
  delete [] pVel;
  delete [] pOpen;
  delete [] pBB;
  delete [] pBC;
  delete [] force;
  delete pUnits;

  std::cout<<"LB finished successfully"<<std::endl;
}

void LB::printInfor(){
  using namespace std;
  cout<<"********************************"<<endl;
  cout<<"IBLB program"<<endl;
  cout<<"--------------Fluid-------------"<<endl;
  cout<<"length L = "<<lx*pUnits->dx<<" m, width W = "<<ly*pUnits->dx<<" m"<<endl;
  cout<<"density = "<<pUnits->rho<<" kg/m^3"<<endl;
  cout<<"Kinetic viscosity  = "<<pUnits->vis<<" m^2/s"<<endl;
  cout<<"Characteristic L0 = "<<pUnits->dx<<endl;
  cout<<"Re = "<<pUnits->Re<<endl;
  cout<<"---------Lattice Boltzman-------"<<endl;
  cout<<"scheme = "<<collisionScheme<<endl;
  cout<<"tau = "<<pUnits->tau<<endl;
  cout<<"viscosity nu  = "<<pUnits->nu<<endl;
  cout<<"Applied LTX speed u_lb = "<<pUnits->u_lb<<endl;
  cout<<"dx = "<<pUnits->dx<<" m"<<endl;
  cout<<"dt = "<<pUnits->dt<<" s"<<endl;
  cout<<"********************************"<<endl;
}
 
void LB::writeLog(const std::string filename){
  using namespace std;
  ofstream out(filename.c_str(),ios::out | ios::app);
    if (out.is_open()){
      out<<"********************************"<<endl;
      out<<"IBLB program"<<endl;
      out<<"--------------Fluid-------------"<<endl;
      out<<"length L = "<<lx*pUnits->dx<<" m, width W = "<<ly*pUnits->dx<<" m"<<endl;
      out<<"density = "<<pUnits->rho<<" kg/m^3"<<endl;
      out<<"Kinetic viscosity  = "<<pUnits->vis<<" m^2/s"<<endl;
      out<<"Characteristic L0 = "<<pUnits->dx<<endl;
      out<<"Re = "<<pUnits->Re<<endl;
      out<<"---------Lattice Boltzman-------"<<endl;
      out<<"scheme = "<<collisionScheme<<endl;
      out<<"tau = "<<pUnits->tau<<endl;
      out<<"viscosity nu  = "<<pUnits->nu<<endl;
      out<<"Applied LTX speed u_lb = "<<pUnits->u_lb<<endl;
      out<<"dx = "<<pUnits->dx<<" m"<<endl;
      out<<"dt = "<<pUnits->dt<<" s"<<endl;
      out<<"********************************"<<endl;
    }else{
      cout<<"cannot open log file"<<endl;
    }
}

void LB::init(){
  pUnits->calculateLBPara();
 // if (collisionScheme.compare("stokes")==0) 
 //   pUnits->setStokesLBPara(0.0005);//set u_lb = 0.001
  omega = pUnits->omega;
  // initialize tensor Qxx, Qxy, Qyy;
  for (int i=0;i<q;i++){
    Qxx[i]=c[i][0]*c[i][0]-cs2;
    Qxy[i]=c[i][0]*c[i][1];
    Qyy[i]=c[i][1]*c[i][1]-cs2;
  }
  int st = 0;
  for (int i = 0; i < nf; i++){
    for (int j=0;j<d;j++){
      force[i*d+j]=0.0;
      v[i*d+j]=0.0;
    }
    //force[i*d] = 1.e-6;
    //force[i*d+1] =0.0;
    //cout<<"force"<<"i="<<i<<" "<<force[i*d]<<" "<<force[i*d+1]<<endl;
    st = i*q;
    for (int j=0; j<q; j++){
      f[st+j] = computeEquilibrium(j, 1.0, 0.0, 0.0, 0.0); //initialize with zero velocity 
      ft[st+j]=f[st+j];
    }
  }
  //force[d*(50*lx+50)]=1e-4;
  //force[d*(20*lx+20)+1]=1e-4;
 

  for(int i=nf;i<nt;i++){
    st = i*q;
    for(int j=0;j<q;j++){
      f[st+j]=0.0;
      ft[st+j]=0.0;
    }
  }
 
  if (collisionScheme.compare("bgk")==0)
    collisionFun=&LB::bgk;//fluid scheme
  else if (collisionScheme.compare("regularized")==0)
    collisionFun=&LB::regularized;//fluid scheme
  else if (collisionScheme.compare("stokes")==0)
    collisionFun=&LB::stokes;//fluid scheme
  
//  using namespace std;
 //--------rescale the velocity to LB units--------- //
  if( nv ){
    double velscale = pUnits->dx/pUnits->dt;
    for(int i=0;i<nv;i++){
      for (int j=0;j<pVel[i].size;j++){
        pVel[i].ux[j] /= velscale;
        pVel[i].uy[j] /= velscale;
      }
    }
  }
}

void LB::computeMacros(int id, double *rho, double * ux, double *uy){
  int st = id*q; //save computation time
  double upperLine = f[st+2] + f[st+5] + f[st+6];
  double mediumLine = f[st] + f[st+1] + f[st+3];
  double lowerLine = f[st+4] + f[st+7] + f[st+8];
  *rho = upperLine + mediumLine + lowerLine;
  *ux = (f[st+1]+f[st+5]+f[st+8] -(f[st+3] +f[st+6] + f[st + 7]))/(*rho);
  *uy = (upperLine - lowerLine)/(*rho);
}

void LB::computeVelocity(){
  double rho;
  for (int i=0;i<nf;i++){
    computeMacros(i,&rho,&v[2*i],&v[2*i+1]);
  }
}

double LB::computeEquilibrium(int idx, double rho,double ux, double uy,double uSqr){
  double c_u = c[idx][0]*ux + c[idx][1]*uy;
  return rho*w[idx]*(1. + 3.*c_u + 4.5*c_u*c_u - 1.5*uSqr);
}

double LB::computeEqStokes(int idx, double rho,double ux, double uy){
  double c_u = c[idx][0]*ux + c[idx][1]*uy;
  return rho*w[idx]*(1. + 3.*c_u );
}

void LB::bgk(int id, void* selfData){
  //cout<<"bgk"<<endl;
  //double omega = *((double *)selfData);
  double rho, ux, uy;
  computeMacros(id,&rho, &ux, &uy);
  double uSqr = ux*ux + uy*uy;
  int st = id*q;
  for (int i=0;i<q;i++){
    f[st+i] *= (1-omega);
    f[st+i] += omega*computeEquilibrium(i,rho, ux, uy, uSqr);
  }
  v[2*id]=ux;
  v[2*id+1]=uy;
}

void LB::regularized(int id, void* selfData){
  double omega = *((double *)selfData);
  double rho, ux, uy;
  //double feq[9], fNeq[9],Qxx[9],Qxy[9],Qyy[9];
  //double neqPixx(0.0), neqPixy(0.0), neqPiyy(0.0);
  neqPixx =0.0;
  neqPixy =0.0;
  neqPiyy =0.0;
  computeMacros(id,&rho, &ux, &uy);
  double uSqr = ux*ux + uy*uy;
  int st = id*9;
  for (int i=0;i<9;i++){
    feq[i]=computeEquilibrium(i,rho, ux, uy, uSqr);
    fNeq[i]=f[st+i]-feq[i];
    neqPixx += Qxx[i]*fNeq[i];
    neqPixy += Qxy[i]*fNeq[i];
    neqPiyy += Qyy[i]*fNeq[i];
  }
  for (int i=0;i<9;i++){
    f[st+i] = feq[i]+(1-omega)*4.5*w[i]*(Qxx[i]*neqPixx + 2*Qxy[i]*neqPixy + Qyy[i]*neqPiyy);
  }
  v[2*id]=ux;
  v[2*id+1]=uy;
}

void LB::stokes(int id, void* selfData){
  double rho, ux, uy;
  computeMacros(id,&rho, &ux, &uy);
  int st = id*q;
  for (int i=0;i<q;i++){
    f[st+i] *= (1-omega);
    f[st+i] += omega*computeEqStokes(i,rho, ux, uy);
  }
  v[2*id]=ux;
  v[2*id+1]=uy;
}
//void LB::applyForce(int id, double fx, double fy){
void LB::applyForce(){
  double F[9];
  double ux, uy;
  double fx,fy;
  int st=0;
  int id;
  for (id=0;id<nf;id++){
    st = id*9;
    fx = force[id*d];
    fy = force[id*d+1];
    //computeMacros(id,&rho, &ux, &uy);
    ux=v[2*id];
    uy=v[2*id+1];
    F[0]=(1.-0.5*omega)*w[0]*(3.*((   -ux)*fx +(   -uy)*fy));
    F[1]=(1.-0.5*omega)*w[1]*(3.*(( 1.-ux)*fx +(   -uy)*fy)+9.*(ux*fx));
    F[2]=(1.-0.5*omega)*w[2]*(3.*((   -ux)*fx +( 1.-uy)*fy)+9.*(uy*fy));
    F[3]=(1.-0.5*omega)*w[3]*(3.*((-1.-ux)*fx +(   -uy)*fy)+9.*(ux*fx));
    F[4]=(1.-0.5*omega)*w[4]*(3.*((   -ux)*fx +(-1.-uy)*fy)+9.*(uy*fy));
    F[5]=(1.-0.5*omega)*w[5]*(3.*(( 1.-ux)*fx +( 1.-uy)*fy)+9.*((ux+uy)*(fx+fy)));
    F[6]=(1.-0.5*omega)*w[6]*(3.*((-1.-ux)*fx +( 1.-uy)*fy)+9.*((ux-uy)*(fx-fy)));
    F[7]=(1.-0.5*omega)*w[7]*(3.*((-1.-ux)*fx +(-1.-uy)*fy)+9.*((ux+uy)*(fx+fy)));
    F[8]=(1.-0.5*omega)*w[8]*(3.*(( 1.-ux)*fx +(-1.-uy)*fy)+9.*((ux-uy)*(fx-fy)));
    for (int i=0;i<9;i++)
      f[st+i] += F[i];
  }
}
/*swap algorithm, takes time*/
void LB::stream(){
  int tgt=0; 
  int st=0;
  //int ix,iy,tix,tiy,tmp;
  for (int i=0;i<nf;i++){
    st = i*q;
    for (int j=1;j<q;j++){
      ft[st+j]=f[st+j];
     // tgt = nbList[i*(q-1)+j-1];
     // ft[tgt*q+j]=f[tgt*q+j];
     // f[idx*q+j]=ftmp[st+j];
    }
  }
  for(int i=0;i<nf;i++){
    st = i*q;
    for (int j=1;j<q;j++){
      tgt = nbList[i*(q-1)+j-1];
      f[tgt*q+j]=ft[st+j];
      /*ix = coor[2*i];
      iy = coor[2*i+1];
      tix = ix + c[j][0];
      tiy = iy + c[j][1];
      tmp = tiy*lx+tix;

      tgt = xy2idx[tiy*lx+tix];
      //std::cout<<"tgt "<<tgt<<" lx*ly "<<lx*ly<<" ix "<<ix<<" "<<tix<<" iy "<<iy<<" "<<tiy<<std::endl;
      f[tgt*q+j]=ft[st+j];*/
    }   
  }
}

void LB::collide(){
  for (int i=0;i<nf;i++){
    (this->*collisionFun)(i,(void *)&omega);
  }
  /*dynType pf;
  for (int i=0;i<nf;i++)  
    {
      pf = pDyn[i].dyn;
      (this->*pf)(i,pDyn[i].selfData);
    }
	*/
}

void LB::collideSwap(){
  int st=0;
  for (int i=0;i<nf;i++){
    (this->*collisionFun)(i,(void *)&omega);
    st = i*q;
    swap(f[st+1],f[st+3]);
    swap(f[st+2],f[st+4]);
    swap(f[st+5],f[st+7]);
    swap(f[st+6],f[st+8]);
  }
}

void LB::streamSwap(){
  int st=0;
  int tgt=0;
  for (int i=0;i<nf;i++){
    st = i*q;
    tgt = nbList[i*(q-1)+2];
    swap(f[st+1],f[tgt*q+3]);
    //if(tgt<nf) swap(f[st+1],f[tgt*q+3]);
    tgt = nbList[i*(q-1)+3];
    swap(f[st+2],f[tgt*q+4]);
    //if(tgt<nf) swap(f[st+2],f[tgt*q+4]);
    tgt = nbList[i*(q-1)+6];
    swap(f[st+5],f[tgt*q+7]);
    //if(tgt<nf) swap(f[st+5],f[tgt*q+7]);
    tgt = nbList[i*(q-1)+7];
    swap(f[st+6],f[tgt*q+8]);
    //if(tgt<nf) swap(f[st+6],f[tgt*q+8]);
  }
}
  
void LB::applyBC(){
  for (int i=0;i<nbc;i++){
    pBC[i].dynf(f,pBC[i].selfData);
  }
}

void LB::writeVelocity(const std::string filename){
  using namespace std;
  ofstream out(filename.c_str(),ios::out | ios::app);
    if (out.is_open()){
      //double rho;
      double ux;
      double uy;
      for (int i=0; i<nf; i++){
        //computeMacros(i, &rho, &ux, &uy);
        ux=v[2*i];
        uy=v[2*i+1];
        out<<setw(6)<<ux<<" "<<setw(6)<<uy<<endl;
      }
    }else{
      cout<<"cannot open output file"<<endl;
    }
}

void LB::writeGeometry(const std::string filename){
  using namespace std;
  ofstream out(filename.c_str(),ios::out | ios::app);
    if (out.is_open()){
      out<<"lx"<<endl;
      out<<setw(10)<<lx<<endl;
      out<<"ly"<<endl;
      out<<setw(10)<<ly<<endl;
      out<<"nf"<<endl;
      out<<setw(10)<<nf<<endl;
      out<<"nbList"<<endl;
        int nadj = q - 1;
        for (int i = 0; i < nf; i++){
            out<<setw(10)<<i;//current node number
            for (int j = 0; j < nadj; j++)
              out << setw(10)<< nbList[i*nadj + j];
            out<<setw(10)<< coor[i*d]<<setw(10)<< coor[i*d + 1]; //x,y
            out<<endl;
          }
        out<<"xy2idx"<<endl;
        for (int row = 0; row < ly; row++){
          for (int col = 0; col < lx; col++)
            {out<<setw(10)<< xy2idx[row*lx + col];}
          out<<endl;
        }
        out<<"boundaries"<<endl;
        out<<setw(10)<<nbc<<endl;
        out<<"velocity"<<endl;
        out<<setw(10)<<nv<<endl;
        for (int i=0;i<nv;i++){
          velData *data =(velData*)pBC[i].selfData; 
          out<<setw(10)<<data->size<<endl;
          out<<setw(10)<<data->norm[0]<<setw(10)<<data->norm[1]<<endl;
          for(int j=0;j<data->size;j++)
            out<<setw(10)<<data->node[j]<<" "<<setw(10)<<data->ux[j]
              <<" "<<setw(10)<<data->uy[j]<<endl;
        }
        out<<"open"<<endl;
        out<<setw(10)<<no<<endl;
        for (int i=0;i<no;i++){
          openData *data =(openData*)pBC[nv+i].selfData; 
          out<<setw(10)<<data->size<<endl;
          out<<setw(10)<<data->norm[0]<<setw(10)<<data->norm[1]<<endl;
          for(int j=0;j<data->size;j++)
            out<<setw(10)<<data->node[j]<<setw(10)<<data->nodeNext[j]<<endl;
        }
        out<<"bounceback"<<endl;
        out<<setw(10)<<nb<<endl;
        for (int i=0;i<nb;i++){
          bbData *data =(bbData*)pBC[nv+no+i].selfData; 
          out<<setw(10)<<data->size<<endl;
          for(int j=0;j<data->size;j++)
            out<<setw(10)<<data->node[j]<<endl;
        }
    }else{
      cout<<"cannot open output file"<<endl;
    }
}

void LB::writeForce(const std::string filename){
  using namespace std;
  ofstream out(filename.c_str(),ios::out | ios::app);
    if (out.is_open()){
      for (int i=0; i<nf; i++){
        out<<setw(6)<<force[2*i]<<" "<<setw(6)<<force[2*i+1]<<endl;
      }
    }else{
      cout<<"cannot open fluid force output file"<<endl;
    }
}
