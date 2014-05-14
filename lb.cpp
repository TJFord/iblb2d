# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>
# include "lb.h"
# include "units.h"
# include "boundary.h"

using namespace std;

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
  lx=0;ly=0;nf=0;
  IJidx=NULL;
  nbList=NULL;
  f=NULL;
  coor=NULL;
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
    delete [] IJidx;
    delete [] nbList;
    delete [] f;
    delete [] coor;
    delete [] pBC;
    delete [] pVel;
    delete [] pOpen;
    delete [] pBB;
    delete [] force;
    delete pUnits;

    lx=rhs.lx; ly=rhs.ly; nf=rhs.nf;
    omega = rhs.omega;
    nbc=rhs.nbc; nv=rhs.nv; no=rhs.no; nb=rhs.nb;
    collisionFun = rhs.collisionFun;
   
    int tmp;
    IJidx = new int[lx*ly];
    tmp = lx*ly;
    for (int i=0;i<tmp;i++)
      IJidx[i]=rhs.IJidx[i];
    
    nbList = new int[nf*(q-1)];
    tmp = nf*(q-1);
    for (int i=0;i<tmp;i++)
      nbList[i]=rhs.nbList[i];

    f = new double[(nf+1)*q];
    tmp = (nf+1)*q;
    for (int i=0;i<tmp;i++)
      f[i]=rhs.f[i];

    coor = new int[nf*d];
    tmp = nf*d;
    for (int i=0;i<tmp;i++)
     coor[i]=rhs.coor[i];
    
    force = new double[nf*d];
    tmp = nf*d;
    for (int i=0;i<tmp;i++)
     force[i]=rhs.force[i];
    
    /*
    for (int i=0;i<nbc;i++){
      pBC[i].dynf=rhs.pBC[i].dynf;
      *(pBC[i].selfData) = *(rhs.pBC[i].selfData);
    }*/
    
    pVel = new velData[nv];
    for(int i=0;i<nv;i++){
      pVel[i].size=rhs.pVel[i].size;
      pVel[i].norm[0]=rhs.pVel[i].norm[0];
      pVel[i].norm[1]=rhs.pVel[i].norm[1];
      for(int j=0;j<pVel[i].size;j++){
        pVel[i].node[j]=rhs.pVel[i].node[j];
        pVel[i].ux[j]=rhs.pVel[i].ux[j];
        pVel[i].uy[j]=rhs.pVel[i].uy[j];
      }
    }
    pOpen = new openData[no];
    for(int i=0;i<no;i++){
      pOpen[i].size=rhs.pOpen[i].size;
      pOpen[i].norm[0]=rhs.pOpen[i].norm[0];
      pOpen[i].norm[1]=rhs.pOpen[i].norm[1];
      for(int j=0;j<pOpen[i].size;j++){
        pOpen[i].node[j]=rhs.pOpen[i].node[j];
        pOpen[i].nodeNext[j]=rhs.pOpen[i].nodeNext[j];
      }
    }
    pBB = new bbData[nb];
    for(int i=0;i<nb;i++){
      pBB[i].size=rhs.pBB[i].size;
      for(int j=0;j<pBB[i].size;j++){
        pBB[i].node[j]=rhs.pBB[i].node[j];
      }
    }
    pBC = new Dynamics[nbc];
    tmp=0;
    for (int i=0;i<nv;i++){
      pBC[tmp].selfData = (void*)&pVel[i];
      pBC[tmp].dynf = rhs.pBC[tmp].dynf;
      tmp++;
    } 
    for (int i=0;i<no;i++){
      pBC[tmp].selfData = (void*)&pOpen[i];
      pBC[tmp].dynf = rhs.pBC[tmp].dynf;
      tmp++;
    } 
    for (int i=0;i<nb;i++){
      pBC[tmp].selfData = (void*)&pBB[i];
      pBC[tmp].dynf = rhs.pBC[tmp].dynf;
      tmp++;
    } 
    pUnits = new Units;
    *pUnits = *(rhs.pUnits);
  }  
}

void LB::readInput(const std::string filename)
{	
  ifstream in(filename.c_str(),ios::in);
  if (in.is_open())
    {
     // cout<<"I get here"<<endl;
      string str; 
      int tmp,size(0);
      pUnits = new Units;
      while (!in.eof()){
        in>>str;
        if(str.compare("dx")==0){
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
        else if (str.compare("nf") == 0)
        { in >> nf; f = new double[q*(nf+1)];}
        else if (str.compare("nbList") == 0){
          nbList = new int[nf*(q-1)];
          coor = new int[nf*d];
          force = new double[nf*d];
          int nadj = q - 1;
          for (int i = 0; i < nf; i++){
            in >> tmp;//current node number
            for (int j = 0; j < nadj; j++){
              in >> nbList[i*nadj + j];
            }
            in >> coor[i*d] >> coor[i*d + 1]; //x,y
          }
        }else if(str.compare("IJidx") == 0){
          IJidx = new int[lx*ly];
          for (int j = 0; j < lx; j++)
            for (int i = 0; i < ly; i++)
              in >> IJidx[j*ly + i];
        }else if(str.compare("boundaries")==0){
          in >> tmp;
          for (int m=0;m<tmp;m++){
            in >> str;
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
              //in >> pBB[0].size;
              in>>size; pBB[0]=bbData(size);
              nbc++;
              for (int j = 0; j < pBB[0].size; j++){
                in >> pBB[0].node[j];
              }
            }
          }
        }else{
          if (in.eof())
            cout<<"succesfully read file"<<endl;
          else{ 
            cout<<"wrong keywords in input! "<<str<<endl;
            break;
          }
        }
      }

      //pUnits->calculateLBPara();

      pBC = new Dynamics[nbc];
      tmp = 0;
      for(int i=0;i<nv;i++){
        pBC[tmp].selfData = (void*)&pVel[i];
        if (pVel[i].norm[0]==0){
          if (pVel[i].norm[1]==1){
            pBC[tmp].dynf = &upperZouHe;
            tmp++;
          }else if(pVel[i].norm[1]==-1){
            pBC[tmp].dynf = &lowerZouHe;
            tmp++;
          }else
            cout<<"wrong norm"<<endl;
        }else if(pVel[i].norm[0]==1){
          pBC[tmp].dynf = &rightZouHe;
          tmp++;
        }else if(pVel[i].norm[0]==-1){
          pBC[tmp].dynf = &leftZouHe;
          tmp++;cout<<"leftZouHe"<<endl;
        }
      }
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
          tmp++;
        }
      }

      pBC[tmp].selfData = (void*)&pBB[0];
      pBC[tmp].dynf = &bounceBack;

      in.close();
     } else {
      cout<<"cannot open input file"<<endl;
    }
}

LB::~LB(){
  delete [] IJidx;
  delete [] nbList;
  delete [] f;
  delete [] coor;
  delete [] pVel;
  delete [] pOpen;
  delete [] pBB;
  delete [] pBC;
  delete [] force;
  delete pUnits;

  cout<<"Destroy LB"<<endl;
}

void LB::printInfor(){
  cout<<"********************************"<<endl;
  //cout.precision(5);
  cout<<"IBLB program"<<endl;
  cout<<"fluid:"<<endl;
  cout<<"length L = "<<lx*pUnits->dx<<" m, width W = "<<ly*pUnits->dx<<" m"<<endl;
  cout<<"density rho = "<<pUnits->rho<<" kg/m^3"<<endl;
  cout<<"viscosity nu = "<<pUnits->nu<<" m^2/s"<<endl;
  cout<<"Re = "<<pUnits->Re<<endl;
  cout<<"LB:"<<endl;
  cout<<"tau = "<<pUnits->tau<<endl;
  cout<<"u_lb = "<<pUnits->u_lb<<endl;
  cout<<"dx = "<<pUnits->dx<<" m"<<endl;
  cout<<"dt = "<<pUnits->dt<<" s"<<endl;
  if (pUnits->dt == 0.0) cout<<"oh shit"<<endl;
  cout<<"********************************"<<endl;
}
  
void LB::init(){
  pUnits->calculateLBPara();
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
    }
    //force[i*d] = 1.e-6;
    //force[i*d+1] =0.0;
    //cout<<"force"<<"i="<<i<<" "<<force[i*d]<<" "<<force[i*d+1]<<endl;
    st = i*q;
    for (int j=0; j<q; j++){
      f[st+j] = computeEquilibrium(j, 1.0, 0.0, 0.0, 0.0); //initialize with zero velocity          
    }
  }
  st = nf*q;
  for (int j=0;j<q;j++)
    f[st+j]=0.0; //extra node, doesn't matter

  //collisionFun=&LB::bgk;//fluid scheme
  //collisionFun=&LB::regularized;//fluid scheme
  collisionFun=&LB::stokes;//fluid scheme
  
  if( nv !=0 ){
    velData *data =(velData*) pBC[0].selfData;
    for (int i=0;i<data->size;i++)
      cout<<"ux"<<data->ux[i]<<endl;

    cout<<"after"<<endl;
    double velscale = pUnits->dx/pUnits->dt;
    for (int i=0;i<pVel[0].size;i++)
    pVel[0].ux[i] /= velscale;

    velData *data1 =(velData*) pBC[0].selfData;
    for (int i=0;i<data1->size;i++)
      cout<<"ux"<<data1->ux[i]<<endl;
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
}

void LB::stokes(int id, void* selfData){
  double rho, ux, uy;
  computeMacros(id,&rho, &ux, &uy);
  int st = id*q;
  for (int i=0;i<q;i++){
    f[st+i] *= (1-omega);
    f[st+i] += omega*computeEqStokes(i,rho, ux, uy);
  }
}
//void LB::applyForce(int id, double fx, double fy){
void LB::applyForce(){
  double F[9];
  double rho, ux, uy;
  double fx,fy;
  int st=0;
  int id;
  for (id=0;id<nf;id++){
    st = id*9;
    fx = force[id*d];
    fy = force[id*d+1];
    computeMacros(id,&rho, &ux, &uy);
    F[0]=(1.-0.5*omega)*w[0]*(3.*((  -ux)*fx +(  -uy)*fy));
    F[1]=(1.-0.5*omega)*w[1]*(3.*(( 1.-ux)*fx +(  -uy)*fy)+9.*(ux*fx));
    F[2]=(1.-0.5*omega)*w[2]*(3.*((  -ux)*fx +( 1.-uy)*fy)+9.*(uy*fy));
    F[3]=(1.-0.5*omega)*w[3]*(3.*((-1.-ux)*fx +(  -uy)*fy)+9.*(ux*fx));
    F[4]=(1.-0.5*omega)*w[4]*(3.*((  -ux)*fx +(-1.-uy)*fy)+9.*(uy*fy));
    F[5]=(1.-0.5*omega)*w[5]*(3.*(( 1.-ux)*fx +( 1.-uy)*fy)+9.*((ux+uy)*(fx+fy)));
    F[6]=(1.-0.5*omega)*w[6]*(3.*((-1.-ux)*fx +( 1.-uy)*fy)+9.*((ux-uy)*(fx-fy)));
    F[7]=(1.-0.5*omega)*w[7]*(3.*((-1.-ux)*fx +(-1.-uy)*fy)+9.*((ux+uy)*(fx+fy)));
    F[8]=(1.-0.5*omega)*w[8]*(3.*(( 1.-ux)*fx +(-1.-uy)*fy)+9.*((ux-uy)*(fx-fy)));
    for (int i=0;i<9;i++)
      f[st+i] += F[i];
  }
}
/*swap algorithm, takes time*/
void LB::stream()
{/*
  
  double *ftmp;
  ftmp = new double[nn*q];
  for (int i=0;i<nn*q;i++)
          ftmp[i]=f[i];
  for (int i=0;i<nf;i++)
  {
    int st = i*q;
    for (int j=1;j<q;j++)
    {
      int a = i*(q-1)+j-1;
      int idx = nList[a];
      if (idx > nn)
        cout<<"i"<<i<<" "<<"j"<<j<<" "<<"a"<<a<<" "<< nList[4]<<" "<<idx<<endl;
            
      else
        f[idx*q+j]=ftmp[st+j];
    }
  }
  delete [] ftmp;
*/
  /*
int st = 0;
int tgt = 0;
  for (int i=0; i< nf; i++)
  { 
    st = i*q;
    swap(f[st+1], f[st+3]);
    swap(f[st+2], f[st+4]);
    swap(f[st+5], f[st+7]);
    swap(f[st+6], f[st+8]);    
  } 
  for (int i=0;i<nf;i++)
  {
	  st = i*q; //f position, not node
	
	  {

		tgt = nList[i*(q-1)+2];
		if (tgt < nf) swap(f[st+1],f[tgt*q+3]);
		tgt = nList[i*(q-1)+3];
		if (tgt < nf) swap(f[st+2],f[tgt*q+4]);
		tgt = nList[i*(q-1)+6];
		if (tgt < nf) swap(f[st+5],f[tgt*q+7]);
		tgt = nList[i*(q-1)+7];
		if (tgt < nf) swap(f[st+6],f[tgt*q+8]);
		
		
	   }    
  }*/
}

void LB::collide()
{
  /*dynType pf;
  for (int i=0;i<nf;i++)  
    {
      pf = pDyn[i].dyn;
      (this->*pf)(i,pDyn[i].selfData);
    }
	*/
  for (int i=0;i<nf;i++)
  {
    bgk(i,(void *)&omega);
  }
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
    if(tgt<nf) swap(f[st+1],f[tgt*q+3]);
    tgt = nbList[i*(q-1)+3];
    if(tgt<nf) swap(f[st+2],f[tgt*q+4]);
    tgt = nbList[i*(q-1)+6];
    if(tgt<nf) swap(f[st+5],f[tgt*q+7]);
    tgt = nbList[i*(q-1)+7];
    if(tgt<nf) swap(f[st+6],f[tgt*q+8]);
  }
}
  
void LB::applyBC(){
  for (int i=0;i<nbc;i++){
    pBC[i].dynf(f,pBC[i].selfData);
  }
}

void LB::output(const std::string filename)
{
  ofstream out(filename.c_str(),ios::out | ios::app);
    if (out.is_open())
    {/*
      out<<setw(10)<<"lx"<<setw(10)<<lx<<endl;
      out<<setw(10)<<"ly"<<setw(10)<<ly<<endl;
      out<<setw(10)<<"nf"<<setw(10)<<nf<<endl;
      out<<setw(10)<<"nbList"<<endl;
        int nadj = q - 1;
        for (int i = 0; i < nf; i++){
            out<<setw(10)<<i;//current node number
            for (int j = 0; j < nadj; j++)
              out << setw(10)<< nbList[i*nadj + j];
            out<<setw(10)<< coor[i*d]<<setw(10)<< coor[i*d + 1]; //x,y
            out<<endl;
          }
        out<<setw(10)<<"IJidx"<<endl;
        for (int j = 0; j < lx; j++){
          for (int i = 0; i < ly; i++)
            {out<<setw(10)<< IJidx[j*ly + i];}
          out<<endl;
        }
        out<<"boundaries"<<setw(10)<<nbc<<endl;
        out<<"velocity"<<setw(10)<<nv<<endl;
        for (int i=0;i<nv;i++){
          velData *data =(velData*)pBC[i].selfData; 
          out<<"size"<<setw(10)<<data->size<<endl;
          out<<setw(10)<<data->norm[0]<<setw(10)<<data->norm[1]<<endl;
          for(int j=0;j<data->size;j++)
            out<<setw(10)<<data->node[j]<<setw(10)<<data->ux[j]
              <<setw(10)<<data->uy[j]<<endl;
        }
        out<<"open"<<setw(10)<<no<<endl;
        for (int i=0;i<no;i++){
          openData *data =(openData*)pBC[nv+i].selfData; 
          out<<"size"<<setw(10)<<data->size<<endl;
          out<<setw(10)<<data->norm[0]<<setw(10)<<data->norm[1]<<endl;
          for(int j=0;j<data->size;j++)
            out<<setw(10)<<data->node[j]<<setw(10)<<data->nodeNext[j]<<endl;
        }
        out<<"bounceback"<<setw(10)<<nb<<endl;
        for (int i=0;i<nb;i++){
          bbData *data =(bbData*)pBC[nv+no+i].selfData; 
          out<<"size"<<setw(10)<<data->size<<endl;
          for(int j=0;j<data->size;j++)
            out<<setw(10)<<data->node[j]<<endl;
        }*/
// this is used to save velocity

      double rho;
      double ux;
      double uy;
      for (int i=0; i<nf; i++)    
//	  for (int i=0; i<ninlet; i++)
        {
      //    computeMacros(inList[i], &rho, &ux, &uy);
        //  cout<<setw(10)<<ux<<" "<<rho<<endl;
	   computeMacros(i, &rho, &ux, &uy);
          out<<setw(6)<<ux<<" "<<setw(6)<<uy<<endl;
        }

        out.close();
    }
    else{
      cout<<"cannot open output file"<<endl;
    }
		
}
