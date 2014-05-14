# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>
# include "units.h"

using namespace std;

Units::Units(){
  dx=0;dt=0;dm=0;
  vis=0; u=0;

  omega=0; tau=0;
  nu=0; Re=0; u_lb=0;
}

Units::Units(const double dx_, const double tau_, const double rho_, const double vis_,
const double u_)
{
  dx  = dx_;
  tau = tau_;
  rho = rho_;
  vis = vis_;
  u   = u_;
}

void Units::setParameters(const std::string filename){
  ifstream in(filename.c_str(),ios::in);
  string str;
  if (in.is_open()){
    while(!in.eof()){
      in>>str;
      if (str.compare("dx")==0)
        in >> dx;
      else if (str.compare("dt")==0)
        in >> dt;
      else if (str.compare("dm")==0)
        in >> dm;
      else if (str.compare("tau")==0)
        in >> tau;
      else if (str.compare("vis")==0)
        in >> vis;
      else 
        cout<<"unknown keywords in units input"<<endl;
    }
    in.close();
  }else{
    cout<<"cannot open Units input file"<<endl;
  }
}
void Units::calculateLBPara()
{
  if (tau !=0.){
    omega = 1./tau;
    nu = (tau - 0.5)/3.;
    dt = nu*dx*dx/vis;
    //cout<<"Units:dt="<<dt<<endl;
    u_lb = u*dt/dx;
    Re = u_lb/nu;
    dm = rho*dx*dx*dx; 
  }else if(dt != 0.){
    nu = vis*dx*dx/dt;
    tau = 3*nu + 0.5;
    u_lb = u*dt/dx;
    Re = u_lb/nu;
    dm = rho*dx*dx*dx;
    omega = 1./tau;
  }
  //cout<<"Re"<<dx*u/vis<<endl;
}
