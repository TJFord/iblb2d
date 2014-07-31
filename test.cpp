#include <iostream>
//#include "boundary.h"
//#include <random>
#include "random_mars.h"

using namespace std;
/*
struct A{
  int *f;
public:
  A(){
    f = new int;    
  }
  void init(int sz)
  {
    //int sz=5;
    *f = sz;
  }
  void output(){cout<<"f= "<<*f<<endl;}
  ~ A(){delete f;}
};*/

class Base{
  public:
    int a;
    int *x;
  public:
    Base(int a_){a=a_;x=NULL;}
    virtual ~Base(){//delete x;
      cout<<"base destructor"<<endl;
    }
    virtual void out(){cout<<a*a<<endl;}
};

class Child: public Base{
  public:
    int b;
  public:
    Child(int a_, int b_):Base(a_){
      b=b_; 
      x=new int[5];
      for (int i=0;i<5;i++)
        x[i]=i;
    }
    ~Child(){delete x;
      cout<<"child destructor"<<endl;
    }
    void out(){cout<<b*b<<endl;}
};

int main()
{
  /*Base *pb;
  Child child(1,2);
  pb = &child;
  pb->out();
  for (int i=0;i<5;i++)
    cout<<"x "<<pb->x[i]<<endl;
*/
  /*
  const int nrolls=10000;
  const int nstars=100;

  std::default_random_engine generator;
  std::normal_distribution<double> distribution(5.0,2.0);
  RanMars *random;
  random = new RanMars(5);
  int p[10]={};
  double min, max;
  min=0;
  max=0;
  for (int i=0; i<nrolls;++i){
    //double number = distribution(generator);
    double number = random->gaussian();//distribution(generator);
    if (number > max) max=number;
    if (number < min) min=number;
    if ((number >=0.0) &&(number<10.0)) ++p[int(number)];
  }
  std::cout<<"normal distribution (5.0,2.0):"<<std::endl;
  std::cout<<"max, min"<<max<<" "<<min<<std::endl;
  for (int i=0;i<10;++i){
    std::cout<<i<<"-"<<(i+1)<<":";
    std::cout<<std::string(p[i]*nstars/nrolls,'*')<<std::endl;
  }
*/
  for (int i=0;i<5;i++)
    for (int j=0;j<5;j++){
      if (i==0 || i==4) {
        if (j==1 || j==2) continue;
      }
      std::cout<<"i,j="<<i<<" "<<j<<std::endl;
    }
  /*
  int a=3;
  double c=3.0;
  int b=2;
  int d;
  std::cout<<"a/b "<<a/b<<std::endl;
  std::cout<<"c/b "<<c/b<<std::endl;
  d = c/b*b;
  cout<<"d "<<d<<endl;*/

  return 0;
}
