#include <iostream>
#include "boundary.h"

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
  Base *pb;
  Child child(1,2);
  pb = &child;
  pb->out();
  for (int i=0;i<5;i++)
    cout<<"x "<<pb->x[i]<<endl;
  return 0;
}
