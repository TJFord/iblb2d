#include <iostream>
#include "boundary.h"

using namespace std;

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
};

int main()
{
  A a;
  a.output();
  a.init(5);
  a.output();
  return 0;
}
