#ifndef BOUNDARY_2D_H
#define BOUNDARY_2D_H

struct bbData{
  int size;
  int *node;
public:
  bbData(){size=0;node=NULL;}
  bbData(int size_){
    size = size_;
    node = new int[size];
  }
  bbData& operator=(const bbData& rhs){
    if (this !=&rhs){
      delete [] node;
      size = rhs.size;
      node = new int[size];
      for (int i=0;i<size;i++)
        node[i]=rhs.node[i];
    }
    return *this;
  }
  ~bbData(){
    delete [] node;
  }
};
struct velData{
  int size;
  int norm[2];
  int *node;
  double *ux, *uy;
public:
  velData(){size=0;
  norm[0]=0;norm[1]=0;
  node =NULL; ux=NULL; uy=NULL;}
  velData(int size_){
    size = size_;
    node = new int[size];
    ux = new double[size];
    uy = new double[size];
  }
  velData& operator=(const velData& rhs){
    if (this !=&rhs){
      delete [] node;
      delete [] ux;
      delete [] uy;
      size = rhs.size;
      norm[0]=rhs.norm[0]; norm[1]=rhs.norm[1];
      node = new int[size];
      ux = new double[size];
      uy = new double[size];
      for (int i=0;i<size;i++){
        node[i]=rhs.node[i];
        ux[i]=rhs.ux[i];
        uy[i]=rhs.uy[i];
      }
    }
    return *this;
  }
  ~velData(){
    delete [] node;
    delete [] ux;
    delete [] uy;
  }
};
struct openData{
  int size;
  int norm[2];
  int *node;
  int *nodeNext;
public:
  openData(){size=0;
    norm[0]=0;norm[1]=0;
    node=NULL;nodeNext=NULL;
  }
  openData(int size_){
    size = size_;
    node = new int[size];
    nodeNext = new int[size];
  }
  openData& operator=(const openData& rhs){
    if (this !=&rhs){
      delete [] node;
      delete [] nodeNext;
      size = rhs.size;
      norm[0]=rhs.norm[0]; norm[1]=rhs.norm[1];
      node = new int[size];
      nodeNext = new int[size];
      for (int i=0;i<size;i++){
        node[i]=rhs.node[i];
        nodeNext[i]=rhs.nodeNext[i];
      }
    }
    return *this;
  }
 
  ~openData(){
    delete [] node;
    delete [] nodeNext;
  }
};

void bounceBack(double* f, void* selfData);

// ZouHe velocity boundaries on upper, lower, left and right
//   boundaries
void leftZouHe (double* f, void* selfData);
void rightZouHe (double* f, void* selfData);
void upperZouHe (double* f, void* selfData);
void lowerZouHe (double* f, void* selfData);

void leftOpen (double* f, void* selfData);
void rightOpen (double* f, void* selfData);
void upperOpen (double* f, void* selfData);
void lowerOpen (double* f, void* selfData);

/*
struct BC{
  int size;
  int norm[2];
  int *node;
  double *ux, *uy;
  funType dynf;
public:
  BC(){
    
    size = 2;
    norm[0]=0;
    norm[1]=0;
    node = new int[size];
    ux = new double[size];
    uy = new double[size];
    for (int i=0;i<size;i++){
      node[i]=0;ux[i]=0;uy[i]=0;
   //   node[i]=i;ux[i]=i;uy[i]=i;
    }
    
  }
  BC(const int size_){
    size = size_;
    norm[0]=0;
    norm[1]=0;
    node = new int[size];
    ux = new double[size];
    uy = new double[size];
    for (int i=0;i<size;i++){
      node[i]=0;ux[i]=0;uy[i]=0;
      //node[i]=size;ux[i]=size;uy[i]=size;
    }
    dynf = 0;
  }
  ~BC(){
      
    delete [] node;
    delete [] ux;
    delete [] uy;
    std::cout<<"destructed"<<std::endl;
  }
  BC& operator=(const BC& rhs){
    if (this !=&rhs){
      delete [] node;
      delete [] ux;
      delete [] uy;
      size = rhs.size;
      node = new int[size];
      ux = new double[size];
      uy = new double[size];
      for (int i=0;i<size;i++){
        node[i]=rhs.node[i];
        ux[i]=rhs.ux[i];
        uy[i]=rhs.uy[i];
      }
      dynf = rhs.dynf;
    }
    return *this;
  }
};
*/
/*
void bounceBack(double* f, int id, void* selfData);

// ZouHe velocity boundaries on upper, lower, left and right
//   boundaries
void leftZouHe (double* f, int id, void* selfData);

void rightZouHe (double* f, int id, void* selfData);

void upperZouHe (double* f, int id, void* selfData);

void lowerZouHe (double* f, int id, void* selfData);

void rightOpen (double* f, int id, void* selfData);
*/
/*template <class T>
void velocityZouHe(double *f, int id, void*selfData);
*/
//void velocityZouHe(double* f, int id, int *n, void* selfData);
#endif
