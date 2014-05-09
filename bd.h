#ifndef BOUNDARY_2D_H
#define BOUNDARY_2D_H

typedef struct{
	double ux, uy;
//	int nx,ny;
} bcData;
void bounceBack(double* f, int id, void* selfData);

  // ZouHe velocity boundaries on upper, lower, left and right
  //   boundaries
void leftZouHe (double* f, int id, void* selfData);

void rightZouHe (double* f, int id, void* selfData);

void upperZouHe (double* f, int id, void* selfData);

void lowerZouHe (double* f, int id, void* selfData);

void rightOpen (double* f, int id, void* selfData);

/*template <class T>
void velocityZouHe(double *f, int id, void*selfData);
*/
//void velocityZouHe(double* f, int id, int *n, void* selfData);
#endif
