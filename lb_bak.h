#ifndef LB_2D_H
#define LB_2D_H

# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>

# include "boundary.h"

using namespace std;

typedef void (*dynfType)(double* f, int id, void* selfData);
//typedef void (*dynfType)(double* f, int id, int *n, void* selfData);
typedef struct
{
	//void (*dynfType) *dynf;
	dynfType dynf;
	void* selfData;
} Dynamics;

class LB
{
public:
/*	
	double ts; //timestep
	double Re; //Reynold's number
	double nu; //viscosity
	*/
	double omega; // relaxation = 1/tau;

	static const int d;//dimension
	static const int q;//discrete velocity component
	static const double cs2; // square of sound speed

	static const double c[9][2]; //velocity vector
	static const double w[9];// weight factor

	int nn; // total node number, include cushion nodes
	int nf, ninlet, noutlet, nbb; // fluid nodes,inlet node

	double *f; //density distribution
	int *nList, *inList, *outList, *bbList; //adjacent,inlet,outlet,bounceback node list
	//int *ijIdx; //for easy index
	int *coor; // LB coordinates for each node
	int inletNormVec[2];
	int outletNormVec[2];
	
	Dynamics *pDyn;

public:
	//LB(fstream &in); //only read nn, nlist, etc.
	LB(const string in); //only read nn, nlist, etc.
	~LB();

	void init(bcData *inletVel);//initialize null pointer
	void stream();
	void collide();

	void collideSwap();
	void streamSwap();
//	template<typename F, typename T>
//	void applyBC(int id, F f, T t); 
	void applyBC();

	void bgk(int id, void* selfData);
//	void regularized();

	void computeMacros(int id, double *rho, double * ux, double *uy);
	double computeEquilibrium(int idx, double rho,double ux, double uy,double uSqr);
	int getBCNum();
	
	void output(const string out);
//	void saveVel(fstream &out);//save file
	//helper function
	inline int* getCoor(int id){return &coor[d*id];}
	inline int getNinlet(){return ninlet;}
	inline int getInList(int id){return inList[id];}
};

#endif