# include <iostream>
# include <fstream>
# include <iomanip>
# include <string>
# include "lb.h"
# include "boundary.h"

using namespace std;


const int LB::d = 2;
const int LB::q = 9;
const double LB::cs2 = 1.0/3.0;

const double LB::c[9][2] = {
    {0,0},
    {1,0}, {0,1}, {-1,0}, {0,-1},
    {1,1}, {-1,1}, {-1,-1}, {1,-1}
};

const double LB::w[9] = { 4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

void swap(double &f1, double &f2)
{
	double tmp = f1;
	f1 = f2;
	f2 = tmp;
}

LB::LB(const string filename)
{	
	ifstream in(filename.c_str(),ios::in);
  if (in.is_open())
    {
      omega = 1.0; //omega initialize
      in >> nf;
      nn = nf + 1;

      int tmp;
      int nadj = q - 1;
		
      f = new double[nn*q];
      nList = new int[nf*(q-1)];
      coor = new int[nf*d];

      for (int i = 0; i < nf; i++)
    	{
    	  in >> tmp;//current node number
    	  for (int j = 0; j < nadj; j++)
    	    in >> nList[i*nadj + j];
    	  in >> coor[i*d] >> coor[i*d + 1]; //x,y
    	}

      string str; 

      /* read inlet node list */
      in >> str;
      if (str.compare("inlet") == 0) {
      	in >> ninlet;
		in >> inletNormVec[0] >>inletNormVec[1]; //boundary surface norm
      	inList = new int[ninlet];
      	for (int i = 0; i < ninlet; i++)
      	  in >> inList[i];
      }
		
      in >> str;
      if (str.compare("outlet") == 0) {
      	in >> noutlet;
		in >> outletNormVec[0] >>outletNormVec[1];//boundary surface norm
      	outList = new int[noutlet];
      	for (int i = 0; i < noutlet; i++)
      	  in >> outList[i];
      }

      in >> str;
      if (str.compare("bounceback") == 0) {
      	in >> nbb;
      	bbList = new int[nbb];
      	for (int i = 0; i < nbb; i++)
      	  in >> bbList[i];
      }

	  pDyn = new Dynamics[nf];

      in.close();

      /* assign variables */
      //pNode = new Node[nf];
      //pDyn = new Dynamics[nf];
      // inletVelData = new VelocityBCData[ninlet];
		
    } else {
      cout<<"cannot open input file"<<endl;
    }
	
}

LB::~LB()
{
  delete [] nList;
  delete [] f;
  delete [] coor;
  delete [] inList;
  delete [] outList;
  delete [] bbList;
  
//  delete [] pNode;
  delete [] pDyn;
  //delete [] inletVelData;

  cout<<"Destroy LB"<<endl;
}

void LB::init(bcData *inletVel)
{
  int st = 0;
  for (int i = 0; i < nf; i++)
    {
		pDyn[i].dynf = 0;
		pDyn[i].selfData = 0;
      st = i*q;
      for (int j=0; j<q; j++) // distribution function
        {
          f[st+j] = computeEquilibrium(j, 1.0, 0.0, 0.0, 0.0); //initialize with zero velocity          
        }
    }

    st = nf*q;
    for (int j=0;j<q;j++)
      f[st+j]=0.0; //extra node, doesn't matter

	if (inletNormVec[0]==0){//y
          if (inletNormVec[1] == 1){//positive y
            for (int i=0;i<ninlet;i++)
            {
              pDyn[inList[i]].dynf = &upperZouHe;
              pDyn[inList[i]].selfData = (void*)&inletVel[i];
            }
          }else{//negative y
            for (int i=0;i<ninlet;i++)
            {
              pDyn[inList[i]].dynf = &lowerZouHe;
              pDyn[inList[i]].selfData = (void*)&inletVel[i];
            }
          }
	}else if(inletNormVec[0]==1){//positive x
          for (int i=0;i<ninlet;i++)
            {
              pDyn[inList[i]].dynf = &rightZouHe;
              pDyn[inList[i]].selfData = (void*)&inletVel[i];
            }
	}else{//negative x
          for (int i=0;i<ninlet;i++)
            {
              pDyn[inList[i]].dynf = &leftZouHe;
              pDyn[inList[i]].selfData = (void*)&inletVel[i];
            }
	}



/*	for (int i=0;i<ninlet;i++)
	{
		pDyn[inList[i]].dynf = &leftZouHe;
		pDyn[inList[i]].selfData = (void*)&inletVel[i];
	}
*/
	/*
	for (int i=0;i<noutlet;i++)
	{
		pDyn[outList[i]].dynf = &rightZouHe;
		pDyn[outList[i]].selfData = (void*)&inletVel[i];
	}
	*/
	
	for (int i=0;i<noutlet;i++)
	{
		pDyn[outList[i]].dynf = &rightOpen;
		pDyn[outList[i]].selfData = (void*)&nList[outList[i]*(q-1)+2];
//		cout<<"left right"<<nList[outList[i]*(q-1)+2]<<endl;
	}
	

	for (int i=0;i<nbb;i++)
	{
		pDyn[bbList[i]].dynf = &bounceBack;
		pDyn[bbList[i]].selfData = 0;
	}

}

void LB::computeMacros(int id, double *rho, double * ux, double *uy)
{
  int st = id*q; //save computation time
  double upperLine = f[st+2] + f[st+5] + f[st+6];
  double mediumLine = f[st] + f[st+1] + f[st+3];
  double lowerLine = f[st+4] + f[st+7] + f[st+8];
  *rho = upperLine + mediumLine + lowerLine;
  *ux = (f[st+1]+f[st+5]+f[st+8] -(f[st+3] +f[st+6] + f[st + 7]))/(*rho);
  *uy = (upperLine - lowerLine)/(*rho);

}


double LB::computeEquilibrium(int idx, double rho,double ux, double uy,double uSqr)
{
  double c_u = c[idx][0]*ux + c[idx][1]*uy;
  return rho*w[idx]*(1. + 3.*c_u + 4.5*c_u*c_u - 1.5*uSqr);

}

void LB::bgk(int id, void* selfData)
{
  //cout<<"bgk"<<endl;
  double omega = *((double *)selfData);
  double rho, ux, uy;
  computeMacros(id,&rho, &ux, &uy);
  double uSqr = ux*ux + uy*uy;
  int st = id*q;
  for (int i=0;i<q;i++)
    {
      f[st+i] *= (1-omega);
      f[st+i] += omega*computeEquilibrium(i,rho, ux, uy, uSqr);
  }
  
}

/*swap algorithm, takes time*/
void LB::stream()
{	
	double *ftmp;
	ftmp = new double[nn*q];
	for (int i=0;i<nn*q;i++)
		ftmp[i]=f[i];
	for (int i=0;i<nf;i++)
	{
		int st = i*q;
		for (int j=1;j<q;j++)
		{
			int idx = nList[i*(q-1)+j-1];
			if (idx > nn)
				cout<<"i"<<i<<" "<<"j"<<j<<" "<<idx<<endl;
			f[idx*q+j]=ftmp[st+j];
//			cout<<f[idx*q+1]<<endl;
		}
	}
	delete [] ftmp;
	

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

void LB::collideSwap()
{
	int st=0;
	for (int i=0;i<nf;i++)
	  {
		bgk(i,(void *)&omega);
		st = i*q;
		swap(f[st+1], f[st+3]);
		swap(f[st+2], f[st+4]);
		swap(f[st+5], f[st+7]);
		swap(f[st+6], f[st+8]);   
	  }
}
void LB::streamSwap()
{
	int st=0;
	int tgt = 0;
	for (int i=0;i<nf;i++)
  {
	  st = i*q; //f position, not node	
	  {

/*		tgt = nList[i*(q-1)+2];
		if (tgt < nf) swap(f[st+1],f[tgt*q+3]);
		tgt = nList[i*(q-1)+3];
		if (tgt < nf) swap(f[st+2],f[tgt*q+4]);
		tgt = nList[i*(q-1)+6];
		if (tgt < nf) swap(f[st+5],f[tgt*q+7]);
		tgt = nList[i*(q-1)+7];
		if (tgt < nf) swap(f[st+6],f[tgt*q+8]);*/
		
		tgt = nList[i*(q-1)+2];
		swap(f[st+1],f[tgt*q+3]);
		tgt = nList[i*(q-1)+3];
		swap(f[st+2],f[tgt*q+4]);
		tgt = nList[i*(q-1)+6];
		swap(f[st+5],f[tgt*q+7]);
		tgt = nList[i*(q-1)+7];
		swap(f[st+6],f[tgt*q+8]);
	   }    
  }
}

int LB::getBCNum()
{
	return ninlet;
}

void LB::applyBC()
{
/*	for (int i=0;i<nf;i++)
		pDyn[i].dynf(f,i,pDyn[i].selfData);
		*/
	for (int i=0;i<ninlet;i++)
	{
		pDyn[inList[i]].dynf(f,inList[i],pDyn[inList[i]].selfData);
	}
	for (int i=0;i<noutlet;i++)
	{
		pDyn[outList[i]].dynf(f,outList[i],pDyn[outList[i]].selfData);
	}
	for (int i=0;i<nbb;i++)
	{
		pDyn[bbList[i]].dynf(f,bbList[i],pDyn[bbList[i]].selfData);
	}
}

//void LB::output(ofstream & out)
void LB::output(const string filename)
{
	ofstream out(filename.c_str(),ios::out | ios::app);
    if (out.is_open())
    {
		
      double rho;
      double ux;
      double uy;
      for (int i=0; i<nf; i++)    
	 // for (int i=0; i<ninlet; i++)
        {
        //  computeMacros(inList[i], &rho, &ux, &uy);
        //  cout<<setw(10)<<ux<<" "<<rho<<endl;
		  computeMacros(i, &rho, &ux, &uy);
          out<<setw(6)<<ux<<" "<<setw(6)<<uy<<endl;
	//	  if (i==181)
	//		  cout<<ux<<" uy "<<uy<<endl;
        }

	/*	int nadj = q -1;
      for (int i = 0; i < nf; i++)			
  		{
  		  out<<setw(10)<< i ;
  		  for (int j = 0; j < nadj; j++)
  			out<<setw(10)<<nList[i*nadj+j];	
  		  out<<setw(10)<<coor[i*d]<<setw(10)<<coor[i*d+1]<<endl;			
  		}

      for (int i = 0; i < ninlet; i++)
		out<<setw(10)<< inList[i];
      out<<endl;

      for (int i = 0; i < noutlet; i++)
		out<<setw(10)<< outList[i];
      out<<endl;

      for (int i = 0; i < nbb; i++)
		out <<setw(10)<< bbList[i];
      out<<endl;

      
*/
	out.close();
    }
    else{
      cout<<"cannot open output file"<<endl;
    }
		
}
