#ifndef RT_H
#define RT_H

#include "primitives.h"
#include "edges.h"
#include "integrateQ1.h"
#include "integrateQ0.h"
#include "laguerre.h"


struct Integration_result2{
	std::list<double> val;
	std::list<int> colInd;
	std::list<int> rowInd;
	int MatrixSize;
	Integration_result2();
	void set_matrix(int col,int row,double value);	
};

/**
 * @brief      Class for CGAL Triangulation.
 */
class Triangulation2{
		Laguerre * Lag;
		IntegrationMethod * Integ;
		bool ImageInstanciated;
		bool CoordInstanciated;
		std::vector<Weighted_point> points;
		double Tolx,Toly;
	public:
		Tesselation Tess;
		int nDir;
		Triangulation2();
		~Triangulation2(){
			delete Lag;
			delete Integ;
		}
		// SETTERS
		void setImage(int Nx,int Ny,bool isQ1, int lenw, double* w,double l,double r,double b,double u);
		void setCoordinates(int lenx, double *x, int leny, double *y);
		// USEFUL FUNCTIONS
		int computeLaguerre(int lenpsi, double * psi);
		double FindWopt(int lenw, double * w, int leny, double *y,double alpha);
		// int computeTesselation(){Lag->preprocessing();this->Tess=Lag->computeTesselation(); return Tess.nnzHessian;};
		int computeTesselation(){this->Tess=Lag->computeTesselation(); return Tess.nnzHessian;};
		std::list<std::list<double>> computeTess(){return Lag->computeTess();};
		std::list<std::list<int>> computeAdjacency(){return Lag->computeAdjacency();};
		void Integration(int n, double *Bar,int m, double *Cost,int k, double *Mass, int p, double *val, int q , long *col, int r, long *row);
		void IntegrationParall(int n, double *Bar,int m, double *Cost,int k, double *Mass, int p, double *val, int q , long *col, int r, long *row);


		void HessianCompute2(int p, double *val, int q , long *col, int r, long *row);

		void HessianComputeParall(int p, double *val, int q , long *col, int r, long *row);
		void ComputeLocalMatrix(int n, double *MomentX2,int m, double *MomentY2,int k, double *MomentXY);
};

#endif
