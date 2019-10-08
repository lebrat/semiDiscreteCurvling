#ifndef INTEGRATE_H
#define INTEGRATE_H


#include "edges.h"


/**
 * @brief      { struct_description }
 */
struct Integration_segment{
	double diffTime;
	double cutPosStartX, cutPosStartY, cutPosEndX, cutPosEndY;
	double normDir,twoNormInOut,xPosInside,yPosInside,xPosOutside,yPosOutside;
	double normInOut;
	int indi, indj;
	double thetaEndX, thetaEndY, thetaStartX, thetaStartY, thetaX, thetaY;
	double direction [2];
};

/**
 * @brief      Class for volume integrand.
 */
class VolIntegrand{
	protected:
	int Nx,Ny;
	double Dx,Dy;
	public:
	VolIntegrand(){};
	virtual ~VolIntegrand(){};
	VolIntegrand(int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t) = 0;
};


/**
 * @brief      Virtual class for integration of <CODE>Integration_segment<\CODE>.
 */
class IntegrationMethod {
	protected:
	double gauss1Point,gauss2Point1,gauss2Point2,gauss3Point1,gauss3Point2,gauss3Point3;
	double gauss1Factor,gauss2Factor1,gauss2Factor2,gauss3Factor1,gauss3Factor2,gauss3Factor3;	
	double oneHalf,oneThird,oneFourth;
	int Nx, Ny, NxImage, NyImage;
	double Dx, Dy;
	const double LeftBoundary,RightBoundary,BottomBoundary,UpperBoundary;
	double * image, * Xpos,* Ypos;
	public:	
	IntegrationMethod() : LeftBoundary(0.0), RightBoundary(1.0), BottomBoundary(0.0), UpperBoundary(1.0){};
	IntegrationMethod(double l,double r,double b,double u);
	virtual ~IntegrationMethod(){
		delete [] image;
	}
	double dx(){return this->Dx;}
	double dy(){return this->Dy;}
	double gauss1pt(Integration_segment * int_seg,VolIntegrand * v);
	double gauss2pt(Integration_segment * int_seg,VolIntegrand * v);
	double gauss3pt(Integration_segment * int_seg,VolIntegrand * v);
	double img_ij(int i , int j);
	void precompute(Integration_segment * int_seg,Segmented_dual_edge &se, double current,double previous);
	std::list<double> sample(Segmented_dual_edge & e,double Tolx,double Toly);
	virtual double IntMass(Integration_segment * int_seg) = 0;
	virtual double IntX1(Integration_segment * int_seg) = 0;
	virtual double IntY1(Integration_segment * int_seg) = 0;
	virtual double IntX2(Integration_segment * int_seg) = 0;
	virtual double IntY2(Integration_segment * int_seg) = 0;
	virtual double IntXY(Integration_segment * int_seg) = 0;
	virtual void intHess(Integration_segment * int_seg, double & HessiiWW, double & HessijWW) = 0;
	virtual void intHess2(Integration_segment * int_seg, 
						 double & HessiiWX, double & HessijWX, double & HessiiWY, double & HessijWY,
						 double & HessiiXX, double & HessijXX, double & HessiiXY, double & HessijXY,double & HessijYX,
						 double & HessiiYY, double & HessijYY) = 0;


};

#endif
