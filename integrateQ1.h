#ifndef INTEGRATEQ1_H
#define INTEGRATEQ1_H

#include "integrate.h"
/**
 * @brief      Class for volumic integration for bilinear by pixel images.
 */
class VolIntegrandQ1 : public VolIntegrand {
	protected:
	double * r0;
	double * r1;
	double * u00;
	double * u10;
	double * u01;
	double * u11;
	
	public:
	VolIntegrandQ1(double * u00, double * u10, double * u01, double * u11, 
				int Nx, int Ny, double Dx, double Dy);
	VolIntegrandQ1(): VolIntegrand(){};			
	virtual ~VolIntegrandQ1(){}
	virtual double evaluate(Integration_segment * int_seg,double t) = 0;
};
/**
 * @brief      Order 0 moment volume integrand
 */
class s0 : public VolIntegrandQ1{
	public:
	s0(): VolIntegrandQ1(){};
	s0(double * u00, double * u10, double * u01, double * u11,int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t);
	~s0(){
		delete r0;
		delete r1;
		}
};

/**
 * @brief  Order 1 moment in the \f$\vec{e}_x\f$ direction volume integrand
 */
class s1x : public VolIntegrandQ1{
	private:
	double dx2;
	const double oneThird = 1.0/3.0;
	public:
	s1x(): VolIntegrandQ1(){};
	s1x(double * u00, double * u10, double * u01, double * u11, int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t);
	~s1x(){
	delete r0;
	delete r1;
	}
};
/**
 * @brief      Order 1 moment in the \f$\vec{e}_y\f$ direction volume integrand
 */
class s1y : public VolIntegrandQ1{
	private:
	const double oneThird = 1.0/3.0;
	double dy2;
	public:
	s1y(): VolIntegrandQ1(){};
	s1y(double * u00, double * u10, double * u01, double * u11, int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t);
	~s1y(){
	delete r0;
	delete r1;
	}
};
/**
 * @brief      Order 2 moment in the f$\vec{e}_x\f$ direction volume integrand
 */
class s2x : public VolIntegrandQ1{
	private: 
	const double oneThird = 1.0/3.0;
	double dx3;
	public:
	s2x(): VolIntegrandQ1(){};
	s2x(double * u00, double * u10, double * u01, double * u11, int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t);
	~s2x(){
	delete r0;
	delete r1;
	}
};
/**
 * @brief      Order 2 moment in the \f$\vec{e}_y\f$ direction volume integrand
 */
class s2y : public VolIntegrandQ1{
	private: 
	const double oneThird = 1.0/3.0;
	double dy3;
	public:
	s2y(): VolIntegrandQ1(){};
	s2y(double * u00, double * u10, double * u01, double * u11, int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t);
	~s2y(){
	delete r0;
	delete r1;
	}
};
class sxy : public VolIntegrandQ1{
	private: 
	const double oneThird = 1.0/3.0;
	double dx2dy;
	double * r2;
	public:
	sxy(): VolIntegrandQ1(){};
	sxy(double * u00, double * u10, double * u01, double * u11, int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t);
	~sxy(){
	delete r0;
	delete r1;
	delete r2;
	}
};

/**
 * @brief      Derived class of Integration method used if \f$\mu\f$ is a Q1 density, i.e. a bilinear by
 *             pixel density.
 */
class IntegrateQ1 : public IntegrationMethod{
	private:
	double * u00,* u10,* u01,* u11;
	public:
	s0 * s0Field;
	s1x * s1xField;
	s1y * s1yField;
	s2x * s2xField;
	s2y * s2yField;
	sxy * sxyField;
	virtual double IntMass(Integration_segment * int_seg);
	virtual double IntX1(Integration_segment * int_seg);
	virtual double IntY1(Integration_segment * int_seg);
	virtual double IntX2(Integration_segment * int_seg);
	virtual double IntY2(Integration_segment * int_seg);
	virtual double IntXY(Integration_segment * int_seg);

	virtual void intHess(Integration_segment * int_seg, double & HessiiWW, double & HessijWW);
	virtual void intHess2(Integration_segment * int_seg, 
						 double & HessiiWX, double & HessijWX, double & HessiiWY, double & HessijWY,
						 double & HessiiXX, double & HessijXX, double & HessiiXY, double & HessijXY,double & HessijYX,
						 double & HessiiYY, double & HessijYY);
	IntegrateQ1() : IntegrationMethod(){};
	IntegrateQ1(double LeftB,double RightB,double BottomB,double UpperB,int NxImage,int NyImage,
				int lenw, double* w);
	~IntegrateQ1(){
		delete s0Field;
		delete s1xField;
		delete s2xField;
		delete s1yField;
		delete s2yField;
		delete sxyField;
		delete u00;delete u10;delete u01;delete u11;
	}

};


#endif