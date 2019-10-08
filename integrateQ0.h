#ifndef INTEGRATEQ0_H
#define INTEGRATEQ0_H

#include "integrate.h"

/**
 * @brief      Class for volumic integration for constant by pixel images.
 */
class VolIntegrandQ0 : public VolIntegrand {
	protected:
	double * r0;
	double * u00;
	
	public:
	VolIntegrandQ0(double * u00, 
				int Nx, int Ny, double Dx, double Dy);
	VolIntegrandQ0(): VolIntegrand(){};			
	virtual ~VolIntegrandQ0(){}
	virtual double evaluate(Integration_segment * int_seg,double t) = 0;
};

/**
 * @brief      Order 0 moment volume integrand
 */
class Q0s0 : public VolIntegrandQ0{
	public:
	Q0s0(): VolIntegrandQ0(){};
	Q0s0(double * u00,int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t);
	~Q0s0(){
		delete r0;
	}
};

/**
 * @brief  Order 1 moment in the \f$\vec{e}_x\f$ direction volume integrand
 */
class Q0s1x : public VolIntegrandQ0{
	private:
	double dx2;
	const double oneThird = 1.0/3.0;
	public:
	Q0s1x(): VolIntegrandQ0(){};
	Q0s1x(double * u00, int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t);
	~Q0s1x(){
		delete r0;
	}
};

/**
 * @brief      Order 1 moment in the \f$\vec{e}_y\f$ direction volume integrand
 */
class Q0s1y : public VolIntegrandQ0{
	private:
	const double oneThird = 1.0/3.0;
	double dy2;
	public:
	Q0s1y(): VolIntegrandQ0(){};
	Q0s1y(double * u00, int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t);
	~Q0s1y(){
		delete r0;
	}
};

/**
 * @brief      Order 2 moment in the f$\vec{e}_x\f$ direction volume integrand
 */
class Q0s2x : public VolIntegrandQ0{
	private: 
	const double oneThird = 1.0/3.0;
	double dx3;
	public:
	Q0s2x(): VolIntegrandQ0(){};
	Q0s2x(double * u00, int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t);
	~Q0s2x(){
		delete r0;
	}
};

/**
 * @brief      Order 2 moment in the \f$\vec{e}_y\f$ direction volume integrand
 */
class Q0s2y : public VolIntegrandQ0{
	private: 
	const double oneThird = 1.0/3.0;
	double dy3;
	public:
	Q0s2y(): VolIntegrandQ0(){std::cout<<"DESTRUCTOR : "<<__LINE__<<" in "<<__FILE__<<std::endl;};
	Q0s2y(double * u00, int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t);
	~Q0s2y(){
		delete r0;
	}
};

/**
 * @brief      Order 2 moment in the \f$\vec{e}_y\f$ direction volume integrand
 */
class Q0sxy : public VolIntegrandQ0{
	private: 
	double dx2dy;
	double * r1;
	public:
	Q0sxy(): VolIntegrandQ0(){std::cout<<"DESTRUCTOR : "<<__LINE__<<" in "<<__FILE__<<std::endl;};
	Q0sxy(double * u00, int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t);
	~Q0sxy(){
		delete r0;
		delete r1;
	}
};


/**
 * @brief      Derived class of Integration method used if \f$\mu\f$ is a Q0 density, i.e. a constant by
 *             pixel density.
 */
class IntegrateQ0 : public IntegrationMethod{
	private:
	double * u00;
	public:
	Q0s0  * s0Field;
	Q0s1x * s1xField;
	Q0s1y * s1yField;
	Q0s2x * s2xField;
	Q0s2y * s2yField;
	Q0sxy * sxyField;

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
	IntegrateQ0() : IntegrationMethod(){};
	IntegrateQ0(double LeftB,double RightB,double BottomB,double UpperB,int NxImage,int NyImage,
				int lenw, double* w);
	~IntegrateQ0(){
		delete s0Field;
		delete s1xField;
		delete s2xField;
		delete s1yField;
		delete s2yField;
		delete sxyField;
		delete u00;
	}

};


#endif