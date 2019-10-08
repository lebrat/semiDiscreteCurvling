#ifndef DENSITYQ1_H
#define DENSITYQ1_H

#include "primitives.h"
#include "edges.h"
#include "density.h"


class IntegrateQ1 : public IntegrationMethod{
	private:
	double * u00,* u10,* u01,* u11;
	const double LeftBoundary, RightBoundary,BottomBoundary,UpperBoundary;
	public:
	s0 * s0Field;
	s1x * s1xField;
	s1y * s1yField;
	s2x * s2xField;
	s2y * s2yField;
	virtual double IntMass(Integration_segment * int_seg);
	virtual double IntX1(Integration_segment * int_seg);
	virtual double IntY1(Integration_segment * int_seg);
	virtual double IntX2(Integration_segment * int_seg);
	virtual double IntY2(Integration_segment * int_seg);
	IntegrateQ1() : IntegrationMethod(), LeftBoundary(.0), RightBoundary(.0),
					BottomBoundary(.0), UpperBoundary(.0){};
	IntegrateQ1(double LeftB,double RightB,double BottomB,double UpperB,int NxImage,int NyImage,
				int lenw, double* w);
	~IntegrateQ1(){
		std::cout<<"le destructeur de IntegrateQ1 est appellé"<<std::endl;
		delete s0; delete s1x; delete s1y; delete s2x; delete s2y;
	}

};


class VectorFieldQ1 : public VectorField {
	protected:
	double * r0;
	double * r1;
	double * u00;
	double * u10;
	double * u01;
	double * u11;
	
	public:
	VectorFieldQ1(double * u00, double * u10, double * u01, double * u11, 
				int Nx, int Ny, double Dx, double Dy);
	VectorFieldQ1(): VectorField(){};			
	virtual ~VectorFieldQ1()
	{
		std::cout<<"DESTRUCTOR : "<<__LINE__<<" in "<<__FILE__<<std::endl;
		delete [] r0;
		delete [] r1;
	}
	virtual double evaluate(Integration_segment * int_seg,double t) = 0;
};
class s0 : public VectorFieldQ1{
	public:
	s0(): VectorFieldQ1(){};
	s0(double * u00, double * u10, double * u01, double * u11,int Nx, int Ny, double Dx, double Dy);
	void printNx(){std::cout<<Nx<<std::endl;};
	virtual double evaluate(Integration_segment * int_seg,double t);
	~s0(){std::cout<<"DESTRUCTOR : "<<__LINE__<<" in "<<__FILE__<<std::endl;}
};
class s1x : public VectorFieldQ1{
	private:
	double dx2;
	const double oneThird = 1.0/3.0;
	public:
	s1x(): VectorFieldQ1(){};
	s1x(double * u00, double * u10, double * u01, double * u11, int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t);
	~s1x(){std::cout<<"DESTRUCTOR : "<<__LINE__<<" in "<<__FILE__<<std::endl;}
};
class s1y : public VectorFieldQ1{
	private:
	const double oneThird = 1.0/3.0;
	double dy2;
	public:
	s1y(): VectorFieldQ1(){};
	s1y(double * u00, double * u10, double * u01, double * u11, int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t);
	~s1y(){std::cout<<"DESTRUCTOR : "<<__LINE__<<" in "<<__FILE__<<std::endl;}
};
class s2x : public VectorFieldQ1{
	private: 
	const double oneThird = 1.0/3.0;
	double dx3;
	public:
	s2x(): VectorFieldQ1(){};
	s2x(double * u00, double * u10, double * u01, double * u11, int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t);
	~s2x(){std::cout<<"DESTRUCTOR : "<<__LINE__<<" in "<<__FILE__<<std::endl;}
};
class s2y : public VectorFieldQ1{
	private: 
	const double oneThird = 1.0/3.0;
	double dy3;
	public:
	s2y(): VectorFieldQ1(){std::cout<<"DESTRUCTOR : "<<__LINE__<<" in "<<__FILE__<<std::endl;};
	s2y(double * u00, double * u10, double * u01, double * u11, int Nx, int Ny, double Dx, double Dy);
	virtual double evaluate(Integration_segment * int_seg,double t);
	~s2y(){}
};