#include "integrateQ1.h"


/**
 * @brief      Constructs the object.
 *
 * @param      u00   \f$\mu \verb?[i,j]?\f$
 * @param      u10   \f$\mu \verb?[i+1,j]?-\mu \verb?[i,j]?\f$
 * @param      u01   \f$\mu \verb?[i,j+1]?-\mu \verb?[[i,j]?\f$
 * @param      u11   \f$\mu \verb?[i+1,j+1]? + \mu \verb?[i,j]? 
 *                    - \left( \mu \verb?[i+1,j]? \mu \verb?[i,j+1]? \right)\f$
 * @param[in]  Nx    { parameter_description }
 * @param[in]  Ny    { parameter_description }
 * @param[in]  Dx    { parameter_description }
 * @param[in]  Dy    { parameter_description }
 */
VolIntegrandQ1::VolIntegrandQ1(double * u00, double * u10, double * u01, double * u11,
						 int Nx, int Ny,double Dx, double Dy) : VolIntegrand(Nx,Ny,Dx,Dy){
	this->u00 = u00;
	this->u10 = u10;
	this->u01 = u01;
	this->u11 = u11;		
}


/**
 * @brief      Constructs the object.
 *
 * @param[in]  LefB  The left boundary
 * @param[in]  RigB  The right boundary
 * @param[in]  BotB  The bottom boundary
 * @param[in]  UppB  The upper boundary
 * @param[in]  NxImage  The nx image
 * @param[in]  NyImage  The ny image
 * @param[in]  lenw     The lenw
 * @param      w        The image
 */
IntegrateQ1::IntegrateQ1(double LefB,double RigB,double BotB,double UppB,
						 int NxImage,int NyImage, int lenw, double* w):
						IntegrationMethod(LefB,RigB,BotB,UppB){
	
	this->NxImage = NxImage;
	this->NyImage = NyImage;
	
	this->Nx= NxImage - 1;
	this->Ny= NyImage - 1;
	
	this->Dx = (this->RightBoundary - this->LeftBoundary  )/((double)this->Nx);	
	this->Dy = (this->UpperBoundary - this->BottomBoundary)/((double)this->Ny);	
	if (!(NxImage*NyImage==lenw)) std::cout<<"Raise ValueError a "<<NxImage<<" "<<NyImage<<" "<<lenw<<std::endl;
	this->image = (double *) calloc(this->NxImage*this->NyImage,sizeof(double));
	std::memcpy(this->image,w,this->NxImage*this->NyImage*sizeof(double));

	this->u00 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->u10 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->u01 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->u11 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	
	for (int j=0; j<this->Ny; j++){
		for(int i = 0; i < this->Nx; i++){
			this->u00[i*(this->Ny)+j] = img_ij(i,j);
			this->u10[i*(this->Ny)+j] = (img_ij(i+1,j)-img_ij(i,j));
			this->u01[i*(this->Ny)+j] = (img_ij(i,j+1)-img_ij(i,j));
			this->u11[i*(this->Ny)+j] = (img_ij(i+1,j+1)-img_ij(i+1,j)-img_ij(i,j+1)+img_ij(i,j));
		}
	}
	

	
	this->s0Field  = new s0 (this->u00, this->u10, this->u01, this->u11, this->Nx, this->Ny, this->Dx, this->Dy);	
	this->s1xField = new s1x(this->u00, this->u10, this->u01, this->u11, this->Nx, this->Ny, this->Dx, this->Dy);
	this->s1yField = new s1y(this->u00, this->u10, this->u01, this->u11, this->Nx, this->Ny, this->Dx, this->Dy);
	this->s2xField = new s2x(this->u00, this->u10, this->u01, this->u11, this->Nx, this->Ny, this->Dx, this->Dy);
	this->s2yField = new s2y(this->u00, this->u10, this->u01, this->u11, this->Nx, this->Ny, this->Dx, this->Dy);
	this->sxyField = new sxy(this->u00, this->u10, this->u01, this->u11, this->Nx, this->Ny, this->Dx, this->Dy);
	
}

/**
 * @brief      Compute the volumic integrate of the density
 *
 * @param      int_seg  The int segment
 *
 * @return      \f$ \int_{int_seg} \vec{s} \cdot \vec{n} ?\f$ where \f$ \nabla \cdot \vec{s} = \mu $\f
 */
double IntegrateQ1::IntMass(Integration_segment * int_seg){
	return (*int_seg).diffTime*(*int_seg).direction[1]*IntegrationMethod::gauss2pt(int_seg,static_cast<VolIntegrand*>(this->s0Field));
}

/**
 * @brief      Compute the volumic integrate of the density
 *
 * @param      int_seg  The int segment
 *
 * @return      \f$ \int_{int_seg} \vec{s} \cdot \vec{n} ?\f$ where \f$ \nabla \cdot \vec{s} = \mu x $\f
 */
double IntegrateQ1::IntX1(Integration_segment * int_seg){
	return (*int_seg).diffTime*(*int_seg).direction[1]*IntegrationMethod::gauss3pt(int_seg,static_cast<VolIntegrand*>(this->s1xField));
}
/**
 * @brief      Compute the volumic integrate of the density
 *
 * @param      int_seg  The int segment
 *
 * @return      \f$ \int_{int_seg} \vec{s} \cdot \vec{n} ?\f$ where \f$ \nabla \cdot \vec{s} = \mu y $\f
 */
double IntegrateQ1::IntY1(Integration_segment * int_seg){
	return -(*int_seg).diffTime*(*int_seg).direction[0]*IntegrationMethod::gauss3pt(int_seg,static_cast<VolIntegrand*>(this->s1yField));
}
/**
 * @brief      Compute the volumic integrate of the density
 *
 * @param      int_seg  The int segment
 *
 * @return      \f$ \int_{int_seg} \vec{s} \cdot \vec{n} ?\f$ where \f$ \nabla \cdot \vec{s} = \mu x^2$\f
 */
double IntegrateQ1::IntX2(Integration_segment * int_seg){
	return (*int_seg).diffTime*(*int_seg).direction[1]*IntegrationMethod::gauss3pt(int_seg,static_cast<VolIntegrand*>(this->s2xField));
}
/**
 * @brief      Compute the volumic integrate of the density
 *
 * @param      int_seg  The int segment
 *
 * @return      \f$ \int_{int_seg} \vec{s} \cdot \vec{n} ?\f$ where \f$ \nabla \cdot \vec{s} = \mu  y^2$\f
 */
double IntegrateQ1::IntY2(Integration_segment * int_seg){
	return -(*int_seg).diffTime*(*int_seg).direction[0]*IntegrationMethod::gauss3pt(int_seg,static_cast<VolIntegrand*>(this->s2yField));
}

double IntegrateQ1::IntXY(Integration_segment * int_seg){
	return (*int_seg).diffTime*(*int_seg).direction[1]*IntegrationMethod::gauss3pt(int_seg,static_cast<VolIntegrand*>(this->sxyField));
}
/**
 * @brief      Compute the sufarce integrate of the density
 *
 * @param      int_seg  The int segment
 *
 * @return     All the integration needed for the Hessian w.to \f$ psi $\f
 */
void IntegrateQ1::intHess(Integration_segment * int_seg, double & HessiiWW, double & HessijWW){
	double theta12, theta22, beta12, beta22;
	theta12 = (1 - gauss2Point1)*int_seg->thetaStartX + gauss2Point1*int_seg->thetaEndX;
	beta12 =  (1 - gauss2Point1)*int_seg->thetaStartY + gauss2Point1*int_seg->thetaEndY;
	theta22 = (1 - gauss2Point2)*int_seg->thetaStartX + gauss2Point2*int_seg->thetaEndX;
	beta22 =  (1 - gauss2Point2)*int_seg->thetaStartY + gauss2Point2*int_seg->thetaEndY;
	int i = int_seg->indi;
	int j = int_seg->indj;
	double f01,f02;
	f01 = u00[i*Ny+j]+theta12*u10[i*Ny+j]+beta12*u01[i*Ny+j]+theta12*beta12*u11[i*Ny+j];
	f02 = u00[i*Ny+j]+theta22*u10[i*Ny+j]+beta22*u01[i*Ny+j]+theta22*beta22*u11[i*Ny+j];
	double f0 = 0.5*(f01+f02);
	HessiiWW += f0*int_seg->diffTime*int_seg->normDir/(int_seg->twoNormInOut);
	HessijWW -= f0*int_seg->diffTime*int_seg->normDir/(int_seg->twoNormInOut);
}


void IntegrateQ1::intHess2(Integration_segment * int_seg,
                         double & HessiiWX, double & HessijWX, double & HessiiWY, double & HessijWY,
                         double & HessiiXX, double & HessijXX, double & HessiiXY, double & HessijXY,
						 double & HessijYX,
                         double & HessiiYY, double & HessijYY){
    double theta12, theta22, beta12, beta22;
    theta12 = (1 - gauss2Point1)*int_seg->thetaStartX + gauss2Point1*int_seg->thetaEndX;
    beta12 = (1 - gauss2Point1)*int_seg->thetaStartY + gauss2Point1*int_seg->thetaEndY;
    theta22 = (1 - gauss2Point2)*int_seg->thetaStartX + gauss2Point2*int_seg->thetaEndX;
    beta22 = (1 - gauss2Point2)*int_seg->thetaStartY + gauss2Point2*int_seg->thetaEndY;
    int i = int_seg->indi;
    int j = int_seg->indj;
    double f01,f02;
    f01 = u00[i*Ny+j]+theta12*u10[i*Ny+j]+beta12*u01[i*Ny+j]+theta12*beta12*u11[i*Ny+j];
    f02 = u00[i*Ny+j]+theta22*u10[i*Ny+j]+beta22*u01[i*Ny+j]+theta22*beta22*u11[i*Ny+j];

    double f1x, f1y;
    f1x = 0.5*Dx*((theta12+i)*f01 + (theta22+i)*f02);
    f1y = 0.5*Dy*((beta12+j)*f01 + (beta22+j)*f02);
    double f0 = 0.5*(f01+f02);
    HessiiWX += (f1x-int_seg->xPosInside*f0)*int_seg->diffTime*int_seg->normDir/(int_seg->normInOut);
    HessijWX -= (f1x-int_seg->xPosOutside*f0)*int_seg->diffTime*int_seg->normDir/(int_seg->normInOut);
    HessiiWY += (f1y-int_seg->yPosInside*f0)*int_seg->diffTime*int_seg->normDir/(int_seg->normInOut);
    HessijWY -= (f1y-int_seg->yPosOutside*f0)*int_seg->diffTime*int_seg->normDir/(int_seg->normInOut);

    double theta13, theta23, theta33, beta13, beta23, beta33;
    theta13 = (1 - gauss3Point1)*int_seg->thetaStartX + gauss3Point1*int_seg->thetaEndX;
    beta13  = (1 - gauss3Point1)*int_seg->thetaStartY + gauss3Point1*int_seg->thetaEndY;
    theta23 = (1 - gauss3Point2)*int_seg->thetaStartX + gauss3Point2*int_seg->thetaEndX;
    beta23  = (1 - gauss3Point2)*int_seg->thetaStartY + gauss3Point2*int_seg->thetaEndY;
    theta33 = (1 - gauss3Point3)*int_seg->thetaStartX + gauss3Point3*int_seg->thetaEndX;
    beta33  = (1 - gauss3Point3)*int_seg->thetaStartY + gauss3Point3*int_seg->thetaEndY;

    double f001,f002,f003;
    f001 = u00[i*Ny+j]+theta13*u10[i*Ny+j]+beta13*u01[i*Ny+j]+theta13*beta13*u11[i*Ny+j];
    f002 = u00[i*Ny+j]+theta23*u10[i*Ny+j]+beta23*u01[i*Ny+j]+theta23*beta23*u11[i*Ny+j];
    f003 = u00[i*Ny+j]+theta33*u10[i*Ny+j]+beta33*u01[i*Ny+j]+theta33*beta33*u11[i*Ny+j];
    
    double f2x,f2y,f1x1y;
    f2x = gauss3Factor1*std::pow(Dx*(theta13+i),2.0)*f001
        + gauss3Factor2*std::pow(Dx*(theta23+i),2.0)*f002
        + gauss3Factor3*std::pow(Dx*(theta33+i),2.0)*f003;
    f2y = gauss3Factor1*std::pow(Dy*(beta13+j),2.0)*f001
        + gauss3Factor2*std::pow(Dy*(beta23+j),2.0)*f002
        + gauss3Factor3*std::pow(Dy*(beta33+j),2.0)*f003;
    f1x1y = gauss3Factor1*(Dx*(theta13+i))*(Dy*(beta13+j))*f001
         + gauss3Factor2*(Dx*(theta23+i))*(Dy*(beta23+j))*f002
         + gauss3Factor3*(Dx*(theta33+i))*(Dy*(beta33+j))*f003;
    
    HessiiXX += (-f2x+2.0*f1x*int_seg->xPosInside-f0*std::pow(int_seg->xPosInside,2.0))*2.0*int_seg->diffTime*int_seg->normDir/(int_seg->normInOut);
    HessijXX -= (-f2x+(int_seg->xPosInside+int_seg->xPosOutside)*f1x-int_seg->xPosInside*int_seg->xPosOutside*f0)*2.0*int_seg->diffTime*int_seg->normDir/(int_seg->normInOut);
    HessiiXY += (2.0/(int_seg->normInOut))*int_seg->diffTime*int_seg->normDir*(-f1x1y+int_seg->xPosInside*f1y+int_seg->yPosInside*f1x-int_seg->xPosInside*int_seg->yPosInside*f0);
    HessijXY -= (2.0/(int_seg->normInOut))*int_seg->diffTime*int_seg->normDir*(-f1x1y+int_seg->xPosInside*f1y+int_seg->yPosOutside*f1x-int_seg->xPosInside*int_seg->yPosOutside*f0);
    HessijYX -= (2.0/(int_seg->normInOut))*int_seg->diffTime*int_seg->normDir*(-f1x1y+int_seg->yPosInside*f1x+int_seg->xPosOutside*f1y-int_seg->yPosInside*int_seg->xPosOutside*f0);
    HessiiYY += (2.0/(int_seg->normInOut))*int_seg->diffTime*int_seg->normDir*(-f2y+2.0*f1y*int_seg->yPosInside-f0*std::pow(int_seg->yPosInside,2.0));
    HessijYY -= (2.0/(int_seg->normInOut))*int_seg->diffTime*int_seg->normDir*(-f2y+(int_seg->yPosInside+int_seg->yPosOutside)*f1y-int_seg->yPosInside*int_seg->yPosOutside*f0);
}



s0::s0(double * u00, double * u10, double * u01, double * u11, int Nx, int Ny, double Dx, double Dy) : 
VolIntegrandQ1(u00, u10, u01, u11, Nx, Ny, Dx,Dy) {
	this->r0 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->r1 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	for (int j=0; j<Ny; j++){
		for(int i = 0; i < Nx-1; i++){
			this->r0[(i+1)*Ny+j] = this->r0[i*Ny+j] + 0.5*u10[i*Ny + j]+u00[i*Ny + j];
			this->r1[(i+1)*Ny+j] = this->r1[i*Ny+j] + 0.5*u11[i*Ny + j]+u01[i*Ny + j];			
		}
	}	
}
/**
 * @brief      Evalutate the vector
 *             @f$ \vec{s}$\f for the computation of <CODE>IntMass<\CODE>
 *
 * @param      int_seg  The int segment
 * @param[in]  t        { parameter_description }
 * @param      u00   The u 00
 *
 * @return     { description_of_the_return_value }
 */
double s0::evaluate(Integration_segment * int_seg,double t){
	int i = int_seg->indi;
	int j = int_seg->indj;
	double theta, beta;
	theta = (1.0-t)*int_seg->thetaStartX + t*int_seg->thetaEndX;
	beta  = (1.0-t)*int_seg->thetaStartY + t*int_seg->thetaEndY;
	double halfThetaSquare = std::pow(theta,2.0)*0.5;
	return this->Dx*( theta*u00[i*Ny + j] + (halfThetaSquare)*(u10[i*Ny + j]) 
					  + theta*beta*(u01[i*Ny + j]) + (beta*halfThetaSquare)
					  *(u11[i*Ny + j]) + this->r0[i*Ny + j] + this->r1[i*Ny + j]*beta);
}

s1x::s1x(double * u00, double * u10, double * u01, double * u11,int Nx, int Ny, double Dx, double Dy) : VolIntegrandQ1(u00, u10, u01, u11, Nx, Ny, Dx,Dy){
	this->r0 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->r1 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->dx2 = std::pow(Dx,2);
	for (int j=0; j<Ny; j++){
		for(int i = 0; i < Nx-1; i++){		
			this->r0[(i+1)*Ny+j] = this->r0[(i)*Ny+j] + i*(this->u00[i*Ny+j]+0.5*this->u10[i*Ny+j]) 
								 + 0.5*this->u00[i*Ny+j]+oneThird*this->u10[i*Ny+j];
			this->r1[(i+1)*Ny+j] = this->r1[i*Ny+j] + i*(this->u01[i*Ny+j]+0.5*this->u11[i*Ny+j]) 
								 + 0.5*this->u01[i*Ny+j]+oneThird*this->u11[i*Ny+j];			
		}
	}			
}
/**
 * @brief      Evalutate the vector
 *             @f$ \vec{s}$\f for the computation of <CODE>IntX1<\CODE>
 *
 * @param      int_seg  The int segment
 * @param[in]  t        { parameter_description }
 * @param      u00   The u 00
 *
 * @return     { description_of_the_return_value }
 */
double s1x::evaluate(Integration_segment * int_seg,double t) {
	int i = int_seg->indi;
	int j = int_seg->indj;
	double theta, beta;
	theta = (1.0-t)*int_seg->thetaStartX + t*int_seg->thetaEndX;
	beta  = (1.0-t)*int_seg->thetaStartY + t*int_seg->thetaEndY;
	double thetaSquare = std::pow(theta,2);	
	double thetaCube = std::pow(theta,3); 
	return dx2*(i*(theta*u00[i*Ny + j] + thetaSquare*0.5*u10[i*Ny + j] + beta*theta*u01[i*Ny + j] + beta*thetaSquare*.5*u11[i*Ny + j])
				+ (thetaSquare*.5)*u00[i*Ny + j] + thetaCube*oneThird*u10[i*Ny + j] + thetaSquare*beta*0.5*u01[i*Ny + j] 
				+ thetaCube*beta*oneThird*u11[i*Ny + j]
				+ this->r0[i*Ny + j]+this->r1[i*Ny + j]*beta);
}

s1y::s1y(double * u00, double * u10, double * u01, double * u11,int Nx, int Ny, double Dx, double Dy) : VolIntegrandQ1(u00, u10, u01, u11, Nx, Ny, Dx,Dy){	
	this->r0 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->r1 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->dy2 = std::pow(this->Dy,2);
	for(int i = 0; i < this->Nx; i++){
		for (int j=0; j<this->Ny-1; j++){			
			this->r0[i*Ny+j+1] = this->r0[i*Ny+j] + j*(this->u00[i*Ny+j]+0.5*this->u01[i*Ny+j]) 
											   + (0.5*this->u00[i*Ny+j]+oneThird*this->u01[i*Ny+j]);
			this->r1[i*Ny+j+1] = this->r1[i*Ny+j] + j*(this->u10[i*Ny+j]+0.5*this->u11[i*Ny+j]) 
											   + (0.5*this->u10[i*Ny+j]+oneThird*this->u11[i*Ny+j]);			
		}
	}
}
/**
 * @brief      Evalutate the vector
 *             @f$ \vec{s}$\f for the computation of <CODE>IntY1<\CODE>
 *
 * @param      int_seg  The int segment
 * @param[in]  t        { parameter_description }
 * @param      u00   The u 00
 *
 * @return     { description_of_the_return_value }
 */
double s1y::evaluate(Integration_segment * int_seg,double t) {
	int i = int_seg->indi;
	int j = int_seg->indj;
	double theta, beta;
	theta = (1.0-t)*int_seg->thetaStartX + t*int_seg->thetaEndX;
	beta  = (1.0-t)*int_seg->thetaStartY + t*int_seg->thetaEndY;
	double halfBetaSquare = 0.5*std::pow(beta,2);
	double betaCube = std::pow(beta,3);
	return dy2*( j*( beta*u00[i*Ny + j] + u10[i*Ny + j]*theta*beta +u01[i*Ny + j]*halfBetaSquare 
				+ u11[i*Ny + j]*theta*halfBetaSquare ) + ( u00[i*Ny + j]*halfBetaSquare 
				+ u10[i*Ny + j]*theta*halfBetaSquare + u01[i*Ny + j]*betaCube*oneThird +
				 u11[i*Ny + j]*theta*betaCube*oneThird )+ (this->r0[i*Ny + j]+this->r1[i*Ny + j]*theta));
}



s2x::s2x(double * u00, double * u10, double * u01, double * u11, int Nx, int Ny, double Dx, double Dy) : VolIntegrandQ1(u00, u10, u01, u11, Nx, Ny, Dx,Dy){
	this->r0 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->r1 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->dx3 = std::pow(Dx,3.0);
	for (int j=0; j<Ny; j++){
		for(int i = 0; i < Nx-1; i++){			
			this->r0[(i+1)*Ny+j] = this->r0[i*Ny+j] + std::pow(i,2)*(this->u00[i*Ny+j]+0.5*this->u10[i*Ny+j]) 
											   + 2*i*(0.5*this->u00[i*Ny+j]+oneThird*this->u10[i*Ny+j]) 
											   + (oneThird*this->u00[i*Ny+j]+0.25*this->u10[i*Ny+j]);
			this->r1[(i+1)*Ny+j] = this->r1[i*Ny+j] + std::pow(i,2)*(this->u01[i*Ny+j]+0.5*this->u11[i*Ny+j]) 
											   + 2*i*(0.5*this->u01[i*Ny+j]+oneThird*this->u11[i*Ny+j]) 
											   + (oneThird*this->u01[i*Ny+j]+0.25*this->u11[i*Ny+j]);
		}
	}
}
/**
 * @brief      Evalutate the vector
 *             @f$ \vec{s}$\f for the computation of <CODE>IntX2<\CODE>
 *
 * @param      int_seg  The int segment
 * @param[in]  t        { parameter_description }
 * @param      u00   The u 00
 *
 * @return     { description_of_the_return_value }
 */
double s2x::evaluate(Integration_segment * int_seg,double t){
	int i = int_seg->indi;
	int j = int_seg->indj;
	int ind=i*Ny+j;
	double theta_1 = (1.0-t)*int_seg->thetaStartX + t*int_seg->thetaEndX;
	double beta    = (1.0-t)*int_seg->thetaStartY + t*int_seg->thetaEndY;
	double theta_2 = 0.5*std::pow(theta_1,2);
	double theta_3 = std::pow(theta_1,3)*oneThird;
	double theta_4 = std::pow(theta_1,4)/4.0;

	return dx3*( std::pow(i,2)*(theta_1*u00[ind] + theta_2*u10[ind] + theta_1*beta*u01[ind] + theta_2*beta*u11[ind])
					   + 2.0*i*(theta_2*u00[ind] + theta_3*u10[ind] + theta_2*beta*u01[ind] + theta_3*beta*u11[ind])
		   					 + (theta_3*u00[ind] + theta_4*u10[ind] + theta_3*beta*u01[ind] + theta_4*beta*u11[ind])
				 + this->r0[ind]+beta*this->r1[ind] );
}

s2y::s2y(double * u00, double * u10, double * u01, double * u11, int Nx, int Ny, double Dx, double Dy) : VolIntegrandQ1(u00, u10, u01, u11, Nx, Ny, Dx,Dy){
	this->r0 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->r1 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->dy3 = std::pow(Dx,3.0);
	for(int i = 0; i < Nx; i++){
		for (int j=0; j<Ny-1; j++){			
			
			this->r0[i*Ny+j+1] = this->r0[i*Ny+j] + std::pow(j,2)*(this->u00[i*Ny+j]+0.5*this->u01[i*Ny+j]) 
											   + 2*j*(0.5*this->u00[i*Ny+j]+oneThird*this->u01[i*Ny+j])
											   + (oneThird*this->u00[i*Ny+j]+0.25*this->u01[i*Ny+j]);
			this->r1[i*Ny+j+1] = this->r1[i*Ny+j] + std::pow(j,2)*(this->u10[i*Ny+j]+0.5*this->u11[i*Ny+j]) 
											   + 2*j*(0.5*this->u10[i*Ny+j]+oneThird*this->u11[i*Ny+j]) 
											   + (oneThird*this->u10[i*Ny+j]+0.25*this->u11[i*Ny+j]);			
		}
	}
}
/**
 * @brief      Evalutate the vector
 *             @f$ \vec{s}$\f for the computation of <CODE>IntY2<\CODE>
 *
 * @param      int_seg  The int segment
 * @param[in]  t        { parameter_description }
 * @param      u00   The u 00
 *
 * @return     { description_of_the_return_value }
 */
double s2y::evaluate(Integration_segment * int_seg,double t){
	int i = int_seg->indi;
	int j = int_seg->indj;
	int ind=i*Ny+j;
	double theta  = (1.0-t)*int_seg->thetaStartX + t*int_seg->thetaEndX;
	double beta_1 = (1.0-t)*int_seg->thetaStartY + t*int_seg->thetaEndY;
	double beta_2 = std::pow(beta_1,2)*.5;
	double beta_3 = std::pow(beta_1,3)*oneThird;
	double beta_4 = std::pow(beta_1,4)*.25;
	return dy3*( std::pow(j,2)*( beta_1*u00[ind] + theta*beta_1*u10[ind] + beta_2*u01[ind] + theta*beta_2*u11[ind])
				+          2*j*( beta_2*u00[ind] + theta*beta_2*u10[ind] + beta_3*u01[ind] + theta*beta_3*u11[ind])
				+              ( beta_3*u00[ind] + theta*beta_3*u10[ind] + beta_4*u01[ind] + theta*beta_4*u11[ind])
				+ this->r0[ind]+theta*this->r1[ind] );
}

sxy::sxy(double * u00, double * u10, double * u01, double * u11, int Nx, int Ny, double Dx, double Dy) : VolIntegrandQ1(u00, u10, u01, u11, Nx, Ny, Dx,Dy){
	this->r0 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->r1 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->r2 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->dx2dy = std::pow(Dx,2)*Dy;
	for (int j=0; j<Ny; j++){
		for(int i = 0; i < Nx-1; i++){
			int ind =i*Ny+j;
			this->r0[(i+1)*Ny+j] = this->r0[ind] + i*j*(this->u00[ind]+     0.5*this->u10[ind]) 
											   + j*(0.5*this->u00[ind]+oneThird*this->u10[ind]);
			this->r1[(i+1)*Ny+j] = this->r1[ind] +   i*(this->u00[ind]+     0.5*this->u10[ind]) 
											   +   (0.5*this->u00[ind]+oneThird*this->u10[ind])
											      +i*j*(this->u01[ind]+     0.5*this->u11[ind])
											    +j*(0.5*this->u01[ind]+oneThird*this->u11[ind]);
			this->r2[(i+1)*Ny+j] = this->r2[ind]   + i*(this->u01[ind]+     0.5*this->u11[ind]) 
											     + (0.5*this->u01[ind]+oneThird*this->u11[ind]);
		}
	}
}
/**
 * @brief      Evalutate the vector
 *             @f$ \vec{s}$\f for the computation of <CODE>IntX2<\CODE>
 *
 * @param      int_seg  The int segment
 * @param[in]  t        { parameter_description }
 * @param      u00   The u 00
 *
 * @return     { description_of_the_return_value }
 */
double sxy::evaluate(Integration_segment * int_seg,double t){
	int i = int_seg->indi;
	int j = int_seg->indj;
	int ind =i*Ny+j;
	double theta_1 = (1.0-t)*int_seg->thetaStartX + t*int_seg->thetaEndX;
	double beta    = (1.0-t)*int_seg->thetaStartY + t*int_seg->thetaEndY;
	double theta_2 = 0.5*std::pow(theta_1,2);
	double theta_3 = oneThird*std::pow(theta_1,3);
	return dx2dy*(
		i*(beta+j)*(theta_1*u00[ind] + theta_2*u10[ind] + beta*theta_1*u01[ind]+beta*theta_2*u11[ind])
		+ (beta+j)*(theta_2*u00[ind] + theta_3*u10[ind] + beta*theta_2*u01[ind]+beta*theta_3*u11[ind])
		+r0[ind]+beta*r1[ind]+std::pow(beta,2)*r2[ind]
				 );
}