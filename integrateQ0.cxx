#include "integrateQ0.h"


/**
 * @brief      Constructs the object.
 *
 * @param      \f$\mu \verb?[[i,j]?\f$
 * @param[in]  Nx    the Number of base Pixel/Rectangle in the \f$e_x\f$ direction
 * @param[in]  Ny    the Number of base Pixel/Rectangle in the \f$e_y\f$ direction
 * @param[in]  Dx    the distance between two consecutive neighbour in the \f$e_x\f$ direction
 * @param[in]  Dy    the distance between two consecutive neighbour in the \f$e_y\f$ direction
 */
VolIntegrandQ0::VolIntegrandQ0(double * u00,
						 int Nx, int Ny,double Dx, double Dy) : VolIntegrand(Nx,Ny,Dx,Dy){
	this->u00 = u00;		
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
IntegrateQ0::IntegrateQ0(double LefB,double RigB,double BotB,double UppB,
						 int NxImage,int NyImage, int lenw, double* w):
						IntegrationMethod(LefB,RigB,BotB,UppB){
	this->NxImage = NxImage;
	this->NyImage = NyImage;
	
	this->Nx= NxImage;
	this->Ny= NyImage;
	
	this->Dx = (this->RightBoundary - this->LeftBoundary)/((double)this->Nx);	
	this->Dy = (this->UpperBoundary - this->BottomBoundary)/((double)this->Ny);	
	if (!(NxImage*NyImage==lenw)) std::cout<<"Raise ValueError a "<<NxImage<<" "<<NyImage<<" "<<lenw<<std::endl;
	this->image = (double *) calloc(this->NxImage*this->NyImage,sizeof(double));
	std::memcpy(this->image,w,this->NxImage*this->NyImage*sizeof(double));
	this->u00 = (double *) calloc(this->NxImage*this->NyImage,sizeof(double));
	std::memcpy(this->u00,w,this->NxImage*this->NyImage*sizeof(double));
	
	this->s0Field  = new  Q0s0(this->image, this->Nx, this->Ny, this->Dx, this->Dy);	
	this->s1xField = new Q0s1x(this->image, this->Nx, this->Ny, this->Dx, this->Dy);
	this->s1yField = new Q0s1y(this->image, this->Nx, this->Ny, this->Dx, this->Dy);
	this->s2xField = new Q0s2x(this->image, this->Nx, this->Ny, this->Dx, this->Dy);
	this->s2yField = new Q0s2y(this->image, this->Nx, this->Ny, this->Dx, this->Dy);
	this->sxyField = new Q0sxy(this->image, this->Nx, this->Ny, this->Dx, this->Dy);
}
/**
 * @brief      Compute the volumic integrate of the density
 *
 * @param      int_seg  The int segment
 *
 * @return      \f$ \int_{int_seg} \vec{s} \cdot \vec{n} ?\f$ where \f$ \nabla \cdot \vec{s} = \mu $\f
 */
double IntegrateQ0::IntMass(Integration_segment * int_seg){
	return (*int_seg).diffTime*(*int_seg).direction[1]*IntegrationMethod::gauss1pt(int_seg,static_cast<VolIntegrand*>(this->s0Field));
}
/**
 * @brief      Compute the volumic integrate of the density
 *
 * @param      int_seg  The int segment
 *
 * @return      \f$ \int_{int_seg} \vec{s} \cdot \vec{n} ?\f$ where \f$ \nabla \cdot \vec{s} = \mu x $\f
 */
double IntegrateQ0::IntX1(Integration_segment * int_seg){
	return (*int_seg).diffTime*(*int_seg).direction[1]*IntegrationMethod::gauss2pt(int_seg,static_cast<VolIntegrand*>(this->s1xField));
}
/**
 * @brief      Compute the volumic integrate of the density
 *
 * @param      int_seg  The int segment
 *
 * @return      \f$ \int_{int_seg} \vec{s} \cdot \vec{n} ?\f$ where \f$ \nabla \cdot \vec{s} = \mu x^2$\f
 */
double IntegrateQ0::IntY1(Integration_segment * int_seg){
	return -(*int_seg).diffTime*(*int_seg).direction[0]*IntegrationMethod::gauss2pt(int_seg,static_cast<VolIntegrand*>(this->s1yField));
}
/**
 * @brief      Compute the volumic integrate of the density
 *
 * @param      int_seg  The int segment
 *
 * @return      \f$ \int_{int_seg} \vec{s} \cdot \vec{n} ?\f$ where \f$ \nabla \cdot \vec{s} = \mu  y^2$\f
 */
double IntegrateQ0::IntX2(Integration_segment * int_seg){
	return (*int_seg).diffTime*(*int_seg).direction[1]*IntegrationMethod::gauss2pt(int_seg,static_cast<VolIntegrand*>(this->s2xField));
}

double IntegrateQ0::IntY2(Integration_segment * int_seg){
	return -(*int_seg).diffTime*(*int_seg).direction[0]*IntegrationMethod::gauss2pt(int_seg,static_cast<VolIntegrand*>(this->s2yField));
}

double IntegrateQ0::IntXY(Integration_segment * int_seg){
	return (*int_seg).diffTime*(*int_seg).direction[1]*IntegrationMethod::gauss2pt(int_seg,static_cast<VolIntegrand*>(this->sxyField));
}
/**
 * @brief      Compute the sufarce integrate of the density
 *
 * @param      int_seg  The int segment
 *
 * @return     All the integration needed for the Hessian w.to \f$ psi $\f
 */
void IntegrateQ0::intHess(Integration_segment * int_seg, double & HessiiWW, double & HessijWW){
	int i = int_seg->indi;
	int j = int_seg->indj;
	double	tmp=u00[i*Ny+j]*int_seg->diffTime*int_seg->normDir/(int_seg->twoNormInOut);
	HessiiWW += tmp;
	HessijWW -= tmp;

}

void IntegrateQ0::intHess2(Integration_segment * int_seg, 
                         double & HessiiWX, double & HessijWX, double & HessiiWY, double & HessijWY,
                         double & HessiiXX, double & HessijXX, double & HessiiXY, double & HessijXY,double & HessijYX,
                         double & HessiiYY, double & HessijYY){
    double theta12, theta22, beta12, beta22;
    theta12 = (1 - gauss2Point1)*int_seg->thetaStartX + gauss2Point1*int_seg->thetaEndX;
    beta12 = (1 - gauss2Point1)*int_seg->thetaStartY + gauss2Point1*int_seg->thetaEndY;
    theta22 = (1 - gauss2Point2)*int_seg->thetaStartX + gauss2Point2*int_seg->thetaEndX;
    beta22 = (1 - gauss2Point2)*int_seg->thetaStartY + gauss2Point2*int_seg->thetaEndY;
    int i = int_seg->indi;
    int j = int_seg->indj;
    double f01,f02;
    f01 = u00[i*Ny+j];
    f02 = u00[i*Ny+j];
    
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
    beta13 = (1 - gauss3Point1)*int_seg->thetaStartY + gauss3Point1*int_seg->thetaEndY;
    theta23 = (1 - gauss3Point2)*int_seg->thetaStartX + gauss3Point2*int_seg->thetaEndX;
    beta23 = (1 - gauss3Point2)*int_seg->thetaStartY + gauss3Point2*int_seg->thetaEndY;
    theta33 = (1 - gauss3Point3)*int_seg->thetaStartX + gauss3Point3*int_seg->thetaEndX;
    beta33 = (1 - gauss3Point3)*int_seg->thetaStartY + gauss3Point3*int_seg->thetaEndY;

    double f001,f002,f003;
    f001 = u00[i*Ny+j];
    f002 = u00[i*Ny+j];
    f003 = u00[i*Ny+j];
    
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
    double factorized1 = 2.0*int_seg->diffTime*int_seg->normDir/(int_seg->normInOut);
	double factorized2 = (2.0/(int_seg->normInOut))*int_seg->diffTime*int_seg->normDir;
    HessiiXX += (-f2x+2.0*f1x*int_seg->xPosInside-f0*std::pow(int_seg->xPosInside,2.0))*factorized1;
    HessijXX -= (-f2x+(int_seg->xPosInside+int_seg->xPosOutside)*f1x-int_seg->xPosInside*int_seg->xPosOutside*f0)*factorized1;
    HessiiXY += factorized2*(-f1x1y+int_seg->xPosInside*f1y+int_seg->yPosInside*f1x-int_seg->xPosInside*int_seg->yPosInside*f0);
    HessijXY -= factorized2*(-f1x1y+int_seg->xPosInside*f1y+int_seg->yPosOutside*f1x-int_seg->xPosInside*int_seg->yPosOutside*f0);
	HessijYX -= factorized2*(-f1x1y+int_seg->yPosInside*f1x+int_seg->xPosOutside*f1y-int_seg->yPosInside*int_seg->xPosOutside*f0);
    HessiiYY += factorized2*(-f2y+2.0*f1y*int_seg->yPosInside-f0*std::pow(int_seg->yPosInside,2.0));
    HessijYY -= factorized2*(-f2y+(int_seg->yPosInside+int_seg->yPosOutside)*f1y-int_seg->yPosInside*int_seg->yPosOutside*f0);
}




Q0s0::Q0s0(double * u00, int Nx, int Ny, double Dx, double Dy) : 
VolIntegrandQ0(u00, Nx, Ny, Dx,Dy) {
	this->r0 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	for (int j=0; j<Ny; j++){
		for(int i = 0; i < Nx-1; i++){
			this->r0[(i+1)*Ny+j] = this->r0[i*Ny+j] +u00[i*Ny + j];			
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
double Q0s0::evaluate(Integration_segment * int_seg,double t){
	int i = int_seg->indi;
	int j = int_seg->indj;
	double theta;
	theta = (1.0-t)*int_seg->thetaStartX + t*int_seg->thetaEndX;
	return this->Dx*(theta*u00[i*Ny + j]  + this->r0[i*Ny + j]);
}

Q0s1x::Q0s1x(double * u00, int Nx, int Ny, double Dx, double Dy) : VolIntegrandQ0(u00, Nx, Ny, Dx,Dy){
	this->r0 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->dx2 = std::pow(Dx,2);
	for (int j=0; j<Ny; j++){
		for(int i = 0; i < Nx-1; i++){		
			this->r0[(i+1)*Ny+j] = this->r0[i*Ny+j] + (i+0.5)*this->u00[i*Ny+j];
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
double Q0s1x::evaluate(Integration_segment * int_seg,double t) {
	int i = int_seg->indi;
	int j = int_seg->indj;
	double theta;
	theta = (1.0-t)*int_seg->thetaStartX + t*int_seg->thetaEndX;
	return this->dx2*((i*theta+std::pow(theta,2)*.5)*u00[i*Ny + j] + this->r0[i*Ny + j]);
}

Q0s1y::Q0s1y(double * u00, int Nx, int Ny, double Dx, double Dy) : VolIntegrandQ0(u00, Nx, Ny, Dx,Dy){	
	this->r0 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->dy2 = std::pow(this->Dy,2);
	for(int i = 0; i < this->Nx; i++){
		for (int j=0; j<this->Ny-1; j++){			
			this->r0[i*Ny+j+1] = this->r0[i*Ny+j] + (j+0.5)*this->u00[i*Ny+j];		
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
double Q0s1y::evaluate(Integration_segment * int_seg,double t) {
	int i = int_seg->indi;
	int j = int_seg->indj;
	double  beta;
	beta  = (1.0-t)*int_seg->thetaStartY + t*int_seg->thetaEndY;
	return this->dy2*((j*beta+std::pow(beta,2)*.5)*u00[i*Ny + j] + this->r0[i*Ny + j]);
}

Q0s2x::Q0s2x(double * u00, int Nx, int Ny, double Dx, double Dy) : VolIntegrandQ0(u00, Nx, Ny, Dx,Dy){
	this->r0 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->dx3 = std::pow(Dx,3.0);
	for (int j=0; j<Ny; j++){
		for(int i = 0; i < Nx-1; i++){			
			this->r0[(i+1)*Ny+j] = this->r0[i*Ny+j] + (std::pow(i,2)+i+oneThird)*(this->u00[i*Ny+j]);
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
double Q0s2x::evaluate(Integration_segment * int_seg,double t){
	int i = int_seg->indi;
	int j = int_seg->indj;
	double theta;
	theta = (1.0-t)*int_seg->thetaStartX + t*int_seg->thetaEndX;
	double ThetaSquare = std::pow(theta,2);
	double thirdThetaCube = std::pow(theta,3)*oneThird;
	return dx3*((std::pow(i,2)*theta+i*ThetaSquare+thirdThetaCube)*u00[i*Ny + j]
		    + this->r0[i*Ny + j] );
}

Q0s2y::Q0s2y(double * u00, int Nx, int Ny, double Dx, double Dy) : VolIntegrandQ0(u00, Nx, Ny, Dx,Dy){
	this->r0 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->dy3 = std::pow(Dx,3.0);
	for(int i = 0; i < Nx; i++){
		for (int j=0; j<Ny-1; j++){			
			this->r0[i*Ny+j+1] = this->r0[i*Ny+j] + (std::pow(j,2)+j+oneThird)*(this->u00[i*Ny+j]);		
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
double Q0s2y::evaluate(Integration_segment * int_seg,double t){
	int i = int_seg->indi;
	int j = int_seg->indj;
	double beta;
	beta  = (1.0-t)*int_seg->thetaStartY + t*int_seg->thetaEndY;
	return dy3*((std::pow(j,2)*beta+j*std::pow(beta,2)+std::pow(beta,3)*oneThird)*u00[i*Ny + j]
		    + this->r0[i*Ny + j] );
}


Q0sxy::Q0sxy(double * u00, int Nx, int Ny, double Dx, double Dy) : VolIntegrandQ0(u00, Nx, Ny, Dx,Dy){
	this->r0 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->r1 = (double *) calloc((this->Nx)*(this->Ny),sizeof(double));
	this->dx2dy = std::pow(Dx,2.0)*Dy;
	for(int i = 0; i < Nx-1; i++){
		for (int j=0; j<Ny; j++){			
			this->r0[(i+1)*Ny+j] = this->r0[i*Ny+j] + (0.5+i)*j*(this->u00[i*Ny+j]);	
			this->r1[(i+1)*Ny+j] = this->r1[i*Ny+j] + (0.5+i)*(this->u00[i*Ny+j]);		
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
double Q0sxy::evaluate(Integration_segment * int_seg,double t){
	int i = int_seg->indi;
	int j = int_seg->indj;
	double theta,beta;
	theta  = (1.0-t)*int_seg->thetaStartX + t*int_seg->thetaEndX;
	beta   = (1.0-t)*int_seg->thetaStartY + t*int_seg->thetaEndY;
	return dx2dy*((0.5*std::pow(theta,2)+theta*i)*(beta+j)*u00[i*Ny + j]
		    + beta*this->r1[i*Ny + j] +this->r0[i*Ny + j]);
}