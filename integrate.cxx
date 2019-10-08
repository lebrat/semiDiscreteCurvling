#include "integrate.h"

/**
 * @brief      Constructs the base VolIntegrand Object. It will help us storing the data to perform 
 * 			   the volumic integration computation for momentum from order 0  up to 2.
 * 
 * @param[in]  Nx    the Number of base Pixel/Rectangle in the \f$e_x\f$ direction
 * @param[in]  Ny    the Number of base Pixel/Rectangle in the \f$e_y\f$ direction
 * @param[in]  Dx    the distance between two consecutive neighbour in the \f$e_x\f$ direction
 * @param[in]  Dy    the distance between two consecutive neighbour in the \f$e_y\f$ direction
 */
VolIntegrand::VolIntegrand(int Nx,int Ny,double Dx,double Dy){
	this->Nx = Nx;
	this->Ny = Ny;
	this->Dx = Dx;
	this->Dy = Dy;
}

/**
 * @brief      Constructs base object integration method, will realize the gauss quadrature method
 * 			   to compute all the derivated needed represented by its abstract function :\n
 * 			   		\t-IntMass\n
 * 			   		\t-IntX1 and IntY1 \n
 * 			   		\t-IntX2 and IntY2 \n
 * 			   		\t-intHess
 * 			   
 *
 * @param[in]  LefB  The left boundary
 * @param[in]  RigB  The right boundary
 * @param[in]  BotB  The bottom boundary
 * @param[in]  UppB  The upper boundary
 */
IntegrationMethod::IntegrationMethod(double LefB,double RigB,double BotB,double UppB) : LeftBoundary(0.0), RightBoundary(1.0), BottomBoundary(0.0), UpperBoundary(1.0){
	this->oneHalf = 0.5;
	this->oneThird = 1.0/3.0;
	this->oneFourth = 0.25;
	this->gauss1Point = 0.5;
	this->gauss2Point1 = 0.5*(1.0-1.0/std::sqrt(3.0));
	this->gauss2Point2 = 0.5*(1.0+1.0/std::sqrt(3.0));
	this->gauss3Point1 = 0.5*(1.0-std::sqrt(3.0/5.0));
	this->gauss3Point2 = 0.5;
	this->gauss3Point3 = 0.5*(1.0+std::sqrt(3.0/5.0));
	this->gauss1Factor = 1.0;
	this->gauss2Factor1 = 0.5;
	this->gauss2Factor2 = 0.5;
	this->gauss3Factor1 = 2.5/9.0;
	this->gauss3Factor2 = 4.0/9.0;
	this->gauss3Factor3 = 2.5/9.0;
}

/**
 * @brief      Compute the value exact value of \f$\int_0^1 f(t) dt\f$ for \f$f(t)\f$ 
 *			   a polynomial with a degree up to 1. 			    
 *
 * @param      int_seg  The integrated segment
 * @param      v        A certain vector field \f$ \vec{s}\f$ s.t. \f$ \nabla \cdot \vec{s}  = \mu \f$
 *
 * @return     the value of the integral
 */
double IntegrationMethod::gauss1pt(Integration_segment * int_seg,VolIntegrand* v){
	return this->gauss1Factor*v->evaluate(int_seg,this->gauss1Point);
}

/**
 * @brief      Compute the value exact value of \f$\int_0^1 f(t) dt\f$ for \f$f(t)\f$ 
 *			   a polynomial with a degree up to 3. 			    
 *
 * @param      int_seg  The integrated segment
 * @param      v        A certain vector field \f$ \vec{s}\f$ s.t. \f$ \nabla \cdot \vec{s}  = \mu \f$
 *
 * @return     the value of the integral
 */
double IntegrationMethod::gauss2pt(Integration_segment * int_seg,VolIntegrand* v){
	return this->gauss2Factor1*v->evaluate(int_seg,this->gauss2Point1) 
		 + this->gauss2Factor2*v->evaluate(int_seg,this->gauss2Point2);
}

/**
 * @brief      Compute the value exact value of \f$\int_0^1 f(t) dt\f$ for \f$f(t)\f$ 
 *			   a polynomial with a degree up to 5. 			    
 *
 * @param      int_seg  The integrated segment
 * @param      v        A certain vector field \f$ \vec{s}\f$ s.t. \f$ \nabla \cdot \vec{s}  = \mu \f$
 *
 * @return     the value of the integral
 */
double IntegrationMethod::gauss3pt(Integration_segment * int_seg,VolIntegrand* v){
	return this->gauss3Factor1*v->evaluate(int_seg,this->gauss3Point1) 
		 + this->gauss3Factor2*v->evaluate(int_seg,this->gauss3Point2)
		 + this->gauss3Factor3*v->evaluate(int_seg,this->gauss3Point3);
}

/**
 * @brief      Return the value of the grided image for indexes resquested
 *
 * @param[in]  i     the column index
 * @param[in]  j     the row index 
 *
 * @return     \f$  \mu[i,j] \f$
 */
double IntegrationMethod::img_ij(int i , int j){
	return this->image[i*this->NyImage + j];
}


/**
 * @brief      Intersect and edge with the below image  \f$  \mu \f$
 *
 * @param      e     The egde
 * @param[in]  Tolx  The toloreance to determine if an edhge is vertical
 * @param[in]  Toly  The toloreance to determine if an edhge is horizontal
 *
 * @return     { description_of_the_return_value }
 */
std::list<double> IntegrationMethod::sample(Segmented_dual_edge & e,double Tolx,double Toly){
	bool is_vertical=(std::abs(e.target[0]-e.source[0]) <Tolx);
	bool is_horizontal=(std::abs(e.target[1]-e.source[1]) <Toly);
	int imin= ((int) ((std::min(e.source[0],e.target[0]) - BottomBoundary)/Dx));
	int imax= ((int) ((std::max(e.source[0],e.target[0]) - BottomBoundary)/Dx))+1;
	int jmin= ((int) ((std::min(e.source[1],e.target[1]) - LeftBoundary)/Dy));
	int jmax= ((int) ((std::max(e.source[1],e.target[1]) - LeftBoundary)/Dy))+1;
	std::list<double> time;
	if (!(is_vertical)){
		for (int i = imin; i < imax; i++)
		{
			time.push_back((i*Dx+BottomBoundary-e.source[0])/(e.target[0]-e.source[0]));	
		}
	}
	if (!(is_horizontal)){
		for (int j = jmin; j < jmax; j++)
		{
			time.push_back((j*Dy+LeftBoundary-e.source[1])/(e.target[1]-e.source[1]));	
		}
	}	
	time.sort();
	while (!time.empty() && !(time.front()>0.)) {time.pop_front();}
	while (!time.empty() && !(time.back()<1.)) {time.pop_back();}
	time.push_front(0.);
	time.push_back(1.);
	return time;
}

/**
 * @brief      Convert a section of the edge on a constant pattern of the underlying measure in local
 *			   coordinate (pixel/Q1) 
 *
 * @param      int_seg   The int segment to fill
 * @param      se        The global segmented dual edge
 * @param[in]  current   The entry time in this the pixel
 * @param[in]  previous  The leaving time in this pixel
 * 
 */
void IntegrationMethod::precompute(Integration_segment * int_seg,Segmented_dual_edge &se, double current,double previous){
	int_seg->diffTime=current - previous;
	int_seg->cutPosStartX=(1-previous)*(se).source[0] + previous*(se).target[0];
	int_seg->cutPosStartY=(1-previous)*(se).source[1] + previous*(se).target[1];
	int_seg->cutPosEndX  =(1-current )*(se).source[0] +  current*(se).target[0];
	int_seg->cutPosEndY  =(1-current )*(se).source[1] +  current*(se).target[1];
	int_seg->indi=std::max(std::min((int)((0.5*(int_seg->cutPosStartX + int_seg->cutPosEndX) - this->LeftBoundary  )/this->Dx),this->Nx-1),0);
	int_seg->indj=std::max(std::min((int)((0.5*(int_seg->cutPosStartY + int_seg->cutPosEndY) - this->BottomBoundary)/this->Dy),this->Ny-1),0);
	int_seg->thetaEndX  =(int_seg->cutPosStartX - this->LeftBoundary  )/this->Dx - int_seg->indi;
	int_seg->thetaEndY  =(int_seg->cutPosStartY - this->BottomBoundary)/this->Dy - int_seg->indj;
	int_seg->thetaStartX=(int_seg->cutPosEndX   - this->LeftBoundary  )/this->Dx - int_seg->indi;
	int_seg->thetaStartY=(int_seg->cutPosEndY   - this->BottomBoundary)/this->Dy - int_seg->indj;
	int_seg->thetaX=0.5*(int_seg->thetaStartX+int_seg->thetaEndX);
	int_seg->thetaY=0.5*(int_seg->thetaStartY+int_seg->thetaEndY);
}