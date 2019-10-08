#include "rt.h"
#include <omp.h>

#include <sys/time.h>


#define NbPROC   4  

Integration_result2::Integration_result2(){
	this->val = std::list<double>();
	this->colInd = std::list<int>();
	this->rowInd = std::list<int>();
}

void Integration_result2::set_matrix(int col,int row,double value){
	this->val.push_back(value);
	this->colInd.push_back(col);
	this->rowInd.push_back(row);
};


double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

//SETTERS

/**
 * @brief      Sets the coordinates.
 *
 * @param[in]  lenx  The lenx
 * @param      x     Vector of X positions
 * @param[in]  leny  The leny
 * @param      y     Vector of Y positions
 */

Triangulation2::Triangulation2(){
	this->ImageInstanciated=false;
	this->CoordInstanciated=false;
}

void Triangulation2::setImage(int Nx,int Ny,bool isQ1, int lenw, double* w,double l,double r,double b,double u){	
	if (isQ1) {
		IntegrateQ1 * integ = new IntegrateQ1(l,r,b,u,Nx,Ny,lenw,w);
		this->Integ = static_cast<IntegrationMethod*>(integ);
	}
	if (!isQ1) {
		IntegrateQ0 * integ = new IntegrateQ0(l,r,b,u,Nx,Ny,lenw,w);
		this->Integ = static_cast<IntegrationMethod*>(integ);
	}
	Tolx=1.e-12*this->Integ->dx();
	Toly=1.e-12*this->Integ->dy();
	Laguerre * L=new Laguerre(l,r,b,u,Tolx,Toly);
	this->Lag = L;
	this->ImageInstanciated=true;

}


void Triangulation2::setCoordinates(int lenx, double *x, int leny, double *y){
	if(leny != lenx){
		throw std::invalid_argument( "X and Y vectors should have the same size");
	}
	this->nDir = lenx;
	this->points.reserve(this->nDir);
	// std::cout<<x[0]<<std::endl;
	for (int i = 0; i < this->nDir; ++i)
	{
		this->points[i] = Weighted_point(Point(x[i],y[i]),0);
	}
	this->CoordInstanciated=true;
}

double Triangulation2::FindWopt(int lenw, double * w, int leny, double *y,double alpha){
	if(lenw != this->nDir){
		throw std::invalid_argument("Number of weigth should match the size of Points");
	}
	std::vector<std::pair<Weighted_point,int> > points_tmp(this->nDir);
	for (int i=0; i<this->nDir; i++){
		points_tmp[i]=std::make_pair(Weighted_point(this->points[i].point(),w[i]), i);
	}
	this->Lag->T.clear();
	this->Lag->T.insert(points_tmp.begin(),points_tmp.end());
	this->Lag->T.infinite_vertex()->info() = -1;
	int hidden = this->Lag->T.number_of_hidden_vertices();
	double * p1 = new double [2];
	double * p2 = new double [2];
	double * p3 = new double [2];
	double * ph = new double [2];
	double * dp = new double [2];
	double * D0 = new double [2];
	double * Q = new double [2];
	double * Am1 = new double [4];


	double w1,w2,w3,wh;
	double norm2p1,detA,buff1,buff2;
	int i1,i2,i3,ih;
	Vertex_handle_RT neigh;
	double localAlphaOpt,globalAlphaOpt;
	globalAlphaOpt = alpha;
	if(hidden){
		RT::Hidden_vertices_iterator vit;		
		for (vit = this->Lag->T.hidden_vertices_begin(); vit != this->Lag->T.hidden_vertices_end(); ++vit)
		{
			RT::Face_handle F = this->Lag->T.locate(vit->point());
			Point dualPoint = this->Lag->T.dual(F);
			dp[0] = dualPoint.x();
			dp[1] = dualPoint.y();
			
			ph[0] = vit->point().x();
			ph[1] = vit->point().y();						
			ih = vit->info();
			wh = vit->point().weight() - alpha*y[ih];
			
			neigh = F->vertex(0);
			i1 = neigh->info();
			p1[0] = neigh->point().x();
			p1[1] = neigh->point().y();
			w1 = neigh->point().weight() - alpha*y[i1];
			neigh = F->vertex(1);
			i2 = neigh->info();
			p2[0] = neigh->point().x();
			p2[1] = neigh->point().y();
			w2 = neigh->point().weight() - alpha*y[i2];
			neigh = F->vertex(2);
			i3 = neigh->info();
			p3[0] = neigh->point().x();
			p3[1] = neigh->point().y();
			w3 = neigh->point().weight() - alpha*y[i3];

			detA = (p2[0] - p1[0])*(p3[1] - p1[1]) - (p2[1] - p1[1])*(p3[0] - p1[0]);
			Am1[0] = 0.5*(p1[1] - p3[1])/detA;
			Am1[1] = 0.5*(p2[1] - p1[1])/detA;
			Am1[2] = 0.5*(p3[0] - p1[0])/detA;
			Am1[3] = 0.5*(p1[0] - p2[0])/detA;
			
			norm2p1 = pow(p1[0],2.0) + pow(p1[1],2.0);
			buff1 = norm2p1 - w1  - (pow(p2[0],2.0) + pow(p2[1],2.0)) + w2;
			buff2 = norm2p1 - w1  - (pow(p3[0],2.0) + pow(p3[1],2.0)) + w3;
			D0[0] = Am1[0]*buff1 + Am1[1]*buff2;
			D0[1] = Am1[2]*buff1 + Am1[3]*buff2;

			buff1 = y[i2] - y[i1];
			buff2 = y[i3] - y[i1];
			Q[0] = Am1[0]*buff1 + Am1[1]*buff2;
			Q[1] = Am1[2]*buff1 + Am1[3]*buff2;

			localAlphaOpt = (2*(D0[0]*(p1[0]-ph[0])+D0[1]*(p1[1]-ph[1])) + w1 - wh - norm2p1 + pow(ph[0],2.0) + pow(ph[1],2.0))/(-2.0*(Q[0]*(p1[0]-ph[0])+Q[1]*(p1[1]-ph[1]))-y[i1]+y[ih]);			
			globalAlphaOpt = std::min(globalAlphaOpt,localAlphaOpt);
		}						
	}
	
	
	if(hidden + this->Lag->T.number_of_vertices() != ((unsigned int)this->nDir)){
		std::cout<<"\033[1;31m Merged vertices "<<this->nDir - this->Lag->T.number_of_vertices()<<"\033[0m\n"<<std::endl;
	}
	
	delete[] p1;delete[] p2;delete[] p3;delete[] ph;delete[] dp;delete[] D0; delete[] Q;delete[] Am1;

	return globalAlphaOpt;
}

int Triangulation2::computeLaguerre(int lenpsi, double * psi){
	if (!this->ImageInstanciated){
		throw std::invalid_argument( "No Image instanciated");
	}	
	if (!this->CoordInstanciated){
		throw std::invalid_argument( "No coordinates instanciated");
	}
	if(lenpsi != this->nDir){
		throw std::invalid_argument("Number of weigth should match the size of Points");
	}
	std::vector<std::pair<Weighted_point, long int> > points_tmp(this->nDir);
	for (long int i=0; i<this->nDir; i++){
		points_tmp[i]=std::make_pair(Weighted_point(this->points[i].point(),psi[i]), i);
	}
	this->Lag->T.clear();
	// double timeExc = get_wall_time();
	this->Lag->T.insert(points_tmp.begin(),points_tmp.end());
	// timeExc = get_wall_time() -timeExc;
	// std::cout<<timeExc<<std::endl;
	this->Lag->T.infinite_vertex()->info() = -1;
	int hidden = this->Lag->T.number_of_hidden_vertices();
	
	if(hidden){
		// std::cout<<"\033[1;31m hidden vertices ("<<hidden<<")\033[0m\n"<<std::endl;
	}
	if(hidden + this->Lag->T.number_of_vertices() != ((unsigned int)this->nDir)){
		std::cout<<"\033[1;31m Merged vertices "<<this->nDir - this->Lag->T.number_of_vertices()<<"\033[0m\n"<<std::endl;
	}

	return hidden;
}


void Triangulation2::Integration(int n, double *Moment1,int m, double *Moment2,int k, double *Mass, int p, double *val, int q , long *col, int r, long *row){
	RT::Finite_vertices_iterator vit;
	std::list<Unsegmented_dual_edge> ListUnsegmentedDualEdges;
	Integration_segment int_seg;
	int inside;
	double mass,moment1X,moment1Y,moment2;
	double current, prev;
	double HessiiWW,HessijWW;
	// double timeExc = get_wall_time();
	for (std::vector<Polygon2>::iterator poly = this->Tess.polyList.begin(); poly != this->Tess.polyList.end(); ++poly){

		inside = poly->VertexIndex;
		if(inside > this->nDir-1) continue;
		if(inside < 0) continue;
		mass = 0.0; moment1X = 0.0; moment1Y = 0.0; moment2 = 0.0;
		HessiiWW=0.;
		int nnzEdge=1;
		for(std::list<Segmented_dual_edge>::iterator se = (*poly).ListEdges.begin(); se != (*poly).ListEdges.end(); se++)
		{	
			int outside = (*se).outside;

			int_seg.direction[0]=(*se).target[0] - (*se).source[0];
			int_seg.direction[1]=(*se).target[1] - (*se).source[1];
			int_seg.normDir = sqrt(pow(int_seg.direction[0],2.0) + pow(int_seg.direction[1],2.0));
			if( outside >= 0 && outside <this->nDir){
			int_seg.twoNormInOut = 2.0*sqrt(pow(this->points[inside].x() - this->points[outside].x(),2.0)
								+pow(this->points[inside].y() - this->points[outside].y(),2.0));
			}
			std::list<double> time = this->Integ->sample(*se,Tolx,Toly);

			std::list<double>::iterator it = time.begin();
			it ++;
			prev = time.front();
			HessijWW=0.;
			for (; it != time.end(); it++)
			{
				current = * it;
				Integ->precompute(&int_seg,*se,current,prev);
				prev = current;	
					
				mass += Integ->IntMass(&int_seg);
				moment1X += Integ->IntX1(&int_seg);
				moment1Y += Integ->IntY1(&int_seg);
				moment2 += Integ->IntX2(&int_seg) + Integ->IntY2(&int_seg);

				if (outside > -1){
					Integ->intHess(&int_seg, HessiiWW, HessijWW);
					
				}

			}	

			if( outside >= 0 && outside <this->nDir){
				if(poly->MatrixIndex+nnzEdge>=p){
				std::cout<<"Raise Pb MatIndex = "<<poly->MatrixIndex<<"+"<<nnzEdge<<" "<<p<<" Bugg"<<std::endl;

				}

				col[poly->MatrixIndex+nnzEdge]=outside;
				row[poly->MatrixIndex+nnzEdge]=inside;
				val[poly->MatrixIndex+nnzEdge]=HessijWW;				
				nnzEdge++;
			}		
		}
		if(poly->MatrixIndex<p){
			col[poly->MatrixIndex]=inside;
			row[poly->MatrixIndex]=inside;
			val[poly->MatrixIndex]=HessiiWW;
		}
		else{
			std::cout<<"Raise Pb MatIndex = "<<poly->MatrixIndex<<" "<<p<<" Bugg"<<std::endl;
			std::abort();
		}
		Mass[inside] = mass;
		Moment1[2*inside]= moment1X;
		Moment1[2*inside+1]= moment1Y;
		Moment2[inside] =  moment2;
	}
	// timeExc = get_wall_time() -timeExc;
	// std::cout<<timeExc<<std::endl;
}


void Triangulation2::IntegrationParall(int n, double *Moment1,int m, double *Moment2,int k, double *Mass, int p, double *val, int q , long *col, int r, long *row){
	
	std::vector<Polygon2> vector_poly{ std::make_move_iterator(std::begin(Tess.polyList)), 
											std::make_move_iterator(std::end(Tess.polyList)) };
	std::vector<int>::size_type vector_poly_size = vector_poly.size();
	#pragma omp parallel for num_threads(NbPROC)
	for (unsigned int it_poly = 0; it_poly<vector_poly_size; it_poly++){
		Polygon2 poly=vector_poly[it_poly];		
		std::list<Unsegmented_dual_edge> ListUnsegmentedDualEdges;
		Integration_segment int_seg;
		int inside;
		double mass,moment1X,moment1Y,moment2;
		double current, prev;
		inside = poly.VertexIndex;
		if(inside > this->nDir-1) continue;
		if(inside < 0) continue;
		mass = 0.0; moment1X = 0.0; moment1Y = 0.0; moment2 = 0.0;
		double HessiiWW=0.;
		int nnzEdge=1;
		for(std::list<Segmented_dual_edge>::iterator se = poly.ListEdges.begin(); se != poly.ListEdges.end(); se++)
		{	
			int outside = (*se).outside;

			int_seg.direction[0]=(*se).target[0] - (*se).source[0];
			int_seg.direction[1]=(*se).target[1] - (*se).source[1];
			int_seg.normDir = sqrt(pow(int_seg.direction[0],2.0) + pow(int_seg.direction[1],2.0));
			if( outside >= 0 && outside <this->nDir){
			int_seg.twoNormInOut = 2.0*sqrt(pow(this->points[inside].x() - this->points[outside].x(),2.0)
								+pow(this->points[inside].y() - this->points[outside].y(),2.0));
			}
			std::list<double> time = this->Integ->sample(*se,Tolx,Toly);

			std::list<double>::iterator it = time.begin();
			it ++;
			prev = time.front();
			double HessijWW=0.;
			for (; it != time.end(); it++)
			{
			    current = * it;
				Integ->precompute(&int_seg,*se,current,prev);
				prev = current;	
				mass += Integ->IntMass(&int_seg);
				moment1X += Integ->IntX1(&int_seg);
				moment1Y += Integ->IntY1(&int_seg);
				moment2 += Integ->IntX2(&int_seg) + Integ->IntY2(&int_seg);

				if (outside > -1){
					Integ->intHess(&int_seg, HessiiWW, HessijWW);			
				}

			}	

			if( outside >= 0 && outside <this->nDir){
				col[poly.MatrixIndex+nnzEdge]=outside;
				row[poly.MatrixIndex+nnzEdge]=inside;
				val[poly.MatrixIndex+nnzEdge]=HessijWW;				
				nnzEdge++;
			}		
		}
		if(poly.MatrixIndex<p){
			col[poly.MatrixIndex]=inside;
			row[poly.MatrixIndex]=inside;
			val[poly.MatrixIndex]=HessiiWW;
		}
		Mass[inside] = mass;
		Moment1[2*inside]= moment1X;
		Moment1[2*inside+1]= moment1Y;
		Moment2[inside] =  moment2;
	}
}



class Integration_result_on_edge{
public :
double HessijXX,HessijYY,HessijXY,HessijYX;
double HessiiXX,HessiiYY,HessiiXY;
double HessijWW,HessiiWW;
double HessijWX,HessijWY;
double HessiiWX,HessiiWY;

void init_on_edge(){
	HessijXX=0.;HessijYY=0.;HessijXY=0.;HessijYX=0.;
	HessijWW=0.;
	HessijWX=0.;HessijWY=0.;
}
void init_on_polygon(){
	HessiiXX=0.;HessiiYY=0.;HessiiXY=0.;
	HessiiWW=0.;
	HessiiWX=0.;HessiiWY=0.;
}
};


void Triangulation2::HessianCompute2(int p,double *val, int q , long *col, int r, long *row){
	RT::Finite_vertices_iterator vit;
	std::list<Unsegmented_dual_edge> ListUnsegmentedDualEdges;
	Integration_segment int_seg;
	Integration_result_on_edge integr;	
	int inside;
	double current, prev;
	int stride=2*this->Tess.nnzHessian;
	for (std::vector<Polygon2>::iterator poly = this->Tess.polyList.begin(); poly != this->Tess.polyList.end(); ++poly){
		double mass=0.0;
		inside = poly->VertexIndex;
		if(inside > this->nDir-1) continue;
		if(inside < 0) continue;
		integr.init_on_polygon();
		int nnzEdge=1;
		for(std::list<Segmented_dual_edge>::iterator se = (*poly).ListEdges.begin(); se != (*poly).ListEdges.end(); se++)
		{	
			integr.init_on_edge();
			int outside = (*se).outside;
			int_seg.direction[0]=(*se).target[0] - (*se).source[0];
			int_seg.direction[1]=(*se).target[1] - (*se).source[1];
			int_seg.normDir = sqrt(pow(int_seg.direction[0],2.0) + pow(int_seg.direction[1],2.0));

			if( outside >= 0 && outside <this->nDir){
				int_seg.twoNormInOut = 2.0*sqrt(pow(this->points[inside].x() - this->points[outside].x(),2.0)
								+pow(this->points[inside].y() - this->points[outside].y(),2.0));
				int_seg.xPosOutside=this->points[outside].x();
				int_seg.yPosOutside=this->points[outside].y();
				int_seg.normInOut = 0.5*int_seg.twoNormInOut;
			}
			int_seg.xPosInside=this->points[inside].x();
			int_seg.yPosInside=this->points[inside].y();	
			std::list<double> time = this->Integ->sample(*se,Tolx,Toly);
			std::list<double>::iterator it = time.begin();
			it ++;
			prev = time.front();

			for (; it != time.end(); it++)
			{
				current = * it;
				Integ->precompute(&int_seg,*se,current,prev);
				prev = current;	
				mass += Integ->IntMass(&int_seg);

				if ( outside >= 0 && outside <this->nDir){
					Integ->intHess2(&int_seg,integr.HessiiWX, integr.HessijWX,
					 						 integr.HessiiWY, integr.HessijWY,integr.HessiiXX, integr.HessijXX,
											 integr.HessiiXY, integr.HessijXY,integr.HessijYX,integr.HessiiYY, integr.HessijYY);
				}
			}
	
			if( outside >= 0 && outside <this->nDir){
				if(2*stride+2*(poly->MatrixIndex+nnzEdge)+1>=p){
				std::cout<<"Raise Pb MatIndex0 = "<<poly->MatrixIndex<<"+"<<nnzEdge<<" "<<stride<<" "<<p<<" Bugg"<<std::endl;
				std::abort();
				}

				row[         2*(poly->MatrixIndex+nnzEdge)  ]=inside;
				col[         2*(poly->MatrixIndex+nnzEdge)  ]=outside;
				val[         2*(poly->MatrixIndex+nnzEdge)  ]=-integr.HessijWX;
				row[         2*(poly->MatrixIndex+nnzEdge)+1]=inside;
				col[         2*(poly->MatrixIndex+nnzEdge)+1]=outside+nDir;
				val[         2*(poly->MatrixIndex+nnzEdge)+1]=-integr.HessijWY;
				row[  stride+2*(poly->MatrixIndex+nnzEdge)  ]=inside;
				col[  stride+2*(poly->MatrixIndex+nnzEdge)  ]=outside;
				val[  stride+2*(poly->MatrixIndex+nnzEdge)  ]=integr.HessijXX;
				row[  stride+2*(poly->MatrixIndex+nnzEdge)+1]=inside;
				col[  stride+2*(poly->MatrixIndex+nnzEdge)+1]=outside+nDir;
				val[  stride+2*(poly->MatrixIndex+nnzEdge)+1]=integr.HessijXY;
				row[2*stride+2*(poly->MatrixIndex+nnzEdge)  ]=inside+nDir;
				col[2*stride+2*(poly->MatrixIndex+nnzEdge)  ]=outside;
				val[2*stride+2*(poly->MatrixIndex+nnzEdge)  ]=integr.HessijYX;
				row[2*stride+2*(poly->MatrixIndex+nnzEdge)+1]=inside+nDir;
				col[2*stride+2*(poly->MatrixIndex+nnzEdge)+1]=outside+nDir;
				val[2*stride+2*(poly->MatrixIndex+nnzEdge)+1]=integr.HessijYY;		
				nnzEdge++;
				}		
		}
		if(2*stride+2*(poly->MatrixIndex)+1>=p){
				std::cout<<"Raise Pb MatIndex1 = "<<poly->MatrixIndex<<" "<<stride<<" "<<p<<" Bugg"<<std::endl;
				std::abort();
		}

		row[         2*(poly->MatrixIndex)  ]=inside;
		col[         2*(poly->MatrixIndex)  ]=inside;
		val[         2*(poly->MatrixIndex)  ]=-integr.HessiiWX;
		row[         2*(poly->MatrixIndex)+1]=inside;
		col[         2*(poly->MatrixIndex)+1]=inside+nDir;
		val[         2*(poly->MatrixIndex)+1]=-integr.HessiiWY;
		row[  stride+2*(poly->MatrixIndex)  ]=inside;
		col[  stride+2*(poly->MatrixIndex)  ]=inside;
		val[  stride+2*(poly->MatrixIndex)  ]=2*mass+integr.HessiiXX;
		row[  stride+2*(poly->MatrixIndex)+1]=inside;
		col[  stride+2*(poly->MatrixIndex)+1]=inside+nDir;
		val[  stride+2*(poly->MatrixIndex)+1]=integr.HessiiXY;
		row[2*stride+2*(poly->MatrixIndex)  ]=inside+nDir;
		col[2*stride+2*(poly->MatrixIndex)  ]=inside;
		val[2*stride+2*(poly->MatrixIndex)  ]=integr.HessiiXY;
		row[2*stride+2*(poly->MatrixIndex)+1]=inside+nDir;
		col[2*stride+2*(poly->MatrixIndex)+1]=inside+nDir;
		val[2*stride+2*(poly->MatrixIndex)+1]=2*mass+integr.HessiiYY;
		}
}

void Triangulation2::HessianComputeParall(int p,double *val, int q , long *col, int r, long *row){
	int stride=2*this->Tess.nnzHessian;
	#pragma omp parallel for num_threads(NbPROC)
	for (int it_poly = 0; it_poly<Tess.nbPolygons; it_poly++){
		Polygon2 poly=Tess.polyList[it_poly];		
		Integration_segment int_seg;
		Integration_result_on_edge integr;	
		int inside;
		double mass=0.0;
		double current, prev;
		inside = poly.VertexIndex;
		if(inside > this->nDir-1) continue;
		if(inside < 0) continue;
		integr.init_on_polygon();
		int nnzEdge=1;
		for(std::list<Segmented_dual_edge>::iterator se = poly.ListEdges.begin(); se != poly.ListEdges.end(); se++)
		{	
			integr.init_on_edge();
			int outside = (*se).outside;
			int_seg.direction[0]=(*se).target[0] - (*se).source[0];
			int_seg.direction[1]=(*se).target[1] - (*se).source[1];
			int_seg.normDir = sqrt(pow(int_seg.direction[0],2.0) + pow(int_seg.direction[1],2.0));
			if( outside >= 0 && outside <this->nDir){
				int_seg.twoNormInOut = 2.0*sqrt(pow(this->points[inside].x() - this->points[outside].x(),2.0)
								+pow(this->points[inside].y() - this->points[outside].y(),2.0));
				int_seg.xPosOutside=this->points[outside].x();
				int_seg.yPosOutside=this->points[outside].y();
				int_seg.normInOut = 0.5*int_seg.twoNormInOut;
			}
			int_seg.xPosInside=this->points[inside].x();
			int_seg.yPosInside=this->points[inside].y();
			
			std::list<double> time = this->Integ->sample(*se,Tolx,Toly);
			std::list<double>::iterator it = time.begin();
			it ++;
			prev = time.front();
			for (; it != time.end(); it++)
			{
				current = * it;
				Integ->precompute(&int_seg,*se,current,prev);
				prev = current;	
				mass += Integ->IntMass(&int_seg);

				if (outside > -1){
					Integ->intHess2(&int_seg,integr.HessiiWX, integr.HessijWX,
					 						 integr.HessiiWY, integr.HessijWY,integr.HessiiXX, integr.HessijXX,
											 integr.HessiiXY, integr.HessijXY,integr.HessijYX,integr.HessiiYY, integr.HessijYY);
				}
			}	
			if( outside >= 0 && outside <this->nDir){
				row[         2*(poly.MatrixIndex+nnzEdge)  ]=inside;
				col[         2*(poly.MatrixIndex+nnzEdge)  ]=outside;
				val[         2*(poly.MatrixIndex+nnzEdge)  ]=-integr.HessijWX;
				row[         2*(poly.MatrixIndex+nnzEdge)+1]=inside;
				col[         2*(poly.MatrixIndex+nnzEdge)+1]=outside+nDir;
				val[         2*(poly.MatrixIndex+nnzEdge)+1]=-integr.HessijWY;
				row[  stride+2*(poly.MatrixIndex+nnzEdge)  ]=inside;
				col[  stride+2*(poly.MatrixIndex+nnzEdge)  ]=outside;
				val[  stride+2*(poly.MatrixIndex+nnzEdge)  ]=integr.HessijXX;
				row[  stride+2*(poly.MatrixIndex+nnzEdge)+1]=inside;
				col[  stride+2*(poly.MatrixIndex+nnzEdge)+1]=outside+nDir;
				val[  stride+2*(poly.MatrixIndex+nnzEdge)+1]=integr.HessijXY;
				row[2*stride+2*(poly.MatrixIndex+nnzEdge)  ]=inside+nDir;
				col[2*stride+2*(poly.MatrixIndex+nnzEdge)  ]=outside;
				val[2*stride+2*(poly.MatrixIndex+nnzEdge)  ]=integr.HessijYX;
				row[2*stride+2*(poly.MatrixIndex+nnzEdge)+1]=inside+nDir;
				col[2*stride+2*(poly.MatrixIndex+nnzEdge)+1]=outside+nDir;
				val[2*stride+2*(poly.MatrixIndex+nnzEdge)+1]=integr.HessijYY;		
				nnzEdge++;
				}		
		}
		row[         2*(poly.MatrixIndex)  ]=inside;
		col[         2*(poly.MatrixIndex)  ]=inside;
		val[         2*(poly.MatrixIndex)  ]=-integr.HessiiWX;
		row[         2*(poly.MatrixIndex)+1]=inside;
		col[         2*(poly.MatrixIndex)+1]=inside+nDir;
		val[         2*(poly.MatrixIndex)+1]=-integr.HessiiWY;
		row[  stride+2*(poly.MatrixIndex)  ]=inside;
		col[  stride+2*(poly.MatrixIndex)  ]=inside;
		val[  stride+2*(poly.MatrixIndex)  ]=2*mass+integr.HessiiXX;
		row[  stride+2*(poly.MatrixIndex)+1]=inside;
		col[  stride+2*(poly.MatrixIndex)+1]=inside+nDir;
		val[  stride+2*(poly.MatrixIndex)+1]=integr.HessiiXY;
		row[2*stride+2*(poly.MatrixIndex)  ]=inside+nDir;
		col[2*stride+2*(poly.MatrixIndex)  ]=inside;
		val[2*stride+2*(poly.MatrixIndex)  ]=integr.HessiiXY;
		row[2*stride+2*(poly.MatrixIndex)+1]=inside+nDir;
		col[2*stride+2*(poly.MatrixIndex)+1]=inside+nDir;
		val[2*stride+2*(poly.MatrixIndex)+1]=2*mass+integr.HessiiYY;	
		}
}

void Triangulation2::ComputeLocalMatrix(int n, double *MomentX2,int m, double *MomentY2,int k, double *MomentXY){
	RT::Finite_vertices_iterator vit;
	std::list<Unsegmented_dual_edge> ListUnsegmentedDualEdges;
	Integration_segment int_seg;
	int inside;
	double mass,moment1X,moment1Y,moment2X,moment2Y,momentXY;
	double current, prev;
	for (std::vector<Polygon2>::iterator poly = this->Tess.polyList.begin(); poly != this->Tess.polyList.end(); ++poly){
		inside = poly->VertexIndex;
		if(inside > this->nDir-1) continue;
		if(inside < 0) continue;
		mass = 0.0; moment1X = 0.0; moment1Y = 0.0; moment2X = 0.0; moment2Y = 0.0;momentXY = 0.0;
		for(std::list<Segmented_dual_edge>::iterator se = (*poly).ListEdges.begin(); se != (*poly).ListEdges.end(); se++)
		{
			int_seg.direction[0]=(*se).target[0] - (*se).source[0];
			int_seg.direction[1]=(*se).target[1] - (*se).source[1];
			int_seg.normDir = sqrt(pow(int_seg.direction[0],2.0) + pow(int_seg.direction[1],2.0));
			std::list<double> time = this->Integ->sample(*se,Tolx,Toly);
			std::list<double>::iterator it = time.begin();
			it ++;
			prev = time.front();
			for (; it != time.end(); it++)
			{
				current = * it;
				Integ->precompute(&int_seg,*se,current,prev);
				prev = current;
				
				mass += Integ->IntMass(&int_seg);
				moment1X += Integ->IntX1(&int_seg);
				moment1Y += Integ->IntY1(&int_seg);
				moment2X += Integ->IntX2(&int_seg);
				moment2Y += Integ->IntY2(&int_seg);
				momentXY += Integ->IntXY(&int_seg);
			}	
		}
		if (mass > 0.)
		{
			MomentX2[inside]= moment2X-std::pow(moment1X,2)/mass;
			MomentY2[inside]= moment2Y-std::pow(moment1Y,2)/mass;
			MomentXY[inside]= momentXY-   moment1X*moment1Y/mass;
		}
	}
}
