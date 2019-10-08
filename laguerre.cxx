#include "laguerre.h"


Vertex_handle_RT get_source2(RT::Edge edge){
	return edge.first->vertex( RT::ccw(edge.second) );
}    
    

Vertex_handle_RT get_target2(RT::Edge edge){
	return edge.first->vertex( RT::cw(edge.second) );        
}
Polygon2::Polygon2(){
this->NbTrueEdges=0;
}

/**
 * @brief      Constructs the object.
 *
 * @param[in]  l     { parameter_description }
 * @param[in]  r     { parameter_description }
 * @param[in]  b     { parameter_description }
 * @param[in]  u     { parameter_description }
 * @param[in]  tolx  The tolx
 * @param[in]  toly  The toly
 */
Laguerre::Laguerre(double l, double r, double b, double u, double tolx=1.0e-13, double toly=1.0e-13) :
					Tolx(tolx), Toly(toly),LeftBoundary(l),RightBoundary(r),BottomBoundary(b),UpperBoundary(u)
					{
}


/**
 * @brief      Calculates the tesselation.
 *
 * @return     The tesselation.
 */
Tesselation Laguerre::computeTesselation(void) {
	Tesselation result;
	RT::Finite_vertices_iterator vit;
	result.nnzHessian=0;

	for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit){
		Polygon2 poly=this->computePolygon(vit);
		poly.MatrixIndex=result.nnzHessian;
		result.nnzHessian+=poly.NbTrueEdges+1;	
		result.polyList.push_back(poly);
	}
	result.nbPolygons=result.polyList.size();
	
	return result; 
};

std::list<std::list<double>> Laguerre::computeTess(){
	std::list<std::list<double>> ret;
	RT::Finite_vertices_iterator vit;
	for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit){
		Polygon2 poly=this->computePolygon(vit);
		std::list<double> loc;
		loc.push_back((double) poly.VertexIndex);
		for(std::list<Segmented_dual_edge>::iterator edgit = poly.ListEdges.begin(); edgit != poly.ListEdges.end(); ++edgit){
			loc.push_back(edgit->source[0]);
			loc.push_back(edgit->source[1]);
		}  
		ret.push_back(loc);
	}
	return ret; 
};

std::list<std::list<int>> Laguerre::computeAdjacency(){
	std::list<std::list<int>> ret;
	RT::Finite_vertices_iterator vit;
	for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit){
		Polygon2 poly=this->computePolygon(vit);
		std::list<int> voisin;
		voisin.push_back(poly.VertexIndex);
		for(std::list<Segmented_dual_edge>::iterator edgit = poly.ListEdges.begin(); edgit != poly.ListEdges.end(); ++edgit){
			if(edgit->outside >-1)
				voisin.push_back(edgit->outside);			
		}  
		ret.push_back(voisin);
	}
	return ret; 
};

void Laguerre::preprocessing(){
	RT::Finite_faces_iterator fit;
	int toRemove;
	long int index;
	double d1x,d2x,d3x,d1y,d2y,d3y,determinant,d1square,d2square,d3square,distancemax;
	double wRemove, xRemove, yRemove;
	for (fit = T.finite_faces_begin(); fit != T.finite_faces_end(); ++fit){		
		
		d1x = fit->vertex(1)->point().x() - fit->vertex(0)->point().x();
		d1y = fit->vertex(1)->point().y() - fit->vertex(0)->point().y();

		d2x = fit->vertex(2)->point().x() - fit->vertex(0)->point().x();
		d2y = fit->vertex(2)->point().y() - fit->vertex(0)->point().y();		

		determinant = d1x*d2y - d2x*d1y;

		if (std::abs(determinant) < 1.0e-32){
			d3x = fit->vertex(2)->point().x() - fit->vertex(1)->point().x();
			d3y = fit->vertex(2)->point().y() - fit->vertex(1)->point().y();

			d1square = std::pow(d1x,2.0) + std::pow(d1y,2.0);
			d2square = std::pow(d2x,2.0) + std::pow(d2y,2.0);
			d3square = std::pow(d3x,2.0) + std::pow(d3y,2.0);
			distancemax = std::max(d1square,std::max(d2square,d3square));
			if(distancemax == d1square){
				toRemove = 2;
			}
			else if (distancemax == d2square){
				toRemove = 1;
			}
			else if (distancemax == d3square){
				toRemove = 0;
			}
			else{
				std::cout<<"\033[1;31m Outch unexpected case in preprocessing"<<"\033[0m"<<std::endl;
				toRemove = 1;
			}
			wRemove = fit->vertex(toRemove)->point().weight();
			xRemove = fit->vertex(toRemove)->point().x();
			yRemove = fit->vertex(toRemove)->point().y();
			index = fit->vertex(toRemove)->info();
			T.remove(fit->vertex(toRemove));
			std::vector<std::pair<Weighted_point, long int> > points_tmp(1);
			points_tmp.push_back(std::make_pair(Weighted_point(Point(xRemove,yRemove),wRemove),index));
			T.insert(points_tmp.begin(),points_tmp.end());			
			throw std::invalid_argument( "Non-admissible triangulation" );
		}
	}
}


/**
 * @brief      Calculates the polygon.
 *
 * @param[in]  vit   The vit
 *
 * @return     The polygon.
 */
Polygon2 Laguerre::computePolygon(RT::Finite_vertices_iterator vit){
	Polygon2 poly;
	std::list<Unsegmented_dual_edge> ListUnsegmentedDualEdges;
	std::list<Segmented_dual_edge> prepolygon;
	RT::Edge_circulator ec = T.incident_edges(vit);
	RT::Edge_circulator eend = ec;
	CGAL_For_all(ec, eend)
	{
		RT::Edge edge = *ec;
		Vertex_handle_RT source=get_source2(edge);
		Vertex_handle_RT target=get_target2(edge);
		int inside = target->info();
		int outside = source->info();		
		CGAL::Object o ;

		if(T.is_infinite(edge)){
				continue;
		}
		try{		
			o=T.dual(edge);
		}
		catch(...){
			std::cout<<"3 points are aligned in the triangle (iside = "<<inside<<" outside = "<<outside<<std::endl;
		}
		Kernel::Segment_2 s;
		Kernel::Ray_2 r;
		if (CGAL::assign(s,o)){
			Unsegmented_dual_edge e(s,target->info(),source->info());	         
			ListUnsegmentedDualEdges.push_back(e);
			Segmented_dual_edge se=this->clamp(e);
			if (!(se.is_degenerate)){				
				prepolygon.push_back(se);
			};
		}
		if (CGAL::assign(r,o))
		{
			int where=source->info();
			if (!(where==-1))
			{								
				double normale=-r.direction().dy()*(target->point().point().x()-source->point().point().x() );
				normale+=r.direction().dx()*(target->point().point().y()-source->point().point().y() );				
				bool is_reversed = (normale <0);				
				Unsegmented_dual_edge e(r,target->info(),source->info(),is_reversed);
				ListUnsegmentedDualEdges.push_back(e);
				Segmented_dual_edge se=this->clamp(e);
				if (!(se.is_degenerate)){
					prepolygon.push_back(se);
				}
			}
		}
	}
	if ( prepolygon.size())
	{
		poly.ListEdges=closePolygon(prepolygon);
		poly.VertexIndex=vit->info();
		for (std::list<Segmented_dual_edge>::const_iterator e = poly.ListEdges.begin(), end = poly.ListEdges.end(); e != end; ++e) {
			if ((*e).outside >-1) {
				poly.NbTrueEdges++;
			}
		}
	}
	return poly;
}

/**
 * @brief      { function_description }
 *
 * @param[in]  e     { parameter_description }
 *
 * @return     { description_of_the_return_value }
 */
Segmented_dual_edge Laguerre::clamp(Unsegmented_dual_edge e){
	// if (e.inside==67){
	// 	e._print();
	// }
	Segmented_dual_edge se;
	se.inside=e.inside;
	se.outside=e.outside;
	se.is_degenerate=false;
	double tmin=0.0;
	double tmax=1.;
	bool is_infinite_tmax=e.is_infinite_tmax;
	bool is_vertical=(std::abs(e.direction[0]) <Tolx);
	bool is_horizontal=(std::abs(e.direction[1]) <Toly);
	if  (is_vertical && is_horizontal){
		se.is_degenerate=true;
		return se;		
	}
	else{
		if (is_vertical) {
			if ((e.point[0]<LeftBoundary) || (e.point[0]>RightBoundary)){
				se.is_degenerate=true;
				return se;
			}
		}
		else{
			double t1=(LeftBoundary-e.point[0])/e.direction[0];
			double t2=(RightBoundary-e.point[0])/e.direction[0];
			tmin=std::max(tmin,std::min(t1,t2));
			if (is_infinite_tmax){
				tmax=std::max(t1,t2);
				is_infinite_tmax=false;
			}
			else{
				tmax=std::min(tmax,std::max(t1,t2));
			}
		}
		if (is_horizontal) {
			if ((e.point[1]<BottomBoundary) || (e.point[1]>UpperBoundary)) {
				se.is_degenerate=true;
				return se;
			}
		}
		else{
			double t1=(BottomBoundary-e.point[1])/e.direction[1];
			double t2=(UpperBoundary-e.point[1])/e.direction[1];
			tmin=std::max(tmin,std::min(t1,t2));
			if (is_infinite_tmax){
				tmax=std::max(t1,t2);
				is_infinite_tmax=false;
			}
			else{
				tmax=std::min(tmax,std::max(t1,t2));
			}
		}
		if (tmax<=tmin){
			se.is_degenerate=true;
			return se;
		}
	}
	
	if (!(se.is_degenerate)){
		std::vector<double> source(2);
		std::vector<double> target(2);
		for (int i=0;i<2;i++){
			source[i]=e.point[i]+tmin*e.direction[i];
			target[i]=e.point[i]+tmax*e.direction[i];
		}
		bool SourceOnBoundary = !(tmin ==0.);
		bool TargetOnBoundary = ((is_infinite_tmax) || !(tmax==1.));
		if (e.is_reversed){
			se.source=target;
			se.target=source;
			se.source_on_boundary=TargetOnBoundary;
			se.target_on_boundary=SourceOnBoundary;
		}
		else{
			se.source=source;
			se.target=target;
			se.source_on_boundary=SourceOnBoundary;
			se.target_on_boundary=TargetOnBoundary;
		}
	}
	// if (se.inside==67){
	// 	std::cout<<se.outside<<" "<<tmin<<" "<<tmax<<std::endl;
	// 	se._print();
	// }
	return se;
}

/**
 * @brief      Closes a polygon.
 *
 * @param[in]  prepolygon  The prepolygon
 *
 * @return     { description_of_the_return_value }
 */
std::list<Segmented_dual_edge> Laguerre::closePolygon(std::list<Segmented_dual_edge> prepolygon){
	//throw std::invalid_argument( "Last side is new Side" );
	std::list<Segmented_dual_edge> polygon;
	Segmented_dual_edge last_e;
	last_e = prepolygon.back();
	for (std::list<Segmented_dual_edge>::const_iterator e = prepolygon.begin(), end = prepolygon.end(); e != end; ++e) {
		Segmented_dual_edge actual_e=*e;
		if (actual_e.source_on_boundary != last_e.target_on_boundary){
			if (!actual_e.source_on_boundary){actual_e.source_on_boundary=this->is_exactly_on_boundary(actual_e.source) ;}
			if (!last_e.target_on_boundary){  last_e.target_on_boundary=this->is_exactly_on_boundary(  last_e.target) ;}
			if (actual_e.source_on_boundary != last_e.target_on_boundary){
				std::cout<<"RAISING A VALUE ERROR, only one of the source of the actual edge or the target of the previous edge are on the boundary."<<std::endl;
				last_e._print();
				actual_e._print();
				std::abort();
			}
		}
		if (actual_e.source_on_boundary){
			int last_side = this->locate(last_e.target);
			int actual_side = this->locate(actual_e.source);
			if (last_side==-1 || actual_side==-1){
				std::cout<<"MUST RAISE AN ERROR HERE"<<std::endl;
				throw std::invalid_argument( "Last side is new Side" );
			}
			std::vector<double> last_point=last_e.target;
			while (last_side != actual_side){
				last_side = (last_side + 1) % 4;
				std::vector<double>  corner = this->computeCorner(last_side);
				Segmented_dual_edge new_e(last_point,corner,actual_e.inside,-1,1,1);
				polygon.push_back(new_e);
				last_point=corner;          	
			}
			Segmented_dual_edge new_e(last_point,actual_e.source,actual_e.inside,-1,1,1);
			polygon.push_back(new_e);
		}
		polygon.push_back(actual_e);
		last_e=actual_e;	
	}
	return polygon;
}

/**
 * @brief      Calculates the corner.
 *
 * @param[in]  side  The side
 *
 * @return     The corner.
 */
std::vector<double> Laguerre::computeCorner(int side){
	std::vector<double> vec(2);
	if (side==0){
		static const double arr[] = {LeftBoundary,BottomBoundary};
		vec.assign(&arr[0],&arr[0]+2);
		return vec;
	}	
	if (side==3){
		static const double arr[] = {LeftBoundary,UpperBoundary};
		vec.assign(&arr[0],&arr[0]+2);		
		return vec;
	}
	if (side==2){
		static const double arr[] = {RightBoundary,UpperBoundary};
		vec.assign(&arr[0],&arr[0]+2);		
		return vec;
	}
	if (side==1){
		static const double arr[] = {RightBoundary,BottomBoundary};
		vec.assign(&arr[0],&arr[0]+2);
		return vec;
	}
	static const double arr[] = {0.0,0.0};
	vec.assign(&arr[0],&arr[0]+2);
	return vec;
}

/**
 * @brief      Locate which edge to add
 *
 * @param[in]  p     { parameter_description }
 *
 * @return     { description_of_the_return_value }
 */
int Laguerre::locate(std::vector<double> p){
	int side = -1;
	if (std::abs(p[0] -   LeftBoundary) < Tolx) side = 3;
	if (std::abs(p[1] -  UpperBoundary) < Toly) side = 2;
	if (std::abs(p[0] -  RightBoundary) < Tolx) side = 1;
	if (std::abs(p[1] - BottomBoundary) < Toly) side = 0;
	return side;
}

bool Laguerre::is_exactly_on_boundary(std::vector<double> p){
	bool Isboundary=false;
	if  (std::abs(p[0] -   LeftBoundary) < Tolx) Isboundary=true;
	if  (std::abs(p[1] -  UpperBoundary) < Toly) Isboundary=true;
	if  (std::abs(p[0] -  RightBoundary) < Tolx) Isboundary=true;
	if  (std::abs(p[1] - BottomBoundary) < Toly) Isboundary=true;
	return Isboundary;
}
