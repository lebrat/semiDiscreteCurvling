#include "edges.h"

//********************************************** Constructors / Destructors for Unsegmented_dual_edge
Unsegmented_dual_edge::Unsegmented_dual_edge() {
}

/**
 * @brief      Constructs the <CODE>Unsegmented_dual_edge<\CODE>
 *
 * @param[in]  s        segment from cGal
 * @param[in]  inside   Inside
 * @param[in]  outside  The outside
 */
Unsegmented_dual_edge::Unsegmented_dual_edge(Kernel::Segment_2 s,int inside,int outside) {
	this->is_infinite_tmax=false;
	this->inside=inside;
	this->outside=outside;
	this->is_reversed=false;
	std::vector<double> point(2);
	std::vector<double> direction(2);
	
	if (std::abs(s.source().x())+std::abs(s.source().y())>1e6*(std::abs(s.target().x())+std::abs(s.target().y()))){
		this->is_reversed=true;
		point[0]=s.target().x();
		point[1]=s.target().y();
		direction[0]=-s.to_vector().x();
		direction[1]=-s.to_vector().y();
	}
	else{
		point[0]=s.source().x();
		point[1]=s.source().y();
		direction[0]=s.to_vector().x();
		direction[1]=s.to_vector().y();
	}
	this->point=point;
	this->direction=direction;
}

/**
 * @brief      Constructs the <CODE>Unsegmented_dual_edge<\CODE>
 *
 * @param[in]  r            ray from cGal
 * @param[in]  inside       Inside
 * @param[in]  outside      The outside
 * @param[in]  is_reversed  Indicates if reversed
 */
Unsegmented_dual_edge::Unsegmented_dual_edge(Kernel::Ray_2 r,int inside,int outside,bool is_reversed) {
	this->is_infinite_tmax=true;
	this->inside=inside;
	this->outside=outside;
	this->is_reversed=is_reversed;
	std::vector<double> point(2);
	point[0]=r.source().x();
	point[1]=r.source().y();
	this->point=point;
	std::vector<double> direction(2);
	direction[0]=r.direction().dx();
	direction[1]=r.direction().dy();
	this->direction=direction;
}


void Unsegmented_dual_edge::_print(){
	std::cout<<"From ("<<this->point[0]<<" , "<<this->point[1]<< ")";
	std::cout<<" To ("<<this->point[0]+this->direction[0]<<" , "<<this->point[1]+this->direction[1]<< ")";
	std::cout<<" in/out "<<this->inside<<"/"<<this->outside<< " infty? "<<this->is_infinite_tmax<<" reversed? "<<this->is_reversed;
	std::cout<<std::endl;
}

Segmented_dual_edge::Segmented_dual_edge(){
}

/**
 * @brief      Constructs the <CODE>Segmented_dual_edge<\CODE>.
 *
 * @param[in]  source              The source
 * @param[in]  target              The target
 * @param[in]  inside              Inside
 * @param[in]  outside             The outside
 * @param[in]  source_on_boundary  The source on boundary
 * @param[in]  target_on_boundary  The target on boundary
 */
Segmented_dual_edge::Segmented_dual_edge(std::vector<double> source,std::vector<double> target, int inside,int outside, bool source_on_boundary,bool target_on_boundary){
	this->source=source;
	this->target=target;
	this->inside=inside;
	this->outside=outside;
	this->source_on_boundary=source_on_boundary;
	this->target_on_boundary=target_on_boundary;
}

void Segmented_dual_edge::_print(){
	std::cout<<"From ("<<this->source[0]<<" , "<<this->source[1]<< ")["<<this->source_on_boundary <<"]";
	std::cout<<" To ("<<this->target[0]<<" , "<<this->target[1]<< ")["<<this->target_on_boundary <<"]";
	std::cout<<" in/out "<<this->inside<<"/"<<this->outside;
	std::cout<<std::endl;
}