#ifndef EDGES_H
#define EDGES_H
#include "primitives.h"


/**
 * @brief      In Segmented_dual_edge we cast the <CODE>Unsegmented_dual_edge<\CODE> structure, taking care of intersecting
 * 			   the infinite edges with the box of the image \f$\mu\f$.
 */
struct Segmented_dual_edge{
	std::vector<double> source;
	std::vector<double> target;
	int inside;
	int outside;
	bool source_on_boundary;
	bool target_on_boundary;
	bool is_degenerate;
	~Segmented_dual_edge(){
		source = std::vector<double>();
		target =std::vector<double>();
	};
	void _print();
	Segmented_dual_edge();
	Segmented_dual_edge(std::vector<double> source,std::vector<double> target, int inside,int outside,bool source_on_boundary,bool target_on_boundary);
};


/**
 * @brief      In Unsegmented_dual_edge, we extract from <CODE>CGAL<\CODE> the dual of the Delaunay triangulation,
 * 			   the Laguerre tesselation, it can be expressed as a set of segment/line/ray, this structure this data.\n
 * 			   Note : the object can be infinite in this class
 * 			   
 * 
 */
struct Unsegmented_dual_edge{
	std::vector<double> point;
	std::vector<double> direction;
	bool is_infinite_tmax;
	int inside;
	int outside;
	bool is_reversed;
	Unsegmented_dual_edge();
	Unsegmented_dual_edge(Kernel::Segment_2,int inside,int outside);
	Unsegmented_dual_edge(Kernel::Ray_2,int inside,int outside,bool is_reversed);
	void _print();
	// void display(void);
};

typedef std::list<Segmented_dual_edge> Polygon;
#endif