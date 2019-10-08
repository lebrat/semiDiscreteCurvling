#ifndef LAGUERRE_H
#define LAGUERRE_H

#include "primitives.h"
#include "edges.h"

/**
 * @brief      Storage for the Polygon of Laguerre tesselation
 */
struct Polygon2{
	std::list<Segmented_dual_edge> ListEdges;
	int VertexIndex;
	int NbTrueEdges;
    int MatrixIndex;
    Polygon2();
};

/**
 * @brief      Laguerre tesselation and info of adjacency
 */
struct Tesselation{
std::vector<Polygon2> polyList;
int nnzHessian;
int nbPolygons;
};


/**
 * @brief      Class for the computation of the Laguerre tesselation (Power diagram)
 */
class Laguerre{
    private:
        Polygon2 computePolygon(RT::Finite_vertices_iterator);
        std::list<Segmented_dual_edge> closePolygon(std::list<Segmented_dual_edge>);
     	Segmented_dual_edge clamp(Unsegmented_dual_edge);
        std::vector<double> computeCorner(int side);
        int locate(std::vector<double>);
        bool is_exactly_on_boundary(std::vector<double>);
        const double Tolx, Toly;
        const double LeftBoundary,RightBoundary,BottomBoundary,UpperBoundary;
    public:
        ~Laguerre(){
        }
        Laguerre() : Tolx(1e-13),Toly(1e-13),LeftBoundary(0.0),RightBoundary(0.0),BottomBoundary(0.0),UpperBoundary(0.0){std::cout<<"Watch out calling empty constructor in"<<__FILE__<<std::endl;};
        Laguerre(double l, double r, double b, double u, double tolx, double toly);           
        void preprocessing();
        Tesselation computeTesselation();
        std::list<std::list<double>> computeTess();
        std::list<std::list<int>> computeAdjacency();
        RT T;
};


#endif