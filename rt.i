/* rt.i */

%module rt 


%{
#define SWIG_FILE_WITH_INIT
#include "edges.h"
#include "rt.h"
%}

%ignore Unsegmented_dual_edge::Unsegmented_dual_edge(Kernel::Segment_2,int inside,int outside);
%ignore Unsegmented_dual_edge::Unsegmented_dual_edge(Kernel::Ray_2,int inside,int outside);

%include "numpy.i"
%include "std_list.i"
%include "std_vector.i"



%template(ListUnsegmentedEdge) std::list<Unsegmented_dual_edge>;
%template(ListEdge) std::list<Segmented_dual_edge>;
%template(ListPolygon) std::list<std::list<Segmented_dual_edge>>;
%template(DoubleVector) std::vector<double>;
%template(DoubleList) std::list<double>;
%template(IntList) std::list<int>;
%template(IntVector) std::vector<int>;
%template(DoubleVectorPoint) std::vector<double*>;
%template(Tesselation) std::list<std::list<double>>;
%template(adja) std::list<std::list<int>>;

%init %{
    import_array();
%}

%apply (double* IN_ARRAY1, int DIM1) {(double* vec, int n)}
%apply (int DIM1, double* IN_ARRAY1) {(int lenx, double* x), (int leny, double* y), (int lenw, double* w), (int lenpsi, double* psi)}
%apply (int DIM1, double* INPLACE_ARRAY1) {(int p, double *val)};
%apply (int DIM1, long* INPLACE_ARRAY1) {(int q, long *col),(int r, long *row)};
%apply (int DIM1, double* INPLACE_ARRAY1) {(int n, double *Bar),(int m, double *Cost),(int k, double *Mass)};
%apply (int DIM1, double* INPLACE_ARRAY1) {(int n, double *MomentX2),(int m, double *MomentY2),(int k, double *MomentXY)};



%typemap(out) double * {
 	npy_intp shape[1] = { 2 };
	$result = PyArray_SimpleNewFromData(1,shape,NPY_DOUBLE,$1);	
}
%include "edges.h"

%include "rt.h"