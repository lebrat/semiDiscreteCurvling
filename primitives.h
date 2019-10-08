#ifndef PRIMITIVES_H
#define PRIMITIVES_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_vertex_base_2.h>
// #include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <cstring>
// CF
// https://doc.cgal.org/latest/Triangulation_2/Triangulation_2_2info_insert_with_pair_iterator_regular_2_8cpp-example.html#_a1
typedef CGAL::Exact_predicates_inexact_constructions_kernel                 Kernel;
typedef CGAL::Regular_triangulation_vertex_base_2<Kernel>                   Vbase;
typedef CGAL::Triangulation_vertex_base_with_info_2<long int, Kernel,Vbase> Vb;
typedef CGAL::Regular_triangulation_face_base_2<Kernel>                     Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                         Tds;
typedef CGAL::Regular_triangulation_2<Kernel, Tds>                          RT;
typedef Kernel::Point_2                                                     Point;
typedef Kernel::Weighted_point_2                                            Weighted_point;
typedef RT::Vertex_handle Vertex_handle_RT;

// typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
// typedef CGAL::Triangulation_vertex_base_with_info_2<long int,Kernel> Vb;
// typedef Kernel::FT                                               weight;
// typedef CGAL::Triangulation_face_base_2<Kernel> Fb;
// typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
// typedef CGAL::Regular_Triangulation_2<Kernel,Tds> RT;
// typedef Kernel::Weighted_point_2 Weighted_point;
// typedef CGAL::Point_2<Kernel> point;


#endif