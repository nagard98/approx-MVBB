#include <CGAL/convex_hull_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/min_quadrilateral_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>

#include "Grid.cpp"
#include "ProjectionPlane.cpp"

typedef struct Caliper {
    int index1, index2;
    K::Point_2 pi1, pi2;
    K::Vector_2 direction;
}Caliper;

K::Vector_3 makeCorrectDirVec(K::Point_3 a, K::Point_3 b);

namespace BoundingBox {
	CGAL::Polygon_2<K> computeMABB(Points2 points);

    CGAL::Polyhedron_3<K> computeApproxMVBB(Points3 pts, double epsilon = 0.005, double gridSize = 5);
}