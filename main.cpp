#include <CGAL/Polyhedron_3.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/point_generators_3.h>

#include <iostream>
#include <fstream>

#include "BoundingBox.hpp"

#ifndef QT_FOUND
#define QT_FOUND (0)
#endif // !QT_FOUND

#if QT_FOUND == 1
#include <QApplication>
#include "PointsBox.cpp"
#endif

typedef std::vector<K::Point_3> Points3;
typedef std::unordered_set<K::Point_3> GridPoints3;

typedef std::vector<K::Point_2> Points2;
typedef std::unordered_set<K::Point_2> GridPoints2;

int main(int argc, char** argv){

    CGAL::Point_set_3<K::Point_3> point_set;
    Points3 points;

    if (argc == 2) {
        // carica modello dragon.ply
        if (strcmp(argv[1], "-dragon") == 0) {
            CGAL::IO::read_PLY("../dragon.ply", point_set);

            for (int i = 0; i < point_set.number_of_points(); i++) {
                points.push_back(point_set.point(i));
            }
        }
        // genera punti casuali in tetraedro
        else if (strcmp(argv[1], "-tetra") == 0) {
            CGAL::Polyhedron_3<K> polyhedron;
            polyhedron.make_tetrahedron(K::Point_3(3.45, -0.14, -6.31), K::Point_3(4.96, 0, 5.26), K::Point_3(-5.7, 0, 2.14), K::Point_3(-2.76, 5.86, 0.48));

            CGAL::Random_points_in_triangle_mesh_3<CGAL::Polyhedron_3<K>> g(polyhedron);
            std::copy_n(g, 300, std::back_inserter(points));

            for (Points3::iterator it = points.begin(); it != points.end(); it++) {
                point_set.insert(*it);
            }
        }
    }
    else {
        // Di default genera punti casuali su tetraedro
        CGAL::Polyhedron_3<K> polyhedron;
        polyhedron.make_tetrahedron(K::Point_3(3.45, -0.14, -6.31), K::Point_3(4.96, 0, 5.26), K::Point_3(-5.7, 0, 2.14), K::Point_3(-2.76, 5.86, 0.48));

        CGAL::Random_points_in_triangle_mesh_3<CGAL::Polyhedron_3<K>> g(polyhedron);
        std::copy_n(g, 300, std::back_inserter(points));

        for (Points3::iterator it = points.begin(); it != points.end(); it++) {
            point_set.insert(*it);
        }
    }

    CGAL::Polyhedron_3<K> oobb = BoundingBox::computeApproxMVBB(points);
    
#if QT_FOUND==1
    QApplication app(argc, argv);
    PointsBox mainwindow(app.activeWindow(), point_set, oobb);
    mainwindow.resize(800, 600);
    mainwindow.show();
    app.exec();
#endif 
    
    return 0;
}

