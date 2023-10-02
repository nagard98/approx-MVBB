#include <CGAL/convex_hull_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/min_quadrilateral_2.h>
#include <CGAL/Point_set_3/IO/OFF.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>

#include <QApplication>

#include "PointsBox.cpp"
#include "Grid.cpp"
#include "ProjectionPlane.cpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Aff_transformation_3<K> Transform3;


K::Vector_3 makeCorrectDirVec(K::Point_3 a, K::Point_3 b){
    K::Point_3 origin = K::Point_3(0,0,0);
    if((a - origin).squared_length() < (b - origin).squared_length()){
        return K::Vector_3(a, b);
    }else{
        return K::Vector_3(b, a);
    }
}

CGAL::Polyhedron_3<K> computeApproxMVBB(Points3 points, double epsilon = 0.005, double gridSize = 5) {
    //Computa axis-aligned bounding box dei punti
    K::Iso_cuboid_3 aabb1 = CGAL::bounding_box(points.begin(), points.end());

    BoxDirections3 bDirs = BoxDirections3(
        K::Vector_3(aabb1[0], aabb1[1]) * (epsilon / (2 * sqrt(3))),
        K::Vector_3(aabb1[0], aabb1[3]) * (epsilon / (2 * sqrt(3))),
        K::Vector_3(aabb1[0], aabb1[5]) * (epsilon / (2 * sqrt(3))),
        K::Point_3(0.0f, 0.0f, 0.0f));

    //Campiona punti su una griglia costruita a partire da direzioni del axis-aligned bounding box
    Grid3 grid = Grid3(bDirs);
    GridPoints3 snappedPoints3D = grid.snapToGrid(points.begin(), points.end());

    //Computa primo diametro sui punti campionati
    std::pair<K::Point_3, K::Point_3> diamEndpoints = grid.getDiameter(snappedPoints3D.begin(), snappedPoints3D.begin(), snappedPoints3D.end());
    K::Vector_3 scaledDiam = makeCorrectDirVec(diamEndpoints.first, diamEndpoints.second);

    //Proietta punti su piano perpendicolare a vettore diametro
    ProjectionPlane hPlane = ProjectionPlane(K::Point_3(0, 0, 0) + scaledDiam, scaledDiam);
    Proj2DRes resProj1 = hPlane.projectTo2D(points.begin(), points.end());

    //Computa axis-aligned bounding box dei punti proiettati su piano
    K::Iso_rectangle_2 aabb2 = CGAL::bounding_box(resProj1._projPoints.begin(), resProj1._projPoints.end());

    BoxDirections2 bDirs2 = BoxDirections2(
        K::Vector_2(aabb2[0], aabb2[1]) * (epsilon / (2 * sqrt(2))),
        K::Vector_2(aabb2[0], aabb2[3]) * (epsilon / (2 * sqrt(2))),
        K::Point_2(0.0f, 0.0f));

    //Campiona punti su una griglia costruita a partire da direzioni del axis-aligned bounding box
    Grid2 grid2 = Grid2(bDirs2);
    std::unordered_set<K::Point_2> snappedPts2D = grid2.snapToGrid(resProj1._projPoints.begin(), resProj1._projPoints.end());

    K::Point_2 testPoint = resProj1._projPoints[0];

    //Computa secondo diametro su proiezione
    K::Segment_2 diam2 = grid2.getDiameter(snappedPts2D);

    K::Vector_2 perpVec = diam2.to_vector().perpendicular(CGAL::POSITIVE);

    //Riproietta in 3D il secondo diametro computato
    K::Point_3 unpSrcDiam2 = hPlane.unprojectTo3D(diam2.source()); //hPlane.unprojectTo3D(srcDiam2Scaled);
    K::Point_3 unpEndDiam2 = hPlane.unprojectTo3D(diam2.end()); //hPlane.unprojectTo3D(endDiam2Scaled);
    K::Vector_3 scaledDiam2 = K::Vector_3(unpSrcDiam2, unpEndDiam2);


    //Proietta punti su linea perpendicolare a secondo diametro
    K::Line_2 hLine = K::Line_2(K::Point_2(0, 0) + K::Vector_2(diam2.source(), diam2.end()), perpVec);
    Points2 lineProjPoints = grid2.projectToLine(snappedPts2D.begin(), snappedPts2D.end(), hLine);

    //computa terzo diametro su punti proiettati su linea
    K::Point_2 tmpMinPoint = lineProjPoints[0], tmpMaxPoint = lineProjPoints[0];
    for (Points2::iterator it = lineProjPoints.begin(); it != lineProjPoints.end(); it++) {
        if (tmpMinPoint > *it) tmpMinPoint = *it;
        else if (tmpMaxPoint < *it) tmpMaxPoint = *it;
    }

    K::Vector_2 hLineNormalized = (hLine.to_vector() / CGAL::approximate_sqrt(hLine.to_vector().squared_length()));
    K::Vector_2 tmpDiam2 = K::Vector_2(diam2.source(), diam2.end());
    K::Vector_2 tmpDiam3 = K::Vector_2(tmpMinPoint, tmpMaxPoint);


    //Riproietta in 3D il terzo diametro computato
    K::Point_3 unpSrcDiam3 = hPlane.unprojectTo3D(tmpMinPoint);
    K::Point_3 unpEndDiam3 = hPlane.unprojectTo3D(tmpMaxPoint);
    K::Vector_3 scaledDiam3 = K::Vector_3(unpSrcDiam3, unpEndDiam3);


    //Exhaustive grid search---------------------------------------
    BoxDirections3 searchDirs = BoxDirections3(
        (epsilon / 5) * epsilon * scaledDiam,
        (epsilon / 5) * epsilon * scaledDiam2,
        (epsilon / 5) * epsilon * scaledDiam3, K::Point_3(0, 0, 0));
    Grid3 searchGrid = Grid3(searchDirs);

    std::cout << "Initial approx Volume = " << CGAL::approximate_sqrt(scaledDiam.squared_length()) * CGAL::approximate_sqrt(scaledDiam2.squared_length()) * CGAL::approximate_sqrt(scaledDiam3.squared_length()) << std::endl;

    std::vector<K::Vector_3> candidates = searchGrid.getSearchCandidates(gridSize);

    Points3 unpMinAreaPts;
    K::Vector_3 bestCandidate;
    K::Vector_3 offsetPlane;
    double finalDiam = 0;
    double bestVolume = 10000;

    for (std::vector<K::Vector_3>::iterator it = candidates.begin(); it != candidates.end(); it++) {
        ProjectionPlane candPlane = ProjectionPlane(K::Point_3(0, 0, 0), *it);

        Proj2DRes resultProj = candPlane.projectTo2D(points.begin(), points.end());

        Points2 outp;
        CGAL::ch_graham_andrew(resultProj._projPoints.begin(), resultProj._projPoints.end(), std::back_inserter(outp));

        CGAL::Polygon_2<K> minAreaRect;
        CGAL::min_rectangle_2(outp.begin(), outp.end(), std::back_inserter(minAreaRect));

        K::Segment_2 base2 = K::Segment_2(minAreaRect.vertices()[0], minAreaRect.vertices()[1]);
        K::Segment_2 height2 = K::Segment_2(minAreaRect.vertices()[1], minAreaRect.vertices()[2]);
        Points3 tmpUnpMinAreaPts = candPlane.unprojectTo3D(minAreaRect.vertices_begin(), minAreaRect.vertices_end());

        K::Segment_3 base = K::Segment_3(tmpUnpMinAreaPts[0], tmpUnpMinAreaPts[1]);
        K::Segment_3 height = K::Segment_3(tmpUnpMinAreaPts[1], tmpUnpMinAreaPts[2]);

        double tmpVol = CGAL::approximate_sqrt(base.squared_length()) * CGAL::approximate_sqrt(height.squared_length()) * resultProj._diam;
        std::cout << "Volume = " << tmpVol << std::endl;

        if (tmpVol < bestVolume) {
            bestVolume = tmpVol;
            bestCandidate = *it;
            unpMinAreaPts = tmpUnpMinAreaPts;
            offsetPlane = resultProj._offsetPlane;
            finalDiam = resultProj._diam;
        }
    }


    K::Vector_3 kv = (bestCandidate / CGAL::approximate_sqrt(bestCandidate.squared_length())) * finalDiam;
    std::cout << "Best Volume is = " << bestVolume << std::endl;

    CGAL::Polyhedron_3<K> oobb;
    K::Point_3 org = K::Point_3(0, 0, 0);
    CGAL::make_hexahedron(
        unpMinAreaPts[0] + offsetPlane,
        unpMinAreaPts[1] + offsetPlane,
        unpMinAreaPts[2] + offsetPlane,
        unpMinAreaPts[3] + offsetPlane,
        unpMinAreaPts[3] + offsetPlane + kv,
        unpMinAreaPts[0] + offsetPlane + kv,
        unpMinAreaPts[1] + offsetPlane + kv,
        unpMinAreaPts[2] + offsetPlane + kv,
        oobb);
    /*CGAL::make_hexahedron(
        org,
        org + scaledDiam,
        org + scaledDiam2 + scaledDiam,
        org + scaledDiam2,
        org + scaledDiam2 + scaledDiam3,
        org + scaledDiam3,
        org + scaledDiam + scaledDiam3,
        org + scaledDiam + scaledDiam2 + scaledDiam3,
        oobb);*/

    return oobb;
}

int main(int argc, char** argv){

    CGAL::Point_set_3<K::Point_3> point_set;
    Points3 points;

    if (argc == 2) {
        //load dragon model
        if (strcmp(argv[1], "dragon") == 0) {
            CGAL::IO::read_PLY("C:/Users/Draga/Documents/GitHub/approx-MVBB/build/Debug/dragon3.ply", point_set);

            for (int i = 0; i < point_set.number_of_points(); i++) {
                points.push_back(point_set.point(i));
            }
        }
        //generate random points in tetrahedron
        else if (strcmp(argv[1], "tetra") == 0) {
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
        CGAL::Polyhedron_3<K> polyhedron;
        polyhedron.make_tetrahedron(K::Point_3(3.45, -0.14, -6.31), K::Point_3(4.96, 0, 5.26), K::Point_3(-5.7, 0, 2.14), K::Point_3(-2.76, 5.86, 0.48));

        CGAL::Random_points_in_triangle_mesh_3<CGAL::Polyhedron_3<K>> g(polyhedron);
        std::copy_n(g, 300, std::back_inserter(points));

        for (Points3::iterator it = points.begin(); it != points.end(); it++) {
            point_set.insert(*it);
        }
    }

    CGAL::Polyhedron_3<K> oobb = computeApproxMVBB(points);

    QApplication app(argc, argv);
    PointsBox mainwindow(app.activeWindow(), point_set, oobb);
    mainwindow.show();
    app.exec();
    
    return 0;
}

