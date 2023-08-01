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


int main(int argc, char** argv){

    std::ofstream f("test.off");
    CGAL::Point_set_3<K::Point_3> point_set;

    Points3 points, result;
    CGAL::Polyhedron_3<K> polyhedron;
    polyhedron.make_tetrahedron(K::Point_3(3.45, -0.14, -6.31), K::Point_3(4.96,0,5.26), K::Point_3(-5.7,0,2.14), K::Point_3(-2.76,5.86,0.48));

    CGAL::Random_points_in_triangle_mesh_3<CGAL::Polyhedron_3<K>> g(polyhedron);
    std::copy_n(g, 200, std::back_inserter(points));

    for(Points3::iterator it = points.begin(); it!= points.end(); it++){
        point_set.insert(*it);
    }

    //CGAL::IO::write_OFF(f, point_set);
    CGAL::IO::write_OFF(f, polyhedron);

    double const epsilon = 0.005;
    double const gridSize = 5;

    K::Iso_cuboid_3 aabb1 = CGAL::bounding_box(points.begin(), points.end());

    //CGAL::IO::write_OFF(f, aabbPolyhedron);

    BoxDirections3 bDirs = BoxDirections3(
        K::Vector_3(aabb1[0], aabb1[1]) * (epsilon/(2*sqrt(3))),
        K::Vector_3(aabb1[0], aabb1[3]) * (epsilon/(2*sqrt(3))),
        K::Vector_3(aabb1[0], aabb1[5]) * (epsilon/(2*sqrt(3))),
        K::Point_3(0.0f,0.0f,0.0f));

    Grid3 grid = Grid3(bDirs);
    GridPoints3 snappedPoints3D = grid.snapToGrid(points.begin(), points.end());

    std::pair<K::Point_3, K::Point_3> diamEndpoints = grid.getDiameter(snappedPoints3D.begin(), snappedPoints3D.begin(), snappedPoints3D.end());
    K::Vector_3 diam = makeCorrectDirVec(diamEndpoints.first, diamEndpoints.second);
    K::Vector_3 scaledDiam = K::Vector_3(
        diam.x() * bDirs._lenX,
        diam.y() * bDirs._lenY, 
        diam.z() * bDirs._lenZ);


    ProjectionPlane hPlane = ProjectionPlane(K::Point_3(0,0,0) + scaledDiam, scaledDiam);
    Proj2DRes resProj1 = hPlane.projectTo2D(points.begin(), points.end());

    //TODO: convert second diameter calculation to 2D
    K::Iso_rectangle_2 aabb2 = CGAL::bounding_box(resProj1._projPoints.begin(), resProj1._projPoints.end());

    BoxDirections2 bDirs2 = BoxDirections2(
        K::Vector_2(aabb2[0], aabb2[1]) * (epsilon/(2*sqrt(2))),
        K::Vector_2(aabb2[0], aabb2[3]) * (epsilon/(2*sqrt(2))), 
        K::Point_2(0.0f, 0.0f));

    Grid2 grid2 = Grid2(bDirs2);
    std::unordered_set<K::Point_2> snappedPts2D = grid2.snapToGrid(resProj1._projPoints.begin(), resProj1._projPoints.end());
    K::Segment_2 diam2 = grid2.getDiameter(snappedPts2D);
    
    K::Point_2 srcDiam2Scaled = K::Point_2(diam2.source().x() * bDirs2._lenX, diam2.source().y() * bDirs2._lenY);
    K::Point_2 endDiam2Scaled = K::Point_2(diam2.end().x() * bDirs2._lenX, diam2.end().y() * bDirs2._lenY);

    K::Vector_2 perpVec = diam2.to_vector().perpendicular(CGAL::POSITIVE);
    K::Point_3 unpSrcDiam2 = hPlane.unprojectTo3D(srcDiam2Scaled);
    K::Point_3 unpEndDiam2 = hPlane.unprojectTo3D(endDiam2Scaled);
    K::Vector_3 scaledDiam2 = K::Vector_3(unpSrcDiam2, unpEndDiam2);

    K::Line_2 hLine = K::Line_2(K::Point_2(0,0) + K::Vector_2(diam2.source(), diam2.end()) , perpVec);

    Points2 lineProjPoints = grid2.projectToLine(snappedPts2D.begin(), snappedPts2D.end(), hLine);


    K::Point_2 tmpMinPoint = lineProjPoints[0], tmpMaxPoint = lineProjPoints[0];
    for(Points2::iterator it = lineProjPoints.begin(); it != lineProjPoints.end(); it++){
        if(tmpMinPoint > *it) tmpMinPoint = *it;
        else if(tmpMaxPoint < *it) tmpMaxPoint = *it;
    }

    std::cout << "angle in 2d: " << hLine.to_vector() * K::Vector_2(diam2.source(), diam2.end()) << std::endl;
    std::cout << "angle in 2d: " <<K::Vector_2(tmpMinPoint, tmpMaxPoint) * diam2.to_vector() << std::endl; // Per ora risulta 90°
    std::cout << "angle in 2d scaled: " << K::Vector_2(K::Point_2(tmpMinPoint.x() * bDirs2._lenX , tmpMinPoint.y() * bDirs2._lenY), K::Point_2(tmpMaxPoint.x() * bDirs2._lenX, tmpMaxPoint.y()  * bDirs2._lenY)) * K::Vector_2(srcDiam2Scaled, endDiam2Scaled) << std::endl;
    K::Vector_2 tmp3 = K::Vector_2(K::Point_2(tmpMinPoint.x() * bDirs2._lenX , tmpMinPoint.y() * bDirs2._lenY), K::Point_2(tmpMaxPoint.x() * bDirs2._lenX, tmpMaxPoint.y()  * bDirs2._lenY));
    std::cout << "Length in 2d: " << tmp3.squared_length() << std::endl;


    K::Point_3 unpSrcDiam3 = hPlane.unprojectTo3D( K::Point_2(tmpMinPoint.x() * bDirs2._lenX , tmpMinPoint.y() * bDirs2._lenY));
    K::Point_3 unpEndDiam3 = hPlane.unprojectTo3D( K::Point_2(tmpMaxPoint.x() * bDirs2._lenX, tmpMaxPoint.y()  * bDirs2._lenY));
    double scaledDiam3Len = CGAL::approximate_sqrt(K::Vector_3(unpSrcDiam3, unpEndDiam3).squared_length());
    K::Vector_3 diam3 = CGAL::cross_product(scaledDiam,scaledDiam2);
    //TODO: fix third diam angle to diam angle 2
    K::Vector_3 scaledDiam3 = /* K::Vector_3(unpSrcDiam3, unpEndDiam3); */(diam3 / CGAL::approximate_sqrt(diam3.squared_length())) * scaledDiam3Len;

    std::cout << "angle: " << CGAL::approximate_angle(scaledDiam, scaledDiam2) << std::endl;
    std::cout << "angle: " << CGAL::approximate_angle(scaledDiam, scaledDiam3) << std::endl;
    std::cout << "angle: " << CGAL::approximate_angle(scaledDiam3, scaledDiam2) << std::endl;


    //Exhaustive grid search---------------------------------------
    BoxDirections3 searchDirs = BoxDirections3(
        (epsilon/5) * epsilon * scaledDiam,
        (epsilon/5) * epsilon * scaledDiam2,
        (epsilon/5) * epsilon * scaledDiam3, K::Point_3(0,0,0));
    Grid3 searchGrid = Grid3(searchDirs);

    std::vector<K::Vector_3> candidates = searchGrid.getSearchCandidates(gridSize);

    Points3 unpMinAreaPts;
    K::Vector_3 bestCandidate;
    K::Vector_3 offsetPlane;
    double finalDiam = 0;
    double bestVolume = 10000;

    for(std::vector<K::Vector_3>::iterator it = candidates.begin(); it != candidates.end(); it++){
        ProjectionPlane candPlane = ProjectionPlane(K::Point_3(0,0,0), *it);

        Proj2DRes resultProj = candPlane.projectTo2D(points.begin(), points.end());

        Points2 outp;
        CGAL::ch_graham_andrew(resultProj._projPoints.begin(), resultProj._projPoints.end(), std::back_inserter(outp));

        CGAL::Polygon_2<K> minAreaRect;
        CGAL::min_rectangle_2(outp.begin(), outp.end(), std::back_inserter(minAreaRect));

        Points3 tmpUnpMinAreaPts = candPlane.unprojectTo3D(minAreaRect.vertices_begin(), minAreaRect.vertices_end());

        K::Segment_3 base = K::Segment_3(tmpUnpMinAreaPts[0], tmpUnpMinAreaPts[1]);
        K::Segment_3 height = K::Segment_3(tmpUnpMinAreaPts[1], tmpUnpMinAreaPts[2]);
        double tmpVol = CGAL::approximate_sqrt(base.squared_length()) * CGAL::approximate_sqrt(height.squared_length()) * CGAL::approximate_sqrt(it->squared_length());
        std::cout << "Volume = " <<  tmpVol << std::endl;
        
        if(tmpVol < bestVolume){
            bestVolume = tmpVol;
            bestCandidate = *it;
            unpMinAreaPts = tmpUnpMinAreaPts;
            offsetPlane = resultProj._offsetPlane;
            finalDiam = resultProj._diam;
        }
    }


    K::Vector_3 kv = (bestCandidate / CGAL::approximate_sqrt(bestCandidate.squared_length())) * finalDiam;
    std::cout << "Best Volume is = " <<  bestVolume << std::endl;

    CGAL::Polyhedron_3<K> oobb;
    K::Point_3 org = K::Point_3(0,0,0);
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
/*     CGAL::make_hexahedron(
        org,
        org + scaledDiam, 
        org + scaledDiam2 + scaledDiam, 
        org + scaledDiam2, 
        org + scaledDiam2 + scaledDiam3, 
        org + scaledDiam3, 
        org + scaledDiam + scaledDiam3, 
        org + scaledDiam + scaledDiam2 + scaledDiam3, 
        oobb); */

    QApplication app(argc,argv);
    PointsBox mainwindow(app.activeWindow(), point_set, oobb);
    mainwindow.show();
    app.exec();

    
    return 0;
}

