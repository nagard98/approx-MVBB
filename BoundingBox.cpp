#include "BoundingBox.hpp"

typedef CGAL::Aff_transformation_3<K> Transform3;

#ifdef NDEBUG
#define _LOG(x);
#else
#define _LOG(x) std::cout << x << std::endl;
#endif


K::Vector_3 makeCorrectDirVec(K::Point_3 a, K::Point_3 b) {
    K::Point_3 origin = K::Point_3(0, 0, 0);
    if ((a - origin).squared_length() < (b - origin).squared_length()) {
        return K::Vector_3(a, b);
    }
    else {
        return K::Vector_3(b, a);
    }
}



namespace BoundingBox {
    CGAL::Polygon_2<K> computeMABB(Points2 points) {
        Caliper calipers[4];
        double angles[4] = { 0.0 };
        int pIndx[4] = { 0 };

        double tmpMaxX = points[0].x();
        double tmpMinX = points[0].x();
        double tmpMaxY = points[0].y();
        double tmpMinY = points[0].y();

        double tmpX, tmpY;
        for (int i = 0; i < points.size(); i++) {
            tmpX = points[i].x();
            tmpY = points[i].y();

            if (tmpMaxX <= tmpX) {
                tmpMaxX = tmpX;
                pIndx[2] = i;
            }
            else if (tmpMinX >= tmpX) {
                tmpMinX = tmpX;
                pIndx[0] = i;
            }

            if (tmpMaxY <= tmpY) {
                tmpMaxY = tmpY;
                pIndx[1] = i;
            }
            else if (tmpMinY >= tmpY) {
                tmpMinY = tmpY;
                pIndx[3] = i;
            }
        }

        for (int i = 0; i < 4; i++) {
            calipers[i].index1 = pIndx[i];
            calipers[i].index2 = calipers[i].index1;
            calipers[i].pi1 = points[pIndx[i]];
            calipers[i].pi2 = calipers[i].pi1;
        }
        calipers[0].direction = K::Vector_2(0, -1);
        calipers[1].direction = K::Vector_2(-1, 0);
        calipers[2].direction = K::Vector_2(0, 1);
        calipers[3].direction = K::Vector_2(1, 0);


        CGAL::Polygon_2<K> minRect;
        double minArea = CGAL_IA_MAX_DOUBLE;
        for (int j = 0; j < points.size(); j++) {
            int tmpMinAngleIndx = 0;
            for (int i = 0; i < 4; i++) {
                K::Vector_2 candDir = points[(calipers[i].index2 + 1) % points.size()] - points[calipers[i].index2];
                K::Vector_2 normCandDir = candDir / CGAL::approximate_sqrt(candDir.squared_length());
                angles[i] = std::acos(normCandDir * calipers[i].direction);
                if (angles[i] < angles[tmpMinAngleIndx]) tmpMinAngleIndx = i;
            }

            double rotAngle = angles[tmpMinAngleIndx];
            K::Aff_transformation_2 rotate(CGAL::ROTATION, sin(rotAngle), cos(rotAngle));

            for (int i = 0; i < 4; i++) {
                if (i == tmpMinAngleIndx) {
                    calipers[i].index1 = calipers[i].index2;
                    calipers[i].index2 = (calipers[i].index2 + 1) % points.size();
                }
                calipers[i].direction = rotate(calipers[i].direction);
            }

            CGAL::Polygon_2<K> tmpRect;
            for (int i = 0; i < 4; i++) {
                //std::cout << "angle directions is : " << calipers[i].direction * calipers[(i + 1) % 4].direction << std::endl;
                K::Line_2 line1 = K::Line_2(points[calipers[i].index2], calipers[i].direction);
                K::Line_2 line2 = K::Line_2(points[calipers[(i + 1) % 4].index2], calipers[(i + 1) % 4].direction);

                //TODO: possono essere allineate le 2 linee
                K::Point_2 p;
                CGAL::assign(p, CGAL::intersection(line1, line2));
                tmpRect.push_back(p);
            }

            double baseLen = CGAL::approximate_sqrt((tmpRect.vertices()[1] - tmpRect.vertices()[0]).squared_length());
            double heightLen = CGAL::approximate_sqrt((tmpRect.vertices()[2] - tmpRect.vertices()[1]).squared_length());

            double tmpArea = baseLen * heightLen; //tmpRect.area();
            if (tmpArea < minArea) {
                minArea = tmpArea;
                minRect = tmpRect;
            }
            //std::cout << "Current minArea: " << minArea << std::endl;
        }
        //std::cout << "min area determined---------------------------------------------------------------" << std::endl;
        return minRect;
    }
       
    CGAL::Polyhedron_3<K> computeApproxMVBB(Points3 pts, double epsilon, double gridSize) {
        //Computa axis-aligned bounding box dei punti
        K::Iso_cuboid_3 aabb1 = CGAL::bounding_box(pts.begin(), pts.end());

        CGAL::Surface_mesh<K::Point_3> convexHull;
        CGAL::convex_hull_3(pts.begin(), pts.end(), convexHull);
        Points3 points;
        points.reserve(convexHull.number_of_vertices());
        for (const K::Point_3& p : convexHull.points()) {
            points.push_back(p);
        }


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

        //Computa secondo diametro su proiezione
        K::Segment_2 diam2 = grid2.getDiameter(snappedPts2D.cbegin(), snappedPts2D.cbegin(), snappedPts2D.cend());

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

        _LOG("Initial approx Volume = " + std::to_string(CGAL::approximate_sqrt(scaledDiam.squared_length()) * CGAL::approximate_sqrt(scaledDiam2.squared_length()) * CGAL::approximate_sqrt(scaledDiam3.squared_length())));

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
            CGAL::convex_hull_2(resultProj._projPoints.begin(), resultProj._projPoints.end(), std::back_inserter(outp));
            //CGAL::ch_graham_andrew(resultProj._projPoints.begin(), resultProj._projPoints.end(), std::back_inserter(outp));

            CGAL::Polygon_2<K> minAreaRect, minAreaRectTest;
            //CGAL::min_rectangle_2(outp.begin(), outp.end(), std::back_inserter(minAreaRectTest));
            minAreaRect = computeMABB(outp);

            K::Segment_2 base2 = K::Segment_2(minAreaRect.vertices()[0], minAreaRect.vertices()[1]);
            K::Segment_2 height2 = K::Segment_2(minAreaRect.vertices()[1], minAreaRect.vertices()[2]);
            Points3 tmpUnpMinAreaPts = candPlane.unprojectTo3D(minAreaRect.vertices_begin(), minAreaRect.vertices_end());

            K::Segment_3 base = K::Segment_3(tmpUnpMinAreaPts[0], tmpUnpMinAreaPts[1]);
            K::Segment_3 height = K::Segment_3(tmpUnpMinAreaPts[1], tmpUnpMinAreaPts[2]);

            double tmpVol = CGAL::approximate_sqrt(base.squared_length()) * CGAL::approximate_sqrt(height.squared_length()) * resultProj._diam;
            _LOG("Volume = " + std::to_string(tmpVol));
            //std::cout << "Volume = " << tmpVol << std::endl;

            if (tmpVol < bestVolume) {
                bestVolume = tmpVol;
                bestCandidate = *it;
                unpMinAreaPts = tmpUnpMinAreaPts;
                offsetPlane = resultProj._offsetPlane;
                finalDiam = resultProj._diam;
            }
        }


        K::Vector_3 kv = (bestCandidate / CGAL::approximate_sqrt(bestCandidate.squared_length())) * finalDiam;
        _LOG("Best Volume is = " + std::to_string(bestVolume));

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
}

