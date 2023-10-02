#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <unordered_set>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef std::vector<K::Point_3> Points3;
typedef std::unordered_set<K::Point_3> GridPoints3;

typedef std::vector<K::Point_2> Points2;
typedef std::unordered_set<K::Point_2> GridPoints2;

typedef struct Coords2 {
    float x, y;
} Coords2;

typedef struct BoxDirections3 {
    BoxDirections3(){}
    BoxDirections3(K::Vector_3 x, K::Vector_3 y, K::Vector_3 z, K::Point_3 o){
        _x = x;
        _y = y;
        _z = z;
        _origin = o;
        _lenX = CGAL::approximate_sqrt(_x.squared_length());
        _lenY = CGAL::approximate_sqrt(_y.squared_length());
        _lenZ = CGAL::approximate_sqrt(_z.squared_length());
    }

    K::Vector_3 _x, _y, _z;
    K::Point_3 _origin;
    double _lenX, _lenY, _lenZ;

} BoxDirections3;

typedef struct BoxDirections2 {

    BoxDirections2(){}
    BoxDirections2(K::Vector_2 x, K::Vector_2 y, K::Point_2 o){
        _x = x;
        _y = y;
        _origin = o;
        _lenX = CGAL::approximate_sqrt(_x.squared_length());
        _lenY = CGAL::approximate_sqrt(_y.squared_length());
    }

    K::Vector_2 _x, _y;
    K::Point_2 _origin;
    double _lenX, _lenY;

} BoxDirections2;

class Grid3{

    public:
        Grid3(BoxDirections3 boxDirections){
            this->boxDirs = boxDirections;
        }

        GridPoints3 snapToGrid(Points3::iterator it, Points3::iterator end){
            GridPoints3 snappedPoints;
            GridPoints3 rescaledPoints;

            for(; it != end; it++){
                K::Vector_3 pointVector = K::Vector_3(boxDirs._origin, *it);
                double snappedX = closestInteger(((pointVector * boxDirs._x) / boxDirs._lenX) / boxDirs._lenX);
                double snappedY = closestInteger(((pointVector * boxDirs._y) / boxDirs._lenY) / boxDirs._lenY);
                double snappedZ = closestInteger(((pointVector * boxDirs._z) / boxDirs._lenZ) / boxDirs._lenZ);

                K::Point_3 gridCoords = K::Point_3(snappedX, snappedY, snappedZ);
                int found = snappedPoints.count(gridCoords);

                if (found == 0) {
                    snappedPoints.insert(gridCoords);
                    rescaledPoints.insert(K::Point_3(snappedX * boxDirs._lenX, snappedY * boxDirs._lenY, snappedY * boxDirs._lenZ));
                }
            }

            return rescaledPoints;
        }


        std::vector<K::Vector_3> getSearchCandidates(double distance){
            std::vector<K::Vector_3> tmpPoints;
            for(double i=0.0f; i <= distance; i += 1.0f ){
                for(double j=0.0f; j <= distance; j += 1.0f ){
                    for(double k=0.0f; k <= distance; k += 1.0f ){
                        K::Vector_3 tmpPoint = i * boxDirs._x + j * boxDirs._y + k * boxDirs._z;
                        if(chebyshevDistance(tmpPoint) <= distance){
                            if(i==0.0f && j==0.0f && k==0.0f) break;
                            else tmpPoints.push_back(tmpPoint);
                        }                        
                    }
                }
            }  
/*            double negDistance = -distance;
            for(double i=-1.0f; i >= negDistance ; i -= 1.0f ){
                for(double j=-1.0f; j >= negDistance; j -= 1.0f ){
                    for(double k=-1.0f; k >= negDistance; k -= 1.0f ){
                        K::Vector_3 tmpPoint = i * boxDirs._x + j * boxDirs._y + k * boxDirs._z;
                        if(chebyshevDistance(tmpPoint) <= distance){
                            tmpPoints.push_back(tmpPoint);
                        }
                    }
                }
            } */

            return tmpPoints;
        }

        std::pair<K::Point_3, K::Point_3> getDiameter(GridPoints3::iterator it1, GridPoints3::iterator it2, GridPoints3::iterator end){
            double tmpMaxLength = 0;
            K::Point_3 first, second;

            for(; it1 != end; it1++){
                for(it2 = it1; it2 != end; it2++){
                    double tmpLength = (*it1 - *it2).squared_length();
                    if(tmpMaxLength < tmpLength){
                        tmpMaxLength = tmpLength;
                        first = *it1;
                        second = *it2;
                    }
                }
            }

            return std::pair<K::Point_3,K::Point_3>(first, second);
        }


    private:
        BoxDirections3 boxDirs;

        double chebyshevDistance(K::Vector_3 a){
            double xDist = CGAL::abs(a.x());
            double yDist = CGAL::abs(a.y());
            double zDist = CGAL::abs(a.z());

            return CGAL::max(CGAL::max(xDist, yDist), zDist);
        }

        double closestInteger(double val){
            if(val < 0.0f){
                return (double)((int) (val - 0.5f));
            }else{
                return (double)((int) (val + 0.5f));
            }
        }


};

class Grid2{

    public:
        Grid2(BoxDirections2 boxDirections){
            this->boxDirs = boxDirections;
        }

        GridPoints2 snapToGrid(Points2::iterator it, Points2::iterator end){
            GridPoints2 snappedPoints;
            GridPoints2 rescaledPoints;

            for(; it != end; it++){
                K::Vector_2 pointVector = K::Vector_2(boxDirs._origin, *it);
                double snappedX = closestInteger(((pointVector * (boxDirs._x / boxDirs._lenX)) / boxDirs._lenX));
                double snappedY = closestInteger(((pointVector * (boxDirs._y / boxDirs._lenY)) / boxDirs._lenY));

                K::Point_2 gridCoords = K::Point_2(snappedX, snappedY);
                int found = snappedPoints.count(gridCoords);

                if (found == 0) {
                    snappedPoints.insert(gridCoords);
                    rescaledPoints.insert(K::Point_2(snappedX * boxDirs._lenX, snappedY * boxDirs._lenY));
                }
            }

            return rescaledPoints;
        }

        K::Segment_2 getDiameter(GridPoints2 snappedPoints){
            double tmpMaxLength = 0;
            K::Segment_2 tmpMaxDiameter;

            for(GridPoints2::iterator it1 = snappedPoints.begin(); it1 != snappedPoints.end(); it1++){
                for (GridPoints2::iterator it2 = it1; it2 != snappedPoints.end(); it2++) {
                    K::Segment_2 tmpDiameter = K::Segment_2(*it1, *it2);
                    double tmpLength = tmpDiameter.squared_length();
                    if(tmpMaxLength < tmpLength){
                        tmpMaxLength = tmpLength;
                        tmpMaxDiameter = tmpDiameter;
                    }
                }
            }

            return tmpMaxDiameter;
        }

        Points2 projectToLine(GridPoints2::iterator it, GridPoints2::iterator end, K::Line_2 hLine){
            Points2 projectedPoints;

            for(; it != end; it++){
                K::Point_2 proj = hLine.projection(*it);
                projectedPoints.push_back(proj);
            }

            return projectedPoints;
        }


    private:
        BoxDirections2 boxDirs;

        double closestInteger(double val){
            if(val < 0.0f){
                return (double)((int) (val - 0.5f));
            }else{
                return (double)((int) (val + 0.5f));
            }
        }


};