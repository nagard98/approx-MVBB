#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Aff_transformation_3<K> Transform3;
typedef std::vector<K::Point_2> Points2;
typedef std::vector<K::Point_3> Points3;

typedef struct Proj2DRes {
    Proj2DRes(Points2 projPoints, K::Vector_3 offsetPlane, double diam){
        _projPoints = projPoints;
        _offsetPlane = offsetPlane;
        _diam = diam;
    }

    Points2 _projPoints;
    K::Vector_3 _offsetPlane;
    double _diam;

} Proj2DRes;

class ProjectionPlane {

    public:
        ProjectionPlane(K::Point_3 origin, K::Vector_3 orth_vec){
            _plane = K::Plane_3(origin, orth_vec);
            setBasis();
            _basisTransform = Transform3(
                _base1.x(), _base1.y(), _base1.z(),
                _base2.x(), _base2.y(), _base2.z(),
                _base3.x(), _base3.y(), _base3.z());
        
            _inverseBasisTransform = _basisTransform.inverse();
        }

        Proj2DRes projectTo2D(Points3::iterator it, Points3::iterator end){
            bool hasNegSide = false;
            bool hasPosSide = false;

            double tmpMinDistNeg = 100000.0f;
            double tmpMinDistPos = 100000.0f;
            double tmpMaxDistNeg = 0.0f;
            double tmpMaxDistPos = 0.0f;
            K::Vector_3 tmpOffsetPlane;
            Points2 projectedPoints;

            for(; it != end; it++){
                K::Point_3 tmpProj = _plane.projection(*it);
                double dist = ((*it)-tmpProj).squared_length();

                if(_plane.has_on_negative_side(*it)){
                    hasNegSide = true;
                    if(!hasPosSide && dist < tmpMinDistNeg) tmpMinDistNeg = dist;
                    if(dist > tmpMaxDistNeg){
                        tmpOffsetPlane = K::Vector_3(tmpProj, *it);
                        tmpMaxDistNeg = dist;
                    }
                }else{
                    hasPosSide = true;
                    if(!hasNegSide && dist < tmpMinDistPos){
                        tmpMinDistPos = dist;
                        tmpOffsetPlane = K::Vector_3(tmpProj, *it);
                    }
                    if(dist > tmpMaxDistPos) {
                        tmpMaxDistPos = dist;
                    }
                }

                K::Point_3 pointInPlaneBas = tmpProj.transform(_basisTransform);

                projectedPoints.push_back(K::Point_2(pointInPlaneBas.x(), pointInPlaneBas.y()));
            }

            return Proj2DRes(projectedPoints, tmpOffsetPlane, CGAL::approximate_sqrt(tmpMaxDistNeg) + CGAL::approximate_sqrt(tmpMaxDistPos));
        }

        Points3 unprojectTo3D(CGAL::Polygon_2<K>::iterator it, CGAL::Polygon_2<K>::iterator end){
            Points3 tmpUnpPoints;
            for(; it != end; it++){
                K::Point_3 tmp = (K::Point_3(it->x(), it->y(), 0)).transform(_inverseBasisTransform);
                tmpUnpPoints.push_back(tmp);
            }
            
            return tmpUnpPoints;
        }

        K::Point_3 unprojectTo3D(K::Point_2 point){
            return (K::Point_3(point.x(), point.y(), 0)).transform(_inverseBasisTransform);
        }


    private:
        void setBasis(){
            _base1 = _plane.base1();
            _base1 = _base1 / CGAL::approximate_sqrt(_base1.squared_length());
            _base2 = _plane.base2();
            _base2 = _base2 / CGAL::approximate_sqrt(_base2.squared_length());
            _base3 = _plane.orthogonal_vector();
            _base3 = _base3 / CGAL::approximate_sqrt(_base3.squared_length());
        };

        K::Plane_3 _plane;
        K::Vector_3 _base1, _base2, _base3;
        K::Aff_transformation_3 _basisTransform, _inverseBasisTransform;

};