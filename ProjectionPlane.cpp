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

        // proietta punti su piano cambiando anche la base una volta proiettati
        Proj2DRes project(Points3::iterator it, Points3::iterator end){
            bool hasLSide = false, hasRSide = false;
            double tmpMinDistLSide = CGAL_IA_MAX_DOUBLE, tmpMinDistRSide = CGAL_IA_MAX_DOUBLE;
            double tmpMaxDistLSide = 0.0f, tmpMaxDistRSide = 0.0f;;
            K::Vector_3 tmpOffsetPlane, projVec, planeNormal = _plane.orthogonal_vector();
            Points2 projectedPoints;

            int k1 = 0, k2 = 0;

            for(; it != end; it++){
                K::Point_3 tmpProj = _plane.projection(*it);
                projVec = (*it) - tmpProj;
                double dist = projVec.squared_length();

                if (planeNormal * projVec < 0) {
                    hasLSide = true;
                    if (!hasRSide && dist < tmpMinDistLSide) tmpMinDistLSide = dist;
                    if (dist > tmpMaxDistLSide) {
                        tmpOffsetPlane = K::Vector_3(tmpProj, *it);
                        tmpMaxDistLSide = dist;
                    }
                }
                else {
                    hasRSide = true;
                    if (!hasLSide && dist < tmpMinDistRSide) {
                        tmpMinDistRSide = dist;
                        tmpOffsetPlane = K::Vector_3(tmpProj, *it);
                    }
                    if (dist > tmpMaxDistRSide) {
                        tmpMaxDistRSide = dist;
                    }
                }

                //Cambio di base
                K::Point_3 pointInPlaneBas = tmpProj.transform(_basisTransform);

                projectedPoints.push_back(K::Point_2(pointInPlaneBas.x(), pointInPlaneBas.y()));
            }

            return Proj2DRes(projectedPoints, tmpOffsetPlane, CGAL::approximate_sqrt(tmpMaxDistLSide) + CGAL::approximate_sqrt(tmpMaxDistRSide));
        }

        // Riporta punti a base originale prima della proiezione
        Points3 restoreBasis(CGAL::Polygon_2<K>::iterator it, CGAL::Polygon_2<K>::iterator end){
            Points3 tmpPoints;
            for(; it != end; it++){
                K::Point_3 tmp = (K::Point_3(it->x(), it->y(), 0)).transform(_inverseBasisTransform);
                tmpPoints.push_back(tmp);
            }
            
            return tmpPoints;
        }

        K::Point_3 restoreBasis(K::Point_2 point) {
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