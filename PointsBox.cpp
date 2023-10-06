#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <unordered_set>
#include <CGAL/Point_set_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/draw_triangulation_3.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Qt/PointsGraphicsItem.h>
#include <CGAL/Qt/GraphicsViewNavigation.h>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QtGui>
#include <QApplication>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

class PointsBox : public CGAL::Basic_viewer_qt {
    public:
        PointsBox(QWidget* parent, const CGAL::Point_set_3<K::Point_3>& pointCloud , CGAL::Polyhedron_3<K> box ) : CGAL::Basic_viewer_qt(parent, "title", true, true, false, false, false, true, true){
            this->_pointCloud = pointCloud;
            this->_box = box;

            computeElements();
        }

    private:
        void computeElements(){

            add_segment(K::Point_3(-10,0,0), K::Point_3(+10,0,0), CGAL::IO::Color(0,255,0));
            add_segment(K::Point_3(0,0,-10), K::Point_3(0,0,10), CGAL::IO::Color(0,255,0));
            add_segment(K::Point_3(0,-10,0), K::Point_3(0,+10,0), CGAL::IO::Color(0,255,0));

            for(CGAL::Point_set_3<K::Point_3>::iterator it = _pointCloud.begin(); it != _pointCloud.end(); it++){
                add_point(_pointCloud.point(*it));
            }
            
            for(auto it = _box.edges_begin(); it != _box.edges_end(); it++){
                add_segment(
                    (it->vertex())->point(), 
                    ((it->opposite())->vertex())->point()
                    );
            } 
            
        }

    protected:
        CGAL::Point_set_3<K::Point_3> _pointCloud;
        CGAL::Polyhedron_3<K> _box; 
};
