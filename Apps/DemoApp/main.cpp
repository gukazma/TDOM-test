#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Projection_traits_xy_3.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/compute_average_spacing.h>

#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Surface_mesh.h>

//#include <CGAL/Polygon_mesh_processing/border.h>
//#include <CGAL/Polygon_mesh_processing/remesh.h>
//#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/polygon_mesh_processing.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <boost/graph/adjacency_list.hpp>
#include <tiny_obj_loader.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <fstream>
#include <queue>
namespace SMS = CGAL::Surface_mesh_simplification;


///////////////////////////////////////////////////////////////////
//! [TIN DS]

using Kernel = CGAL::Simple_cartesian<double>;
using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;
using Segment_3 = Kernel::Segment_3;
using Vector_3 = Kernel::Vector_3;

int main(int argc, char **argv)
{
    using Mesh = CGAL::Surface_mesh<Point_3>;
    Mesh mesh;

    CGAL::Polygon_mesh_processing::IO::read_polygon_mesh("D:/Tile_+000_+000.obj", mesh);


    //std::cout << mesh.number_of_faces() << std::endl;
    //std::cout << mesh.number_of_vertices() << std::endl;
    auto face_normals = mesh.add_property_map<Mesh::Face_index, Vector_3>("f:normals", CGAL::NULL_VECTOR).first;
    //CGAL::Polygon_mesh_processing::calculate_face_normals
    CGAL::Polygon_mesh_processing::compute_face_normals(mesh, face_normals);
    auto face_areas = mesh.add_property_map<Mesh::Face_index, double>("f:areas", 0.0).first;
    CGAL::IO::write_OBJ("D:/www.obj", mesh);
    Vector_3 up = {0.0, 0.0, 1.0};
    std::vector<double> areas;
    for (auto f : mesh.faces())
    {
        auto area = CGAL::Polygon_mesh_processing::face_area(f, mesh);
        face_areas[f] = area;
        areas.push_back(area);
    }
    auto nth = areas.begin() + areas.size() / 2;
    std::nth_element(areas.begin(), nth, areas.end());
    double midienArea = *nth;

    for (auto f : mesh.faces())
    {
        Vector_3 normal = face_normals[f];
        double area = face_areas[f];
        if (std::acos(std::abs(normal * up)) * 180.0 / 3.1415926 > 80.0 &&
            std::acos(std::abs(normal * up)) * 180.0 / 3.1415926 < 90.0 &&  area>midienArea*1.2)
        {
            Mesh::halfedge_index minEdge;
            auto hf = mesh.halfedge(f);
            auto hf0 = mesh.halfedge(f);
            double minDistance = std::numeric_limits<double>::max();
            do
            {
                auto sourcePoint = mesh.point(mesh.source(hf));
                auto targetPoint = mesh.point(mesh.source(hf));
                double distance = CGAL::squared_distance(sourcePoint, targetPoint);
                if (minDistance > distance)
                {
                    minDistance = distance;
                    minEdge = hf;
                }
                
                hf = mesh.next(hf);
            } while (hf!= hf0);
            if (CGAL::Euler::does_satisfy_link_condition(mesh.edge(minEdge), mesh))
            {
                CGAL::Euler::collapse_edge(mesh.edge(minEdge), mesh);
            }
        }
    }

    //
    ////mesh.collect_garbage();
    CGAL::IO::write_OBJ("D:/www2.obj", mesh);

    return EXIT_SUCCESS;
}
