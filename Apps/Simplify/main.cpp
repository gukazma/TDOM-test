#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Simple_cartesian.h>

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

// #include <CGAL/Polygon_mesh_processing/border.h>
// #include <CGAL/Polygon_mesh_processing/remesh.h>
// #include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/polygon_mesh_processing.h>
#include <boost/graph/adjacency_list.hpp>
#include <fstream>
#include <queue>
#include <tiny_obj_loader.h>
#include <filesystem>
// Visitor base
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_distance_placement.h>
double minlength;
double maxLength = std::numeric_limits<double>::min();
namespace SMS = CGAL::Surface_mesh_simplification;
using Kernel = CGAL::Simple_cartesian<double>;
using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;
using Segment_3 = Kernel::Segment_3;
using Vector_3 = Kernel::Vector_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
Vector_3 up = {0.0, 0.0, 1.0};
std::vector<Mesh::face_index> face_tomerge;
struct Stats
{
    std::size_t collected = 0;
    std::size_t processed = 0;
    std::size_t collapsed = 0;
    std::size_t non_collapsable = 0;
    std::size_t cost_uncomputable = 0;
    std::size_t placement_uncomputable = 0;
};
template <class TM>
class Edge_custom_cost
{
  public:
    const TM *sm_ptr;

    Edge_custom_cost()
    {
    }
    typedef SMS::Edge_profile<Mesh> Profile;

    std::optional<typename Profile::FT> operator()(const Profile &profile,
                                                   const std::optional<Mesh::Point>& /*placement*/) const
    {

        typedef std::optional<typename Profile::FT> result_type;
        typename Profile::FT rnt =
            profile.geom_traits().compute_squared_distance_3_object()(profile.p0(), profile.p1());
        Mesh::halfedge_index hf0 = profile.v0_v1();
        Mesh::halfedge_index hf1 = profile.v1_v0();
        auto face0 = sm_ptr->face(hf0);
        auto face1 = sm_ptr->face(hf1);
        if (!is_valid_face_descriptor(face0, *sm_ptr) || !is_valid_face_descriptor(face1, *sm_ptr))
        {
            rnt = 100.0;
            return result_type(rnt);
        }
        
        Vector_3 normal0, normal1;
        normal0 = CGAL::Polygon_mesh_processing::compute_face_normal(face0, *sm_ptr);
        normal1 = CGAL::Polygon_mesh_processing::compute_face_normal(face1, *sm_ptr);

        auto fbox0 = CGAL::Polygon_mesh_processing::face_bbox(face0, *sm_ptr);
        auto fbox1 = CGAL::Polygon_mesh_processing::face_bbox(face1, *sm_ptr);
        double face_lengh0 = std::abs(fbox0.zmax() - fbox0.zmin());
        double face_lengh1 = std::abs(fbox1.zmax() - fbox1.zmin());

        double degree0 = std::acos(std::abs(normal0 * up)) * 180.0 / 3.1415926;
        double degree1 = std::acos(std::abs(normal1 * up)) * 180.0 / 3.1415926;
        if ((degree0 > 60.0 && degree0 < 120.0 && face_lengh0 > maxLength / 5.0) ||
            (degree1 > 60.0 && degree1 < 120.0 && face_lengh1 > maxLength / 5.0))
        {
            /*if (face_lengh0 > 2.0 && face_lengh1 > 2.0)
            {
                rnt = 100.0;
            }*/
            return result_type(rnt);
        }

        rnt = 100.0;
        return result_type(rnt);
    }
};

class Edge_count_stop_custom
{
  public:
    const Mesh *sm_ptr;

    Edge_count_stop_custom()
    {
    }

    template <typename F, typename Profile>
    bool operator()(const F & current_cost, const Profile & profile, std::size_t /*initial_edge_count*/,
                    std::size_t current_edge_count) const
    {
        return current_cost > minlength;
    }

  private:
    std::size_t m_edge_count_threshold;
};
class My_placement
{
  public:

    My_placement()
    {
    }
    typedef SMS::Edge_profile<Mesh> Profile;
    std::optional<Profile::Point> operator()(const Profile &profile) const
    {
        typedef std::optional<typename Profile::Point> result_type;
        auto& mesh = profile.surface_mesh();
        
        auto face0 = mesh.face(profile.v0_v1());
        auto face1 = mesh.face(profile.v1_v0());
        if (!is_valid_face_descriptor(face0, mesh) || !is_valid_face_descriptor(face1, mesh))
        {
            return result_type(profile.geom_traits().construct_midpoint_3_object()(profile.p0(), profile.p1()));
        }
        Vector_3 normal0, normal1;
        normal0 = CGAL::Polygon_mesh_processing::compute_face_normal(face0, mesh);
        normal1 = CGAL::Polygon_mesh_processing::compute_face_normal(face1, mesh);

        double degree0 = std::acos(std::abs(normal0 * up)) * 180.0 / 3.1415926;
        double degree1 = std::acos(std::abs(normal1 * up)) * 180.0 / 3.1415926;
        //Mesh::face_index targetFace;
        Vector_3 targetNormal;
        if (degree0 > 60.0 && degree0 < 120.0)
        {
            targetNormal = normal0;
        }
        else
        {
            targetNormal = normal1;
        }
        auto midpoint = profile.geom_traits().construct_midpoint_3_object()(profile.p0(), profile.p1());
        midpoint += targetNormal*0.1;
        //Profile::Point rnt = halfedgesV0.size() > halfedgesV1.size() ? profile.p1() : profile.p0();

        return result_type(midpoint);
    }
};
// BGL property map which indicates whether an edge is marked as non-removable
struct Border_is_constrained_edge_map
{
    const Mesh *sm_ptr;
    typedef boost::graph_traits<Mesh>::edge_descriptor key_type;
    typedef bool value_type;
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;

    Border_is_constrained_edge_map(const Mesh &sm) : sm_ptr(&sm)
    {
    }

    friend value_type get(const Border_is_constrained_edge_map &m, const key_type &edge)
    {
        auto hf = m.sm_ptr->halfedge(edge);
        auto face0 = m.sm_ptr->face(hf);
        auto face1 = m.sm_ptr->face(m.sm_ptr->opposite(hf));
        if (std::find(face_tomerge.begin(), face_tomerge.end(), face0) != face_tomerge.end());
        {
            return false;
        }
        if (std::find(face_tomerge.begin(), face_tomerge.end(), face1) != face_tomerge.end())
            ;
        {
            return false;
        }
        return true;
    }
};

struct My_visitor : SMS::Edge_collapse_visitor_base<Mesh>
{
    My_visitor(Stats *s) : stats(s)
    {
    }

    // Called during the collecting phase for each edge collected.
    void OnCollected(const Profile &, const std::optional<double> &)
    {
        ++(stats->collected);
        std::cerr << "\rEdges collected: " << stats->collected << std::flush;
    }

    // Called during the processing phase for each edge selected.
    // If cost is absent the edge won't be collapsed.
    void OnSelected(const Profile &, std::optional<double> cost, std::size_t initial, std::size_t current)
    {
        ++(stats->processed);
        if (!cost)
            ++(stats->cost_uncomputable);

        if (current == initial)
            std::cerr << "\n" << std::flush;
        std::cerr << "\r" << current << std::flush;
    }

    // Called during the processing phase for each edge being collapsed.
    // If placement is absent the edge is left uncollapsed.
    void OnCollapsing(const Profile &, std::optional<Point> placement)
    {
        if (!placement)
            ++(stats->placement_uncomputable);
    }

    // Called for each edge which failed the so called link-condition,
    // that is, which cannot be collapsed because doing so would
    // turn the surface mesh into a non-manifold.
    void OnNonCollapsable(const Profile &)
    {
        ++(stats->non_collapsable);
    }

    // Called after each edge has been collapsed
    void OnCollapsed(const Profile &, vertex_descriptor)
    {
        ++(stats->collapsed);
        // CGAL::IO::write_OBJ("D:/www2.obj", *sm_ptr);
    }

    Stats *stats;
    const Mesh *sm_ptr;
};

// Placement class
typedef SMS::Constrained_placement<SMS::Midpoint_placement<Mesh>, Border_is_constrained_edge_map> Placement;

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cout << "Simplify.exe arg1 arg2";
        std::cout << "arg1: root dir";
        std::cout << "arg2: min length";
        return 0;
    }
    std::filesystem::path rootDir = argv[1];
    minlength = std::stod(argv[2]);
    for (auto& fe : std::filesystem::directory_iterator(rootDir))
    {
        if (std::filesystem::is_directory(fe.path()))
        {
            continue;
        }
        Mesh mesh;
        auto fp = fe.path();
        std::string filename = fp.generic_string();
        std::cout << "simplify " << filename << std::endl;
        CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, mesh);
        // CGAL::Polygon_mesh_processing::IO::read_polygon_mesh("D:/cube.obj", mesh);

        // std::cout << mesh.number_of_faces() << std::endl;
        // std::cout << mesh.number_of_vertices() << std::endl;
        auto face_normals = mesh.add_property_map<Mesh::Face_index, Vector_3>("f:normals", CGAL::NULL_VECTOR).first;
        // CGAL::Polygon_mesh_processing::calculate_face_normals
        CGAL::Polygon_mesh_processing::compute_face_normals(mesh, face_normals);
        auto face_lenghs = mesh.add_property_map<Mesh::Face_index, double>("f:lengths", 0.0).first;
        std::vector<double> lengths;
        lengths.reserve(mesh.number_of_faces());
        // CGAL::IO::write_OBJ("D:/www.obj", mesh);
        Vector_3 up = {0.0, 0.0, 1.0};
        double maxLength = std::numeric_limits<double>::min();
        for (auto f : mesh.faces())
        {
            auto fbox = CGAL::Polygon_mesh_processing::face_bbox(f, mesh);

            double face_lengh = std::abs(fbox.zmax() - fbox.zmin());
            face_lenghs[f] = face_lengh;
            lengths.push_back(face_lenghs[f]);
            if (maxLength < face_lengh)
            {
                maxLength = face_lengh;
            }
        }

        Mesh::Halfedge_index hf = *mesh.halfedges_begin();

        for (auto f : mesh.faces())
        {
            Vector_3 normal = face_normals[f];
            auto length = face_lenghs[f];
            double degree = std::acos(std::abs(normal * up)) * 180.0 / 3.1415926;
            if (degree > 60.0 && degree < 120.0 && length > maxLength / 5.0)
            {
                face_tomerge.push_back(f);
            }
            else
            {
            }
        }

        Edge_count_stop_custom stop;
        //SMS::Edge_length_stop_predicate<double> stop(minlength);
        // SMS::Edge_length_cost<Mesh>();
        Border_is_constrained_edge_map bem(mesh);
        Stats stats;
        My_visitor vis(&stats);
        vis.sm_ptr = &mesh;
        // This the actual call to the simplification algorithm.
        // The surface mesh and stop conditions are mandatory arguments.
        // The index maps are needed because the vertices and edges
        // of this surface mesh lack an "id()" field.
        SMS::Edge_length_cost<Mesh>();
        Edge_custom_cost<Mesh> cost;
        cost.sm_ptr = &mesh;
        /*int r = SMS::edge_collapse(mesh, stop,
                                   CGAL::parameters::edge_is_constrained_map(bem)
                                       .get_placement(Placement(bem)));*/
        int r = SMS::edge_collapse(
            mesh, stop,
            CGAL::parameters::get_cost(cost).get_placement(/*SMS::Midpoint_placement<Mesh>()*/ My_placement()));

        mesh.collect_garbage();
        std::filesystem::path outputFilename = filename;
        std::filesystem::path outputDir = outputFilename.parent_path() / "output";
        if (!std::filesystem::exists(outputDir))
        {
            std::filesystem::create_directories(outputDir);
        }
        outputFilename = outputDir / outputFilename.filename();

        CGAL::IO::write_OBJ(outputFilename.generic_string(), mesh);
        std::cout << "complete" << std::endl;
    }

    return EXIT_SUCCESS;
}
