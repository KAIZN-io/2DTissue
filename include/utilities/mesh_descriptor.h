// mesh_descriptor.h

#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>

using Kernel = CGAL::Simple_cartesian<double>;
using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;

namespace _3D {
    using Mesh = CGAL::Surface_mesh<Point_3>;
    using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
    using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
    using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;
    using Seam_edge_pmap = Mesh::Property_map<edge_descriptor, bool>;
    using Seam_vertex_pmap = Mesh::Property_map<vertex_descriptor, bool>;
    using UV_pmap = Mesh::Property_map<halfedge_descriptor, Point_2>;
}

namespace UV {
    using Mesh = CGAL::Seam_mesh<_3D::Mesh,
                                _3D::Seam_edge_pmap,
                                _3D::Seam_vertex_pmap>;
    using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
    using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
}