#pragma once

#include "gx_geometry.h"


#include <set>
#include <algorithm>
#include <cassert>
#include <numeric>

static void print_vector(std::vector<int> v) {
    for (auto x : v) {
        std::cout << x << " ";
    }
    std::cout<<"\n";
}

void writeOBJ(const char* obj, const std::vector<Vector> &vertices,
                const std::vector<TriangleIndices> &indices) {
    FILE* f = fopen(obj, "w+");
    for (int i = 0; i < vertices.size(); i++) {
        fprintf(f, "v %f %f %f\n", vertices[i][0], vertices[i][1], vertices[i][2]);
    }
    for (int i = 0; i < indices.size(); i++) {
        fprintf(f, "f %u %u %u\n", indices[i].vtxi + 1, indices[i].vtxj + 1, indices[i].vtxk + 1);
    }
    fclose(f);
}


std::vector< std::vector<int> > NeighborsList(const TriangleMesh &mesh) {
    const int N = mesh.vertices.size();

    std::vector< std::vector<int> > neighbors(N);

    for (int i = 0; i < mesh.indices.size(); i++) {
        const auto &index_struct = mesh.indices[i];
        int vtxs[3] = {index_struct.vtxi, index_struct.vtxj, index_struct.vtxk};
        for (auto vtx1 : vtxs) {
            for (auto vtx2 : vtxs) {
                if (vtx1 <= vtx2) { continue; }
                neighbors[vtx1].push_back(vtx2);
                neighbors[vtx2].push_back(vtx1);
            }
        }
    }

    return neighbors;
}

/// Returns the indices of the ordered boundary vertices 
std::vector<int> BoundaryVertices(const TriangleMesh &mesh) {
    
    const int N = mesh.vertices.size();
    const auto neighbors = NeighborsList(mesh);

    // An internal neighbor of vis a vertex that appears more than once in the
    // neighbor list neighbors[v]. The next function removes all the internal
    // neighbors.
    auto remove_internal_neighbors = [](std::vector<int> values){
        
        std::sort(values.begin(), values.end());
        std::vector<int> output; 

        const int n = values.size();
        for(int ptr = 0; ptr < n; ptr++) {
            // Loop invariant: values[ptr] != values[ptr-1] at this current point.
            if (ptr == n - 1 || values[ptr] != values[ptr + 1]) {
                output.push_back(values[ptr]);
            }
            else {
                while(ptr < n - 1 && values[ptr] == values[ptr + 1]) { ptr++; }
            }
        }

        return output;
    };

    std::vector< std::vector<int> > internal_neighbors(N);
    std::transform(neighbors.begin(), neighbors.end(), internal_neighbors.begin(),
                    remove_internal_neighbors);

    std::vector<int> output;
    std::set<int> used;

    // first vertex in the boundary
    const auto ptr_first_vtx = std::find_if(internal_neighbors.begin(),
                                            internal_neighbors.end(),
                                            [](const auto &v){ return !v.empty(); });

    assert(ptr_first_vtx != internal_neighbors.end());
    const int first_vtx = std::distance(internal_neighbors.begin(), ptr_first_vtx);

    output.push_back(first_vtx);
    used.insert(first_vtx);

    for(;;) {
        const int curr = output.back();

        const auto next_ptr = std::find_if(internal_neighbors[curr].begin(),
            internal_neighbors[curr].end(), [&used](auto vtx){ return !used.contains(vtx); });
        
        if (next_ptr == internal_neighbors[curr].end()) {
            break;
        }

        const int next = *next_ptr;
        output.push_back(next);
        used.insert(next);
    }

    // After all this, we have output representing the vertices in the boundary
    // in order.

    return output;
}

void TutteEmbedding(const TriangleMesh &mesh, const char *obj) {

    const auto boundary = BoundaryVertices(mesh);
    const int n = boundary.size();

    // Set of boundary points for fast search
    std::set<int> is_in_boundary;
    is_in_boundary.insert(boundary.begin(), boundary.end());

    std::cout << "Boundary Size: " << is_in_boundary.size() << std::endl;

    assert(n >= 3);

    const auto compute_boundary_length = [&]() {
        double res = (mesh.vertices[boundary[0]] - mesh.vertices[boundary[n-1]]).norm();
        for (int i = 1; i < n; i++) {
            res += (mesh.vertices[boundary[i]] - mesh.vertices[boundary[i-1]]).norm();
        }
        return res;
    };

    const double boundary_length = compute_boundary_length();
    

    // Preparing new parametrization
    const int N = mesh.vertices.size();
    std::vector<Vector> vertices(N);

    double cs = 0;
    for (int i = 0; i < n; i++) {
        auto idx = boundary[i];
        const double theta = 2 * M_PI * cs / boundary_length;
        vertices[idx] = Vector(cos(theta), sin(theta), 0.);
        if (i < n - 1) {
            cs += (mesh.vertices[boundary[i+1]] - mesh.vertices[boundary[i]]).norm();
        }
    }

    // To change
    const int nb_iter = 1000;
    const auto neighbors = NeighborsList(mesh);

    const auto update_vertices = [N, &is_in_boundary, &neighbors](std::vector<Vector> vertices_t) {
        auto vertices_tp1 = vertices_t;

        for (int idx = 0; idx < N; idx++) {
            if (!is_in_boundary.contains(idx)) {
                vertices_tp1[idx] = std::accumulate(neighbors[idx].begin(),
                    neighbors[idx].end(), Vector(), [&](Vector acc, int other){
                        return acc + vertices_t[other];
                    }) / (double) neighbors[idx].size();
            }
        }

        return vertices_tp1;
    };

    for (int iter = 0; iter < nb_iter; iter++) {
        vertices = update_vertices(vertices);
    }

    writeOBJ(obj, vertices, mesh.indices);
}