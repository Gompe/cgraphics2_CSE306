#pragma once

#include "gx_polygon.h"

#define AUTO_SHIFT_SCALE 1

static std::vector<Polygon> scale_polygons(const std::vector<Polygon> &polygons) {
    double coord_max = std::numeric_limits<double>::lowest();
    double coord_min = std::numeric_limits<double>::max();

    for (const auto &polygon : polygons) {
        for (const auto &P : polygon.vertices) {
            coord_max = std::max(coord_max, std::max(P[0], P[1]));
            coord_min = std::min(coord_min, std::min(P[0], P[1]));
        }
    }

    coord_max = coord_max;
    coord_min = coord_min;

    auto shiftScale = [coord_max, coord_min](const Vector &P) {
        const Vector onesVector = Vector(1, 1, 0); 
        return (P - onesVector * coord_min) / (coord_max - coord_min);
    };

    std::vector<Polygon> scaled_polygons;
    scaled_polygons.reserve(polygons.size());

    for (const auto &polygon : polygons) {
        Polygon output_polygon;
        output_polygon.vertices.reserve(polygon.size());

        for (const auto &P : polygon.vertices) {
            output_polygon.addVertex(shiftScale(P));
        }
        scaled_polygons.push_back(output_polygon);
    }

    return scaled_polygons;
}

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &_polygons, std::string filename, std::string fillcol = "none") {
#if AUTO_SHIFT_SCALE
    const std::vector<Polygon> &polygons = scale_polygons(_polygons);
#else
    const std::vector<Polygon> &polygons = _polygons;
#endif
    FILE* f = fopen(filename.c_str(), "w+"); 
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (int i=0; i<polygons.size(); i++) {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \""); 
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}


// Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
// polygons is a list of polygons, describing the current frame.
// The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg_animated(const std::vector<Polygon> &_polygons, std::string filename, int frameid, int nbframes) {
#if AUTO_SHIFT_SCALE
    const std::vector<Polygon> &polygons = scale_polygons(_polygons);
#else
    const std::vector<Polygon> &polygons = _polygons;
#endif
    FILE* f;
    if (frameid == 0) {
        f = fopen(filename.c_str(), "w+");
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        fprintf(f, "<g>\n");
    } else {
        f = fopen(filename.c_str(), "a+");
    }
    fprintf(f, "<g>\n");
    for (int i = 0; i < polygons.size(); i++) {
        fprintf(f, "<polygon points = \""); 
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000-polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
    }
    fprintf(f, "<animate\n");
    fprintf(f, "    id = \"frame%u\"\n", frameid);
    fprintf(f, "    attributeName = \"display\"\n");
    fprintf(f, "    values = \"");
    for (int j = 0; j < nbframes; j++) {
        if (frameid == j) {
            fprintf(f, "inline");
        } else {
            fprintf(f, "none");
        }
        fprintf(f, ";");
    }
    fprintf(f, "none\"\n    keyTimes = \"");
    for (int j = 0; j < nbframes; j++) {
        fprintf(f, "%2.3f", j / (double)(nbframes));
        fprintf(f, ";");
    }
    fprintf(f, "1\"\n   dur = \"5s\"\n");
    fprintf(f, "    begin = \"0s\"\n");
    fprintf(f, "    repeatCount = \"indefinite\"/>\n");
    fprintf(f, "</g>\n");
    if (frameid == nbframes - 1) {
        fprintf(f, "</g>\n");
        fprintf(f, "</svg>\n");
    }
    fclose(f);
}