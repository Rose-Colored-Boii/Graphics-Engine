//
// Created by Liam on 11/03/2022.
//

#ifndef ENGINE_FIGURE_H
#define ENGINE_FIGURE_H

#include "vector"
#include "Texture.h"
#include "vector3d.h"
#include "easy_image.h"

using namespace std;

class Face {
public:
    vector<int> point_indexes;

    Face(vector<int> point_indexes) {
        this->point_indexes = point_indexes;
    }
};

class Figure {
public:
    Texture texture;
    vector<Vector3D> points;
    vector<Face> faces;
    img::Color color;
    vector<double> ambientReflection = {0, 0, 0};
    vector<double> diffuseReflection = {0, 0, 0};
    vector<double> specularReflection = {0, 0, 0};
    double reflectionCoefficient;

    Figure() {};

    Figure(vector<Vector3D> points, vector<Face> faces, img::Color color) {
        this->points = points;
        this->faces = faces;
        this->color = color;
    }
};


#endif //ENGINE_FIGURE_H
