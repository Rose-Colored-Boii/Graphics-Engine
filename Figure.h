//
// Created by Liam on 11/03/2022.
//

#ifndef ENGINE_FIGURE_H
#define ENGINE_FIGURE_H
#include "vector"
#include "vector3d.h"
#include "easy_image.h"
using namespace std;

class Face{
public:
    vector<int> point_indexes;

    Face(vector<int> point_indexes){
        this->point_indexes = point_indexes;
    }
};

class Figure {
public:
    vector<Vector3D> points;
    vector<Face> faces;
    img::Color color;

    Figure(){};

    Figure(vector<Vector3D> points, vector<Face> faces, img::Color color){
        this->points = points;
        this->faces = faces;
        this->color = color;
    }
};


#endif //ENGINE_FIGURE_H
