//
// Created by liam on 4/29/22.
//

#ifndef ENGINE_CC_LIGHT_H
#define ENGINE_CC_LIGHT_H

#include "easy_image.h"
#include "vector3d.h"

class Light {
public:
    std::vector<double> ambientLight = {1, 1, 1};
    std::vector<double> diffuseLight = {0, 0, 0};
    std::vector<double> specularLight = {0, 0, 0};
    Vector3D direction;
    bool infinity = false;
    double angle = 90;
    Vector3D location;
};


#endif //ENGINE_CC_LIGHT_H
