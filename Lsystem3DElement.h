//
// Created by liam on 3/23/22.
//

#ifndef ENGINE_CC_LSYSTEM3DELEMENT_H
#define ENGINE_CC_LSYSTEM3DELEMENT_H
#include "vector3d.h"

class Lsystem3DElement {
public:
    Vector3D position;
    Vector3D H;
    Vector3D L;
    Vector3D U;
    Lsystem3DElement(Vector3D& position, Vector3D& H, Vector3D& L, Vector3D& U){
        this->position = position;
        this->H = H;
        this->L = L;
        this->U = U;
    }
};


#endif //ENGINE_CC_LSYSTEM3DELEMENT_H
