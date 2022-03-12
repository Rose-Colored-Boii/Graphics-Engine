//
// Created by Liam on 9/03/2022.
//

#ifndef ENGINE_LSYSTEMELEMENT_H
#define ENGINE_LSYSTEMELEMENT_H
#include "Line2D.h"

class LSystemElement {
private:
    Point2D position;
    double angle;
public:
    LSystemElement(Point2D position, double angle){
        this->position = position;
        this->angle = angle;
    }

    const Point2D &getPosition() const {
        return position;
    }

    void setPosition(const Point2D &position) {
        LSystemElement::position = position;
    }

    double getAngle() const {
        return angle;
    }

    void setAngle(double angle) {
        LSystemElement::angle = angle;
    }
};


#endif //ENGINE_LSYSTEMELEMENT_H
