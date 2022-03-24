//
// Created by Liam on 4/03/2022.
//

#ifndef ENGINE_LINE2D_H
#define ENGINE_LINE2D_H
#include "easy_image.h"

class Point2D{
private:
    double x;
    double y;
public:

    Point2D() = default;

    Point2D(double x, double y){
        this->x = x;
        this->y = y;
    }
    //Getters
    double getX() const {
        return x;
    }
    double getY() const {
        return y;
    }
};

class Line2D {
private:
    Point2D p1;
    Point2D p2;
    img::Color color;
    double z1;
    double z2;
public:
    Line2D(Point2D p1, Point2D p2, img::Color color){
        this->p1 = p1;
        this->p2 = p2;
        this->color = color;
    }
    //Getters
    const Point2D &getP1() const {
        return p1;
    }
    const Point2D &getP2() const {
        return p2;
    }
    const img::Color &getColor() const {
        return color;
    }

    void setP1(const Point2D &p1) {
        Line2D::p1 = p1;
    }

    void setP2(const Point2D &p2) {
        Line2D::p2 = p2;
    }

    double getZ1() const {
        return z1;
    }

    void setZ1(double z1) {
        Line2D::z1 = z1;
    }

    double getZ2() const {
        return z2;
    }

    void setZ2(double z2) {
        Line2D::z2 = z2;
    }
};

#endif //ENGINE_LINE2D_H
