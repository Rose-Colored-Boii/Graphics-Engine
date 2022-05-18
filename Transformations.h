//
// Created by liam on 5/11/22.
//

#ifndef ENGINE_TRANSFORMATIONS_H
#define ENGINE_TRANSFORMATIONS_H

#endif //ENGINE_TRANSFORMATIONS_H

#include "vector3d.h"
#include "Figure.h"
#include "cmath"
#include "Line2D.h"

using Figures3D = vector<Figure>;
using Lines2D = vector<Line2D>;

double pVal(const Vector3D &A, const Vector3D &B, const double dNear, const double dFar, const double dVal,
            const double dRight, const string type) {
    double p = 0;
    if (type == "NearFar") {
        p = (dVal - B.z) / (A.z - B.z);
    } else if (type == "LeftRight") {
        p = (B.x * dNear + B.z * dVal) / ((B.x - A.x) * dNear + (B.z - A.z) * dVal);
    } else {
        p = (B.y * dNear + B.z * dVal) / ((B.y - A.y) * dNear + (B.z - A.z) * dVal);
    }
    return p;
}

vector<Face> triangulate(const Face &face) {
    vector<Face> faces;
    for (int i = 1; i < face.point_indexes.size() - 1; i++) {
        vector<int> indexes;
        indexes.push_back(face.point_indexes[0]);
        indexes.push_back(face.point_indexes[i]);
        indexes.push_back(face.point_indexes[i + 1]);
        faces.push_back(Face(indexes));
    }
    return faces;
}

vector<Face>
clipPlane(vector<Face> &faces, vector<Vector3D> &points, const double val1, const double val2, const double dNear,
          const double dFar, const double hfov, const double aspectRatio, const string &type) {
    double right = dNear * tan((hfov / 2) * (pi / 180));
    double top = right / aspectRatio;
    vector<Face> f1;
    for (const auto &face: faces) {
        int count = 0;
        for (auto point: face.point_indexes) {
            double coord = 0;
            if (type == "LeftRight") {
                coord = (points[point].x * dNear) / -points[point].z;
            } else if (type == "BottomUp") {
                coord = (points[point].y * dNear) / -points[point].z;
            } else {
                coord = points[point].z;
            }
            if (coord < val1) {
                count++;
            }
        }
        if (count == 0) {
            f1.push_back(face);
        }
        if (count == 1) {
            Vector3D in1, in2, out1;
            double coord;
            for (int i = 0; i <= 2; i++) {
                if (type == "LeftRight") {
                    coord = (points[face.point_indexes[i]].x * dNear) / -points[face.point_indexes[i]].z;
                } else if (type == "BottomUp") {
                    coord = (points[face.point_indexes[i]].y * dNear) / -points[face.point_indexes[i]].z;
                } else {
                    coord = points[face.point_indexes[i]].z;
                }
                if (coord < val1) {
                    out1 = points[face.point_indexes[i]];
                    in1 = points[face.point_indexes[(i + 2) % 3]];
                    in2 = points[face.point_indexes[(i + 1) % 3]];
                    break;
                }
            }
            double p1;
            double p2;
            Vector3D point1;
            Vector3D point2;
            p1 = pVal(out1, in1, dNear, dFar, val1, right, type);
            point1 = p1 * out1 + (1 - p1) * in1;
            p2 = pVal(out1, in2, dNear, dFar, val1, right, type);
            point2 = p2 * out1 + (1 - p2) * in2;
            points.push_back(in1);
            points.push_back(point1);
            points.push_back(point2);
            points.push_back(in2);
            f1.push_back(Face({(int) points.size() - 4, (int) points.size() - 3, (int) points.size() - 2}));
            f1.push_back(Face({(int) points.size() - 4, (int) points.size() - 2, (int) points.size() - 1}));
        }
        if (count == 2) {
            Vector3D in1, out1, out2;
            double coord;
            for (int i = 0; i <= 2; i++) {
                if (type == "LeftRight") {
                    coord = (points[face.point_indexes[i]].x * dNear) / -points[face.point_indexes[i]].z;
                } else if (type == "BottomUp") {
                    coord = (points[face.point_indexes[i]].y * dNear) / -points[face.point_indexes[i]].z;
                } else {
                    coord = points[face.point_indexes[i]].z;
                }
                if (coord >= val1) {
                    out1 = points[face.point_indexes[(i + 2) % 3]];
                    in1 = points[face.point_indexes[i]];
                    out2 = points[face.point_indexes[(i + 1) % 3]];
                    break;
                }
            }
            double p1;
            double p2;
            Vector3D point1;
            Vector3D point2;
            p1 = pVal(out1, in1, dNear, dFar, val1, right, type);
            point1 = p1 * out1 + (1 - p1) * in1;
            p2 = pVal(out2, in1, dNear, dFar, val1, right, type);
            point2 = p2 * out2 + (1 - p2) * in1;
            points.push_back(in1);
            points.push_back(point1);
            points.push_back(point2);
            f1.push_back(Face({(int) points.size() - 3, (int) points.size() - 2, (int) points.size() - 1}));
        }
    }
    vector<Face> f2;
    for (const auto &face: f1) {
        int count = 0;
        for (auto point: face.point_indexes) {
            double coord;
            if (type == "LeftRight") {
                coord = (points[point].x * dNear) / -points[point].z;
            } else if (type == "BottomUp") {
                coord = (points[point].y * dNear) / -points[point].z;
            } else {
                coord = points[point].z;
            }
            if (coord > val2) {
                count++;
            }
        }
        if (count == 0) {
            f2.push_back(face);
        }
        if (count == 1) {
            Vector3D in1, in2, out1;
            double coord;
            for (int i = 0; i <= 2; i++) {
                if (type == "LeftRight") {
                    coord = (points[face.point_indexes[i]].x * dNear) / -points[face.point_indexes[i]].z;
                } else if (type == "BottomUp") {
                    coord = (points[face.point_indexes[i]].y * dNear) / -points[face.point_indexes[i]].z;
                } else {
                    coord = points[face.point_indexes[i]].z;
                }
                if (coord > val2) {
                    out1 = points[face.point_indexes[i]];
                    in1 = points[face.point_indexes[(i + 2) % 3]];
                    in2 = points[face.point_indexes[(i + 1) % 3]];
                    break;
                }
            }
            double p1;
            double p2;
            Vector3D point1;
            Vector3D point2;
            p1 = pVal(out1, in1, dNear, dFar, val2, right, type);
            point1 = p1 * out1 + (1 - p1) * in1;
            p2 = pVal(out1, in2, dNear, dFar, val2, right, type);
            point2 = p2 * out1 + (1 - p2) * in2;
            points.push_back(in1);
            points.push_back(point1);
            points.push_back(point2);
            points.push_back(in2);
            f2.push_back(Face({(int) points.size() - 4, (int) points.size() - 3, (int) points.size() - 2}));
            f2.push_back(Face({(int) points.size() - 4, (int) points.size() - 1, (int) points.size() - 2}));
        }
        if (count == 2) {
            Vector3D in1, out1, out2;
            double coord;
            for (int i = 0; i <= 2; i++) {
                if (type == "LeftRight") {
                    coord = (points[face.point_indexes[i]].x * dNear) / -points[face.point_indexes[i]].z;
                } else if (type == "BottomUp") {
                    coord = (points[face.point_indexes[i]].y * dNear) / -points[face.point_indexes[i]].z;
                } else {
                    coord = points[face.point_indexes[i]].z;
                }
                if (coord <= val2) {
                    out1 = points[face.point_indexes[(i + 2) % 3]];
                    in1 = points[face.point_indexes[i]];
                    out2 = points[face.point_indexes[(i + 1) % 3]];
                    break;
                }
            }
            double p1;
            double p2;
            Vector3D point1;
            Vector3D point2;
            p1 = pVal(out1, in1, dNear, dFar, val2, right, type);
            point1 = p1 * out1 + (1 - p1) * in1;
            p2 = pVal(out2, in1, dNear, dFar, val2, right, type);
            point2 = p2 * out2 + (1 - p2) * in1;
            points.push_back(in1);
            points.push_back(point1);
            points.push_back(point2);
            f2.push_back(Face({(int) points.size() - 3, (int) points.size() - 2, (int) points.size() - 1}));
        }
    }
    return f2;
}

void clip(Figures3D &figs, const double dNear, const double dFar, const double hfov, const double aspectRatio) {
    for (auto &fig: figs) {
        vector<Vector3D> points = fig.points;
        vector<Face> f1 = clipPlane(fig.faces, points, -dFar, -dNear, dNear, dFar, hfov, aspectRatio, "NearFar");
        double right = dNear * tan((hfov / 2) * (pi / 180));
        vector<Face> f2 = clipPlane(f1, points, -right, right, dNear, dFar, hfov, aspectRatio, "LeftRight");
        double top = right / aspectRatio;
        vector<Face> f3 = clipPlane(f2, points, -top, top, dNear, dFar, hfov, aspectRatio, "BottomUp");
        fig.faces = f3;
        fig.points = points;
    }
}

Matrix scaleFigure(const double scaleFactor) {
    Matrix scalingMatrix;
    scalingMatrix(1, 1) = scaleFactor;
    scalingMatrix(2, 2) = scaleFactor;
    scalingMatrix(3, 3) = scaleFactor;
    scalingMatrix(4, 4) = 1;
    return scalingMatrix;
}

Matrix rotateX(const double angle) {
    Matrix rotationMatrix;
    rotationMatrix(2, 2) = cos(angle * pi / 180);
    rotationMatrix(3, 3) = cos(angle * pi / 180);
    rotationMatrix(3, 2) = -sin(angle * pi / 180);
    rotationMatrix(2, 3) = sin(angle * pi / 180);
    return rotationMatrix;
}

Matrix rotateY(const double angle) {
    Matrix rotationMatrix;
    rotationMatrix(1, 1) = cos(angle * pi / 180);
    rotationMatrix(3, 3) = cos(angle * pi / 180);
    rotationMatrix(1, 3) = -sin(angle * pi / 180);
    rotationMatrix(3, 1) = sin(angle * pi / 180);
    return rotationMatrix;
}

Matrix rotateZ(const double angle) {
    Matrix rotationMatrix;
    rotationMatrix(1, 1) = cos(angle * pi / 180);
    rotationMatrix(2, 2) = cos(angle * pi / 180);
    rotationMatrix(2, 1) = -sin(angle * pi / 180);
    rotationMatrix(1, 2) = sin(angle * pi / 180);
    return rotationMatrix;
}

Matrix translate(const Vector3D &vector) {
    Matrix translationMatrix;
    translationMatrix(4, 1) = vector.x;
    translationMatrix(4, 2) = vector.y;
    translationMatrix(4, 3) = vector.z;
    return translationMatrix;
}

void applyTransformation(Figure &fig, const Matrix &m) {
    for (auto &point: fig.points) {
        point *= m;
    }
}

void applyTransformation(Figures3D &figs, const Matrix &m) {
    for (auto &figure: figs) {
        applyTransformation(figure, m);
    }
}

void toPolar(const Vector3D &point, double &theta, double &phi, double &r) {
    r = sqrt(pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2));
    theta = atan2(point.y, point.x);
    phi = acos(point.z / r);
}

Matrix eyeMatrix(const double theta, const double phi, const double r) {
    Matrix transformationMatrix;
    transformationMatrix(1, 1) = -sin(theta);
    transformationMatrix(2, 1) = cos(theta);
    transformationMatrix(1, 2) = -cos(theta) * cos(phi);
    transformationMatrix(1, 3) = cos(theta) * sin(phi);
    transformationMatrix(2, 2) = -sin(theta) * cos(phi);
    transformationMatrix(2, 3) = sin(theta) * sin(phi);
    transformationMatrix(3, 2) = sin(phi);
    transformationMatrix(3, 3) = cos(phi);
    transformationMatrix(4, 3) = -r;
    return transformationMatrix;
}

Point2D doProjection(const Vector3D &point, const double d = 1) {
    double x = (d * point.x) / (-point.z);
    double y = (d * point.y) / (-point.z);
    Point2D p(x, y);
    return p;
}

Lines2D doProjection(const Figures3D &figures) {
    Lines2D lines;
    for (auto figure: figures) {
        for (const auto &face: figure.faces) {
            vector<Point2D> points;
            vector<double> zCo;
            for (auto index: face.point_indexes) {
                points.push_back(doProjection(figure.points[index]));
                zCo.push_back(figure.points[index].z);
            }
            for (int i = 0; i != points.size(); i++) {
                if (i != points.size() - 1) {
                    Line2D line(points[i], points[i + 1], figure.color);
                    line.setZ1(zCo[i]);
                    line.setZ2(zCo[i + 1]);
                    lines.push_back(line);
                } else {
                    Line2D line(points[i], points[0], figure.color);
                    line.setZ1(zCo[i]);
                    line.setZ2(zCo[0]);
                    lines.push_back(line);
                }
            }
        }
    }
    return lines;
}

void generateShadowMask(Light &light, Figures3D &figures, int size, bool clipping) {
    applyTransformation(figures, light.eye);
    Lines2D lines = doProjection(figures);
    double xMin = min(lines[0].getP1().getX(), lines[0].getP2().getX());
    double xMax = max(lines[0].getP1().getX(), lines[0].getP2().getX());
    double yMin = min(lines[0].getP1().getY(), lines[0].getP2().getY());
    double yMax = max(lines[0].getP1().getY(), lines[0].getP2().getY());
    for (const auto &line: lines) {
        double maxX = max(line.getP1().getX(), line.getP2().getX());
        double minX = min(line.getP1().getX(), line.getP2().getX());
        double maxY = max(line.getP1().getY(), line.getP2().getY());
        double minY = min(line.getP1().getY(), line.getP2().getY());
        if (minX < xMin)
            xMin = minX;
        if (maxX > xMax)
            xMax = maxX;
        if (minY < yMin)
            yMin = minY;
        if (maxY > yMax)
            yMax = maxY;
    }
    double xRange = xMax - xMin;
    double yRange = yMax - yMin;
    double imageX = size * (xRange) / max(xRange, yRange);
    double imageY = size * (yRange) / max(xRange, yRange);
    img::EasyImage image(lround(imageX), lround(imageY), {});
    ZBuffer zBuffer((int) lround(imageX), (int) lround(imageY));
    double d = 0.95 * (imageX) / xRange;
    light.d = d;
    double DCx = d * (xMin + xMax) / 2;
    double DCy = d * (yMin + yMax) / 2;
    double dx = (imageX / 2) - DCx;
    light.dx = dx;
    double dy = (imageY / 2) - DCy;
    light.dy = dy;
    Matrix matrix;
    ZBuffer shadowMask = ZBuffer((int) imageX, (int) imageY);
    for (auto &fig: figures) {
        for (auto &face: fig.faces) {
            Lights3D lights;
            draw_zbuf_triag(shadowMask, image, fig.points[face.point_indexes[0]],
                            fig.points[face.point_indexes[1]],
                            fig.points[face.point_indexes[2]], light.d, light.dx, light.dy, fig.ambientReflection,
                            fig.diffuseReflection, fig.specularReflection, fig.reflectionCoefficient, lights,
                            "LightedZBuffering",
                            clipping, light.eye, false, fig.texture);
        }
    }
    light.shadowMask = shadowMask;
}