//
// Created by liam on 5/11/22.
//

#ifndef ENGINE_CREATE_FIGURES_H
#define ENGINE_CREATE_FIGURES_H

#endif //ENGINE_CREATE_FIGURES_H

#include "Figure.h"
#include "easy_image.h"
#include "cmath"
#include "Line2D.h"
#include "l_parser.h"
#include "fstream"
#include "Lsystem3DElement.h"
#include "LSystemElement.h"
#include "Draw_Lines.h"
#include "Transformations.h"

using Figures3D = vector<Figure>;
using Lights3D = vector<Light>;
using Lines2D = vector<Line2D>;

img::EasyImage
generateLsystem(const string &fileName, const int size, const img::Color &bgColor, const img::Color &lineColor) {
    ifstream input(fileName);
    LParser::LSystem2D LSystem(input);
    input.close();
    Lines2D lines;
    double angle = LSystem.get_angle();
    double currAngle = LSystem.get_starting_angle();
    vector<LSystemElement> stack;
    string drawString = LSystem.get_initiator();
    Point2D p1(0, 0);
    Point2D p2(0, 0);
    for (int iteration = 0; iteration != LSystem.get_nr_iterations(); iteration++) {
        string newDrawString;
        for (auto oldchar: drawString) {
            if (oldchar == '-' or oldchar == '+' or oldchar == '(' or oldchar == ')') {
                newDrawString += oldchar;
            } else {
                newDrawString += LSystem.get_replacement(oldchar);
            }
        }
        drawString = newDrawString;
    }
    for (auto character: drawString) {
        if (character == '+') {
            currAngle += angle;
        } else if (character == '-') {
            currAngle -= angle;
        } else if (character == '(') {
            LSystemElement data(p2, currAngle);
            stack.push_back(data);
        } else if (character == ')') {
            p1 = stack[stack.size() - 1].getPosition();
            p2 = stack[stack.size() - 1].getPosition();
            currAngle = stack[stack.size() - 1].getAngle();
            stack.pop_back();
        } else {
            Point2D temp(p2.getX() + cos(currAngle * pi / 180), p2.getY() + sin(currAngle * pi / 180));
            p2 = temp;
            if (LSystem.draw(character)) {
                Line2D line(p1, p2, lineColor);
                lines.push_back(line);
            }
            p1 = p2;
        }
    }
    Lights3D lights;
    Matrix matrix;
    return draw2DLines({}, lines, size, "", lights, false, matrix, false, bgColor);
}

Figures3D generateFractal(Figure &fig, const int n, const double scale) {
    Matrix s = scaleFigure(1 / scale);
    Figures3D output = {fig};
    for (int i = 0; i < n; i++) {
        Figures3D temp;
        for (auto f: output) {
            for (int j = 0; j < f.points.size(); j++) {
                Figure figure = f;
                applyTransformation(figure, s);
                Matrix t = translate(f.points[j] - figure.points[j]);
                applyTransformation(figure, t);
                temp.push_back(figure);
            }
        }
        output = temp;
    }
    return output;
}

Figure createCube(const img::Color &color) {
    vector<Vector3D> points;
    vector<Face> faces;
    points.push_back(Vector3D::point(1, -1, -1));
    points.push_back(Vector3D::point(-1, 1, -1));
    points.push_back(Vector3D::point(1, 1, 1));
    points.push_back(Vector3D::point(-1, -1, 1));
    points.push_back(Vector3D::point(1, 1, -1));
    points.push_back(Vector3D::point(-1, -1, -1));
    points.push_back(Vector3D::point(1, -1, 1));
    points.push_back(Vector3D::point(-1, 1, 1));
    faces.push_back(Face({0, 4, 2, 6}));
    faces.push_back(Face({4, 1, 7, 2}));
    faces.push_back(Face({1, 5, 3, 7}));
    faces.push_back(Face({5, 0, 6, 3}));
    faces.push_back(Face({6, 2, 7, 3}));
    faces.push_back(Face({0, 5, 1, 4}));
    Figure cube(points, faces, color);
    return cube;
}

Figure createTetrahedron(const img::Color &color) {
    vector<Vector3D> points;
    vector<Face> faces;
    Vector3D p1 = Vector3D::point(1, -1, -1);
    Vector3D p2 = Vector3D::point(-1, 1, -1);
    Vector3D p3 = Vector3D::point(1, 1, 1);
    Vector3D p4 = Vector3D::point(-1, -1, 1);
    points.push_back(p1);
    points.push_back(p2);
    points.push_back(p3);
    points.push_back(p4);
    vector<int> f1 = {0, 1, 2};
    vector<int> f2 = {1, 3, 2};
    vector<int> f3 = {0, 3, 1};
    vector<int> f4 = {0, 2, 3};
    faces.push_back(f1);
    faces.push_back(f2);
    faces.push_back(f3);
    faces.push_back(f4);
    Figure tetrahedron(points, faces, color);
    return tetrahedron;
}

Figure createOctahedron(const img::Color &color) {
    vector<Vector3D> points;
    vector<Face> faces;
    Vector3D p1 = Vector3D::point(1, 0, 0);
    Vector3D p2 = Vector3D::point(0, 1, 0);
    Vector3D p3 = Vector3D::point(-1, 0, 0);
    Vector3D p4 = Vector3D::point(0, -1, 0);
    Vector3D p5 = Vector3D::point(0, 0, -1);
    Vector3D p6 = Vector3D::point(0, 0, 1);
    points.push_back(p1);
    points.push_back(p2);
    points.push_back(p3);
    points.push_back(p4);
    points.push_back(p5);
    points.push_back(p6);
    vector<int> f1 = {0, 1, 5};
    vector<int> f2 = {1, 2, 5};
    vector<int> f3 = {2, 3, 5};
    vector<int> f4 = {3, 0, 5};
    vector<int> f5 = {1, 0, 4};
    vector<int> f6 = {2, 1, 4};
    vector<int> f7 = {3, 2, 4};
    vector<int> f8 = {0, 3, 4};
    faces.push_back(f1);
    faces.push_back(f2);
    faces.push_back(f3);
    faces.push_back(f4);
    faces.push_back(f5);
    faces.push_back(f6);
    faces.push_back(f7);
    faces.push_back(f8);
    Figure octahedron(points, faces, color);
    return octahedron;
}

Figure createIcosahedron(const img::Color &color) {
    vector<Vector3D> points;
    vector<Face> faces;
    Vector3D p1 = Vector3D::point(0, 0, sqrt(5) / 2);
    points.push_back(p1);
    for (int i = 2; i != 7; i++) {
        Vector3D p = Vector3D::point(cos((i - 2) * (2 * pi / 5)), sin((i - 2) * (2 * pi / 5)), 0.5);
        points.push_back(p);
    }
    for (int i = 7; i != 12; i++) {
        Vector3D p = Vector3D::point(cos(pi / 5 + (i - 7) * (2 * pi / 5)), sin(pi / 5 + (i - 7) * (2 * pi / 5)), -0.5);
        points.push_back(p);
    }
    Vector3D p12 = Vector3D::point(0, 0, -sqrt(5) / 2);
    points.push_back(p12);
    vector<int> f1 = {0, 1, 2};
    vector<int> f2 = {0, 2, 3};
    vector<int> f3 = {0, 3, 4};
    vector<int> f4 = {0, 4, 5};
    vector<int> f5 = {0, 5, 1};
    vector<int> f6 = {1, 6, 2};
    vector<int> f7 = {2, 6, 7};
    vector<int> f8 = {2, 7, 3};
    vector<int> f9 = {3, 7, 8};
    vector<int> f10 = {3, 8, 4};
    vector<int> f11 = {4, 8, 9};
    vector<int> f12 = {4, 9, 5};
    vector<int> f13 = {5, 9, 10};
    vector<int> f14 = {5, 10, 1};
    vector<int> f15 = {1, 10, 6};
    vector<int> f16 = {11, 7, 6};
    vector<int> f17 = {11, 8, 7};
    vector<int> f18 = {11, 9, 8};
    vector<int> f19 = {11, 10, 9};
    vector<int> f20 = {11, 6, 10};
    faces.push_back(f1);
    faces.push_back(f2);
    faces.push_back(f3);
    faces.push_back(f4);
    faces.push_back(f5);
    faces.push_back(f6);
    faces.push_back(f7);
    faces.push_back(f8);
    faces.push_back(f9);
    faces.push_back(f10);
    faces.push_back(f11);
    faces.push_back(f12);
    faces.push_back(f13);
    faces.push_back(f14);
    faces.push_back(f15);
    faces.push_back(f16);
    faces.push_back(f17);
    faces.push_back(f18);
    faces.push_back(f19);
    faces.push_back(f20);
    Figure icosahedron(points, faces, color);
    return icosahedron;
}

Figure createDodecahedron(const img::Color &color) {
    Figure icosahedron = createIcosahedron(color);
    vector<Vector3D> points;
    for (int i = 0; i != 20; i++) {
        vector<int> index = icosahedron.faces[i].point_indexes;
        double x =
                (icosahedron.points[index[0]].x + icosahedron.points[index[1]].x + icosahedron.points[index[2]].x) / 3;
        double y =
                (icosahedron.points[index[0]].y + icosahedron.points[index[1]].y + icosahedron.points[index[2]].y) / 3;
        double z =
                (icosahedron.points[index[0]].z + icosahedron.points[index[1]].z + icosahedron.points[index[2]].z) / 3;
        points.push_back(Vector3D::point(x, y, z));
    }
    vector<Face> faces;
    vector<int> f1 = {0, 1, 2, 3, 4};
    vector<int> f2 = {0, 5, 6, 7, 1};
    vector<int> f3 = {1, 7, 8, 9, 2};
    vector<int> f4 = {2, 9, 10, 11, 3};
    vector<int> f5 = {3, 11, 12, 13, 4};
    vector<int> f6 = {4, 13, 14, 5, 0};
    vector<int> f7 = {19, 18, 17, 16, 15};
    vector<int> f8 = {19, 14, 13, 12, 18};
    vector<int> f9 = {18, 12, 11, 10, 17};
    vector<int> f10 = {17, 10, 9, 8, 16};
    vector<int> f11 = {16, 8, 7, 6, 15};
    vector<int> f12 = {15, 6, 5, 14, 19};
    faces.push_back(f1);
    faces.push_back(f2);
    faces.push_back(f3);
    faces.push_back(f4);
    faces.push_back(f5);
    faces.push_back(f6);
    faces.push_back(f7);
    faces.push_back(f8);
    faces.push_back(f9);
    faces.push_back(f10);
    faces.push_back(f11);
    faces.push_back(f12);
    Figure dodecahedron(points, faces, color);
    return dodecahedron;
}

Figure createSphere(const double radius, const int n, const img::Color &color) {
    Figure icosahedron = createIcosahedron(color);
    vector<Face> faces = icosahedron.faces;
    for (int i = 0; i != n; i++) {
        vector<Face> tempFaces;
        for (auto face: faces) {
            Vector3D A = icosahedron.points[face.point_indexes[0]];
            Vector3D B = icosahedron.points[face.point_indexes[1]];
            Vector3D C = icosahedron.points[face.point_indexes[2]];
            icosahedron.points.push_back(Vector3D::point((A.x + B.x) / 2, (A.y + B.y) / 2, (A.z + B.z) / 2));
            icosahedron.points.push_back(Vector3D::point((A.x + C.x) / 2, (A.y + C.y) / 2, (A.z + C.z) / 2));
            icosahedron.points.push_back(Vector3D::point((C.x + B.x) / 2, (C.y + B.y) / 2, (C.z + B.z) / 2));
            tempFaces.push_back(Face({face.point_indexes[0], (int) icosahedron.points.size() - 3,
                                      (int) icosahedron.points.size() - 2}));
            tempFaces.push_back(Face({face.point_indexes[1], (int) icosahedron.points.size() - 1,
                                      (int) icosahedron.points.size() - 3}));
            tempFaces.push_back(Face({face.point_indexes[2], (int) icosahedron.points.size() - 2,
                                      (int) icosahedron.points.size() - 1}));
            tempFaces.push_back(Face({(int) icosahedron.points.size() - 3, (int) icosahedron.points.size() - 1,
                                      (int) icosahedron.points.size() - 2}));
        }
        faces = tempFaces;
    }
    for (auto &point: icosahedron.points) {
        point.normalise();
    }
    return {icosahedron.points, faces, color};
}

Figure createCone(const int n, const double h, const img::Color &color) {
    vector<Vector3D> points;
    vector<Face> faces;
    for (int i = 0; i < n; i++) {
        points.push_back(Vector3D::point(cos((2 * i * pi) / n), sin((2 * i * pi) / n), 0));
        faces.push_back(Face({i, (i + 1) % n, n}));
    }
    points.push_back(Vector3D::point(0, 0, h));
    vector<int> fn;
    for (int i = n - 1; i >= 0; i--) {
        fn.push_back(i);
    }
    faces.push_back(Face(fn));
    return {points, faces, color};
}

Figure createCylinder(const int n, const double h, const img::Color &color, const bool thicccccc = false) {
    vector<Vector3D> points;
    vector<Face> faces;
    for (int i = 0; i < n; i++) {
        points.push_back(Vector3D::point(cos((2 * i * pi) / n), sin((2 * i * pi) / n), 0));
    }
    for (int i = 0; i < n; i++) {
        points.push_back(Vector3D::point(cos((2 * i * pi) / n), sin((2 * i * pi) / n), h));
    }
    for (int i = 0; i < n; i++) {
        faces.push_back(Face({i, (i + 1) % n, (i + 1) % n + n, i + n}));
    }
    if (!thicccccc) {
        vector<int> f1;
        vector<int> f2;
        for (int i = n - 1; i >= 0; i--) {
            f1.push_back(i);
        }
        faces.push_back(Face(f1));
        for (int i = 0; i < n; i++) {
            f2.push_back(i + n);
        }
        faces.push_back(Face(f2));
    }
    return {points, faces, color};
}

Figure createDonut(const double r, const double R, const int n, const int m, const img::Color &color) {
    vector<Vector3D> points;
    vector<Face> faces;
    for (int i = 0; i != n; i++) {
        for (int j = 0; j != m; j++) {
            double u = (2 * i * pi) / n;
            double v = (2 * j * pi) / m;
            points.push_back(Vector3D::point((R + r * cos(v)) * cos(u), (R + r * cos(v)) * sin(u), r * sin(v)));
        }
    }
    for (int i = 0; i != n; i++) {
        for (int j = 0; j != m; j++) {
            faces.push_back(Face({i * m + j, (((i + 1)) % n) * m + j, (((i + 1)) % n) * m + ((j + 1) % m),
                                  i * m + ((j + 1) % m)}));
        }
    }
    return {points, faces, color};
}

Figure generate3DLsystem(const string &fileName, const img::Color &color) {
    ifstream input(fileName);
    LParser::LSystem3D LSystem(input);
    input.close();
    Lines2D lines;
    double angle = LSystem.get_angle();
    Vector3D p1 = Vector3D::point(0, 0, 0);
    Vector3D p2 = Vector3D::point(0, 0, 0);
    Vector3D H = Vector3D::vector(1, 0, 0);
    Vector3D L = Vector3D::vector(0, 1, 0);
    Vector3D U = Vector3D::vector(0, 0, 1);
    vector<Vector3D> points;
    vector<Face> faces;
    vector<Lsystem3DElement> stack;
    string drawString = LSystem.get_initiator();
    for (int iteration = 0; iteration != LSystem.get_nr_iterations(); iteration++) {
        string newDrawString;
        for (auto oldchar: drawString) {
            if (LSystem.get_alphabet().find(oldchar) != LSystem.get_alphabet().end()) {
                newDrawString += LSystem.get_replacement(oldchar);
            } else {
                newDrawString += oldchar;
            }
        }
        drawString = newDrawString;
    }
    for (auto character: drawString) {
        Vector3D tempH = H;
        Vector3D tempL = L;
        Vector3D tempU = U;
        if (character == '+') {
            H = tempH * cos((angle * pi) / 180) + tempL * sin((angle * pi) / 180);
            L = -tempH * sin((angle * pi) / 180) + tempL * cos((angle * pi) / 180);
        } else if (character == '-') {
            H = tempH * cos((-angle * pi) / 180) + tempL * sin((-angle * pi) / 180);
            L = -tempH * sin((-angle * pi) / 180) + tempL * cos((-angle * pi) / 180);
        } else if (character == '^') {
            H = tempH * cos((angle * pi) / 180) + tempU * sin((angle * pi) / 180);
            U = -tempH * sin((angle * pi) / 180) + tempU * cos((angle * pi) / 180);
        } else if (character == '&') {
            H = tempH * cos((-angle * pi) / 180) + tempU * sin((-angle * pi) / 180);
            U = -tempH * sin((-angle * pi) / 180) + tempU * cos((-angle * pi) / 180);
        } else if (character == '\\') {
            L = tempL * cos((angle * pi) / 180) - tempU * sin((angle * pi) / 180);
            U = tempL * sin((angle * pi) / 180) + tempU * cos((angle * pi) / 180);
        } else if (character == '/') {
            L = tempL * cos((-angle * pi) / 180) - tempU * sin((-angle * pi) / 180);
            U = tempL * sin((-angle * pi) / 180) + tempU * cos((-angle * pi) / 180);
        } else if (character == '|') {
            H = -H;
            L = -L;
        } else if (character == '(') {
            Lsystem3DElement element(p2, H, L, U);
            stack.push_back(element);
        } else if (character == ')') {
            p1 = stack[stack.size() - 1].position;
            p2 = stack[stack.size() - 1].position;
            H = stack[stack.size() - 1].H;
            L = stack[stack.size() - 1].L;
            U = stack[stack.size() - 1].U;
            stack.pop_back();
        } else {
            p2 = Vector3D::point(p2.x + H.x, p2.y + H.y, p2.z + H.z);
            if (LSystem.draw(character)) {
                points.push_back(p1);
                points.push_back(p2);
                faces.push_back(Face({(int) points.size() - 2, (int) points.size() - 1}));
            }
            p1 = p2;
        }
    }
    return {points, faces, color};
}

void generateThiccccccFigure(const Figure &lineDrawing, Figures3D &resultingFigures, const double r, const int n,
                             const int m, const img::Color figColor = {}) {
    Matrix scale = scaleFigure(r);
    for (auto &p: lineDrawing.points) {
        Figure fig = createSphere(r, m, figColor);
        applyTransformation(fig, scale);
        Matrix t = translate(p);
        applyTransformation(fig, t);
        resultingFigures.push_back(fig);
    }
    for (auto &f: lineDrawing.faces) {
        for (int i = 0; i < f.point_indexes.size(); i++) {
            Vector3D Pr = Vector3D::vector(lineDrawing.points[f.point_indexes[(i + 1) % f.point_indexes.size()]] -
                                           lineDrawing.points[f.point_indexes[i]]);
            double height = Pr.length() / r;
            Figure fig = createCylinder(n, height, figColor, true);
            applyTransformation(fig, scale);
            double theta;
            double phi;
            double l;
            toPolar(Pr, theta, phi, l);
            Matrix Yrotation = rotateY(phi * 180 / pi);
            applyTransformation(fig, Yrotation);
            Matrix Zrotation = rotateZ(theta * 180 / pi);
            applyTransformation(fig, Zrotation);
            Matrix t = translate(lineDrawing.points[f.point_indexes[i]]);
            applyTransformation(fig, t);
            resultingFigures.push_back(fig);
        }
    }
}