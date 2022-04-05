#include "easy_image.h"
#include "ini_configuration.h"
#include "Line2D.h"
#include "l_parser.h"
#include "LSystemElement.h"
#include "vector3d.h"
#include "Figure.h"
#include <fstream>
#include <iostream>
#include "limits"
#include <string>
#include <cmath>
#include "algorithm"
#include "ZBuffer.h"
#include "Lsystem3DElement.h"
using namespace std;
using Lines2D = vector<Line2D>;
using Figures3D = vector<Figure>;
constexpr double pi = 3.14159265358979323846264338327950288;

void draw_zbuf_line(ZBuffer& zBuffer, img::EasyImage& image, unsigned int x0, unsigned int y0, double z0, unsigned int x1, unsigned int y1, double z1, const img::Color& color){
    int a = 0;
    if (x0 == x1 and y0 == y1){
        double z = 1/z0;
        if (z < zBuffer.matrix[x0][y0]){
            zBuffer.matrix[x0][y0] = z;
            image(x0, y0) = color;
        }
    }
    else if (x0 == x1)
    {
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++){
            a++;
        }
        //special case for x0 == x1
        int j = 0;
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++)
        {
            double z = 0;
            if (y0 < y1){
                z = (j/(double) (a-1))/z1 + (1-(j/(double) (a-1)))/z0;
            }
            else{
                z = (j/(double) (a-1))/z0 + (1-(j/(double) (a-1)))/z1;
            }
            if (z < zBuffer.matrix[x0][i]){
                zBuffer.matrix[x0][i] = z;
                image(x0, i) = color;
            }
            if (j == 0){
                if (y0 < y1){
                    if (z != 1/z0){
                        cerr << "pain" << endl;
                        exit(1);
                    }
                }
                else{
                    if (z != 1/z1){
                        cerr << "pain" << endl;
                        exit(1);
                    }
                }
            }
            if (j == a-1){
                if (y0 < y1){
                    if (z != 1/z1){
                        cerr << "pain" << endl;
                        exit(1);
                    }
                } else{
                    if (z != 1/z0){
                        cerr << "pain" << endl;
                        exit(1);
                    }
                }
            }
            if (j > a-1){
                cerr << "pain" << endl;
                exit(1);
            }
            j++;
        }
    }
    else if (y0 == y1)
    {
        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++)
        {
            a++;
        }
        //special case for y0 == y1
        int j = 0;
        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++)
        {
            double z = 0;
            if (x0 < x1){
                z = (j/(double) (a-1))/z1 + (1-(j/(double) (a-1)))/z0;
            }
            else{
                z = (j/(double) (a-1))/z0 + (1-(j/(double) (a-1)))/z1;
            }
            if (z < zBuffer.matrix[i][y0]) {
                zBuffer.matrix[i][y0] = z;
                image(i, y0) = color;
            }
            if (j == 0){
                if (x0 < x1){
                    if (z != 1/z0){
                        cerr << "pain" << endl;
                        exit(1);
                    }
                }
                else{
                    if (z != 1/z1){
                        cerr << "pain" << endl;
                        exit(1);
                    }
                }
            }
            if (j == a-1){
                if (x0 < x1){
                    if (z != 1/z1){
                        cerr << "pain" << endl;
                        exit(1);
                    }
                } else{
                    if (z != 1/z0){
                        cerr << "pain" << endl;
                        exit(1);
                    }
                }
            }
            if (j > a-1){
                cerr << "pain" << endl;
                exit(1);
            }
            j++;
        }
    }
    else
    {
        if (x0 > x1)
        {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(z0, z1);
            std::swap(y0, y1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0)
        {
            for (unsigned int i = 0; i <= (x1 - x0); i++)
            {
                a++;
            }
            int j = 0;
            for (unsigned int i = 0; i <= (x1 - x0); i++)
            {
                double z = (j/(double) (a-1))/z1 + (1-(j/(double) (a-1)))/z0;
                if (z < zBuffer.matrix[x0+i][lround(y0+m*i)]) {
                    zBuffer.matrix[x0+i][lround(y0+m*i)] = z;
                    image(x0 + i, (unsigned int) round(y0 + m * i)) = color;
                }
                if (j == 0){
                    if (z != 1/z0){
                        cerr << "pain" << endl;
                        exit(1);
                    }
                }
                if (j == a-1){
                    if (z != 1/z1){
                        cerr << "pain" << endl;
                        exit(1);
                    }
                }
                if (j > a-1){
                    cerr << "pain" << endl;
                    exit(1);
                }
                j++;
            }
        }
        else if (m > 1.0)
        {
            for (unsigned int i = 0; i <= (y1 - y0); i++)
            {
                a++;
            }
            int j = 0;
            for (unsigned int i = 0; i <= (y1 - y0); i++)
            {
                double z = (j/(double) (a-1))/z1 + (1-(j/(double) (a-1)))/z0;
                if (z < zBuffer.matrix[lround(x0 + (i / m))][y0+i]) {
                    zBuffer.matrix[lround(x0 + (i / m))][y0+i] = z;
                    image((unsigned int) round(x0 + (i / m)), y0 + i) = color;
                }
                if (j == 0){
                    if (z != 1/z0){
                        cerr << "pain" << endl;
                        exit(1);
                    }
                }
                if (j == a-1){
                    if (z != 1/z1){
                        cerr << "pain" << endl;
                        exit(1);
                    }
                }
                if (j > a-1){
                    cerr << "pain" << endl;
                    exit(1);
                }
                j++;
            }
        }
        else if (m < -1.0)
        {
            for (unsigned int i = 0; i <= (y0 - y1); i++)
            {
                a++;
            }
            int j = 0;
            for (unsigned int i = 0; i <= (y0 - y1); i++)
            {
                double z = (j/(double) (a-1))/z1 + (1-(j/(double) (a-1)))/z0;
                if (z < zBuffer.matrix[lround(x0 - (i / m))][y0-i]) {
                    zBuffer.matrix[lround(x0 - (i / m))][y0-i] = z;
                    image((unsigned int) round(x0 - (i / m)), y0 - i) = color;
                }
                if (j == 0){
                    if (z != 1/z0){
                        cerr << "pain" << endl;
                        exit(1);
                    }
                }
                if (j == a-1){
                    if (z != 1/z1){
                        cerr << "pain" << endl;
                        exit(1);
                    }
                }
                if (j > a-1){
                    cerr << "pain" << endl;
                    exit(1);
                }
                j++;
            }
        }
    }
}

void draw_zbuf_line(ZBuffer& zBuffer, img::EasyImage& image, const unsigned int x0, const unsigned int x1, const unsigned int y, const double xG, const double yG, const double zG,
                    const double dzdx, const double dzdy, const img::Color& color){
    for (unsigned int i = x0; i <= x1; i++){
        double z = 1.0001*zG + (i-xG)*dzdx + (y-yG)*dzdy;
        if (z < zBuffer.matrix[i][y]){
            zBuffer.matrix[i][y] = z;
            image(i, y) = color;
        }
    }
}

void draw_zbuf_triag(ZBuffer& zBuffer, img::EasyImage& image, const Vector3D& A, const Vector3D& B, const Vector3D& C, double d, double dx, double dy, img::Color& color){
    Point2D newA((d*A.x)/(-A.z) + dx, (d*A.y)/(-A.z) + dy);
    Point2D newB((d*B.x)/(-B.z) + dx, (d*B.y)/(-B.z) + dy);
    Point2D newC((d*C.x)/(-C.z) + dx, (d*C.y)/(-C.z) + dy);
    Point2D G((newA.getX()+newB.getX()+newC.getX())/3, (newA.getY()+newB.getY()+newC.getY())/3);
    double zG = 1/(3*A.z)+1/(3*B.z)+1/(3*C.z);
    Vector3D u = B-A;
    Vector3D v = C-A;
    double w1 = u.y*v.z - u.z*v.y;
    double w2 = u.z*v.x - u.x*v.z;
    double w3 = u.x*v.y - u.y*v.x;
    double k = w1*A.x + w2*A.y + w3*A.z;
    double dzdx = w1/(-d*k);
    double dzdy = w2/(-d*k);
    int yMin = lround(min(newA.getY(), min(newB.getY(), newC.getY()))+0.5);
    int yMax = lround(max(newA.getY(), max(newB.getY(), newC.getY()))-0.5);
    vector<Point2D> points = {newA, newB, newC};
    for (int i = yMin; i <= yMax; i++){
        double x1 = numeric_limits<double>::infinity();
        double x2 = numeric_limits<double>::infinity();
        double x3 = numeric_limits<double>::infinity();
        double x4 = -numeric_limits<double>::infinity();
        double x5 = -numeric_limits<double>::infinity();
        double x6 = -numeric_limits<double>::infinity();
        if ((i-points[0].getY())*(i-points[1].getY()) <= 0){
            x1 = x4 = points[1].getX() + (points[0].getX()-points[1].getX())*((i-points[1].getY())/(points[0].getY()-points[1].getY()));
        }
        if ((i-points[1].getY())*(i-points[2].getY()) <= 0){
            x2 = x5 = points[2].getX() + (points[1].getX()-points[2].getX())*((i-points[2].getY())/(points[1].getY()-points[2].getY()));
        }
        if ((i-points[0].getY())*(i-points[2].getY()) <= 0){
            x3 = x6 = points[2].getX() + (points[0].getX()-points[2].getX())*((i-points[2].getY())/(points[0].getY()-points[2].getY()));
        }
        int xL = lround(min(x1, min(x2, x3))+0.5);
        int xR = lround(max(x4, max(x5, x6))-0.5);
        draw_zbuf_line(zBuffer, image, xL, xR, i, G.getX(), G.getY(), zG, dzdx, dzdy, color);
    }
}

img::EasyImage draw2DLines(const Figures3D& figs, const Lines2D& lines, const int size, const string& type, const img::Color& bgColor = img::Color()){
    double xMin = min(lines[0].getP1().getX(), lines[0].getP2().getX());
    double xMax = max(lines[0].getP1().getX(), lines[0].getP2().getX());
    double yMin = min(lines[0].getP1().getY(), lines[0].getP2().getY());
    double yMax = max(lines[0].getP1().getY(), lines[0].getP2().getY());
    for (const auto& line : lines){
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
    double imageX = size*(xRange)/max(xRange, yRange);
    double imageY = size*(yRange)/max(xRange, yRange);
    img::EasyImage image(lround(imageX), lround(imageY), bgColor);
    ZBuffer zBuffer((int) lround(imageX), (int) lround(imageY));
    double d = 0.95*(imageX)/xRange;
    double DCx = d*(xMin + xMax)/2;
    double DCy = d*(yMin + yMax)/2;
    double dx = (imageX/2) - DCx;
    double dy = (imageY/2) - DCy;
    if (type == "ZBuffering"){
        for (auto fig : figs){
            for (auto face : fig.faces){
                draw_zbuf_triag(zBuffer, image, fig.points[face.point_indexes[0]], fig.points[face.point_indexes[1]], fig.points[face.point_indexes[2]], d, dx, dy, fig.color);
            }
        }
    }
    else {
        for (const auto &line: lines) {
            img::Color color(line.getColor().red, line.getColor().green, line.getColor().blue);
            if (type == "ZBufferedWireframe") {
                draw_zbuf_line(zBuffer, image, lround(line.getP1().getX() * d + dx),
                               lround((line.getP1().getY() * d + dy)), line.getZ1(),
                               lround((line.getP2().getX() * d + dx)), lround((line.getP2().getY() * d + dy)),
                               line.getZ2(), color);
            } else {
                image.draw_line(lround(line.getP1().getX() * d + dx), lround((line.getP1().getY() * d + dy)),
                                lround((line.getP2().getX() * d + dx)), lround((line.getP2().getY() * d + dy)), color);
            }
        }
    }
    return image;
}

img::EasyImage generateLsystem(const string& fileName, const int size, const img::Color& bgColor, const img::Color& lineColor){
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
    for (int iteration = 0; iteration != LSystem.get_nr_iterations(); iteration++){
        string newDrawString;
        for (auto oldchar : drawString){
            if (oldchar == '-' or oldchar == '+' or oldchar == '(' or oldchar == ')'){
                newDrawString += oldchar;
            }
            else{
                newDrawString += LSystem.get_replacement(oldchar);
            }
        }
        drawString = newDrawString;
    }
    for (auto character : drawString){
        if (character == '+'){
            currAngle += angle;
        }
        else if (character == '-'){
            currAngle -= angle;
        }
        else if (character == '('){
            LSystemElement data(p2, currAngle);
            stack.push_back(data);
        }
        else if (character == ')'){
            p1 = stack[stack.size()-1].getPosition();
            p2 = stack[stack.size()-1].getPosition();
            currAngle = stack[stack.size()-1].getAngle();
            stack.pop_back();
        }
        else {
            Point2D temp(p2.getX() + cos(currAngle * pi / 180), p2.getY() + sin(currAngle * pi / 180));
            p2 = temp;
            if (LSystem.draw(character)) {
                Line2D line(p1, p2, lineColor);
                lines.push_back(line);
            }
            p1 = p2;
        }
    }
    return draw2DLines({}, lines, size, "", bgColor);
}

Matrix scaleFigure(const double scaleFactor){
    Matrix scalingMatrix;
    scalingMatrix(1, 1) = scaleFactor;
    scalingMatrix(2, 2) = scaleFactor;
    scalingMatrix(3, 3) = scaleFactor;
    scalingMatrix(4, 4) = 1;
    return scalingMatrix;
}

Matrix rotateX(const double angle){
    Matrix rotationMatrix;
    rotationMatrix(2, 2) = cos(angle * pi/180);
    rotationMatrix(3, 3) = cos(angle * pi/180);
    rotationMatrix(3, 2) = -sin(angle * pi/180);
    rotationMatrix(2, 3) = sin(angle * pi/180);
    return rotationMatrix;
}

Matrix rotateY(const double angle){
    Matrix rotationMatrix;
    rotationMatrix(1, 1) = cos(angle * pi/180);
    rotationMatrix(3, 3) = cos(angle * pi/180);
    rotationMatrix(1, 3) = -sin(angle * pi/180);
    rotationMatrix(3, 1) = sin(angle * pi/180);
    return rotationMatrix;
}

Matrix rotateZ(const double angle){
    Matrix rotationMatrix;
    rotationMatrix(1, 1) = cos(angle * pi/180);
    rotationMatrix(2, 2) = cos(angle * pi/180);
    rotationMatrix(2, 1) = -sin(angle * pi/180);
    rotationMatrix(1, 2) = sin(angle * pi/180);
    return rotationMatrix;
}

Matrix translate(const Vector3D& vector){
    Matrix translationMatrix;
    translationMatrix(4, 1) = vector.x;
    translationMatrix(4, 2) = vector.y;
    translationMatrix(4, 3) = vector.z;
    return translationMatrix;
}

void applyTransformation(Figure& fig, const Matrix &m){
    for (auto& point : fig.points){
        point *= m;
    }
}

void applyTransformation(Figures3D &figs, const Matrix &m){
    for (auto& figure : figs){
        applyTransformation(figure, m);
    }
}

void toPolar(const Vector3D& point, double& theta, double& phi, double& r){
    r = sqrt(pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2));
    theta = atan2(point.y, point.x);
    phi = acos(point.z/r);
}

Matrix eyeMatrix(const double theta, const double phi, const double r){
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

Point2D doProjection(const Vector3D& point, const double d = 1){
    double x = (d*point.x)/(-point.z);
    double y = (d*point.y)/(-point.z);
    Point2D p(x, y);
    return p;
}

Lines2D doProjection(const Figures3D& figures){
    Lines2D lines;
    for (auto figure : figures){
        for (const auto& face : figure.faces){
            vector<Point2D> points;
            vector<double> zCo;
            for (auto index : face.point_indexes){
                points.push_back(doProjection(figure.points[index]));
                zCo.push_back(figure.points[index].z);
            }
            for (int i = 0; i != points.size(); i++){
                if (i != points.size()-1) {
                    Line2D line(points[i], points[i + 1], figure.color);
                    line.setZ1(zCo[i]);
                    line.setZ2(zCo[i+1]);
                    lines.push_back(line);
                }
                else{
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

Figure createCube(const img::Color& color){
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

Figure createTetrahedron(const img::Color& color){
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

Figure createOctahedron(const img::Color& color){
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

Figure createIcosahedron(const img::Color& color){
    vector<Vector3D> points;
    vector<Face> faces;
    Vector3D p1 = Vector3D::point(0, 0, sqrt(5)/2);
    points.push_back(p1);
    for (int i = 2; i != 7; i++){
        Vector3D p = Vector3D::point(cos((i-2)*(2*pi/5)), sin((i-2)*(2*pi/5)), 0.5);
        points.push_back(p);
    }
    for (int i = 7; i != 12; i++){
        Vector3D p = Vector3D::point(cos(pi/5 + (i-7)*(2*pi/5)), sin(pi/5 + (i-7)*(2*pi/5)), -0.5);
        points.push_back(p);
    }
    Vector3D p12 = Vector3D::point(0, 0, -sqrt(5)/2);
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

Figure createDodecahedron(const img::Color& color){
    Figure icosahedron = createIcosahedron(color);
    vector<Vector3D> points;
    for (int i = 0; i != 20; i++){
        vector<int> index = icosahedron.faces[i].point_indexes;
        double x = (icosahedron.points[index[0]].x + icosahedron.points[index[1]].x + icosahedron.points[index[2]].x)/3;
        double y = (icosahedron.points[index[0]].y + icosahedron.points[index[1]].y + icosahedron.points[index[2]].y)/3;
        double z = (icosahedron.points[index[0]].z + icosahedron.points[index[1]].z + icosahedron.points[index[2]].z)/3;
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

Figure createSphere(const double radius, const int n, const img::Color& color){
    Figure icosahedron = createIcosahedron(color);
    vector<Face> faces = icosahedron.faces;
    for (int i = 0; i != n; i++){
        vector<Face> tempFaces;
        for (auto face : faces){
            Vector3D A = icosahedron.points[face.point_indexes[0]];
            Vector3D B = icosahedron.points[face.point_indexes[1]];
            Vector3D C = icosahedron.points[face.point_indexes[2]];
            icosahedron.points.push_back(Vector3D::point((A.x+B.x)/2, (A.y+B.y)/2, (A.z+B.z)/2));
            icosahedron.points.push_back(Vector3D::point((A.x+C.x)/2, (A.y+C.y)/2, (A.z+C.z)/2));
            icosahedron.points.push_back(Vector3D::point((C.x+B.x)/2, (C.y+B.y)/2, (C.z+B.z)/2));
            tempFaces.push_back(Face({face.point_indexes[0], (int) icosahedron.points.size()-3, (int) icosahedron.points.size()-2}));
            tempFaces.push_back(Face({face.point_indexes[1], (int) icosahedron.points.size()-1, (int) icosahedron.points.size()-3}));
            tempFaces.push_back(Face({face.point_indexes[2], (int) icosahedron.points.size()-2, (int) icosahedron.points.size()-1}));
            tempFaces.push_back(Face({(int) icosahedron.points.size()-3, (int) icosahedron.points.size()-1, (int) icosahedron.points.size()-2}));
        }
        faces = tempFaces;
    }
    for (auto& point : icosahedron.points){
        point.normalise();
    }
    return {icosahedron.points, faces, color};
}

Figure createCone(const int n, const double h, const img::Color& color){
    vector<Vector3D> points;
    vector<Face> faces;
    for (int i = 0; i < n; i++){
        points.push_back(Vector3D::point(cos((2*i*pi)/n), sin((2*i*pi)/n), 0));
        faces.push_back(Face({i, (i+1)%n, n}));
    }
    points.push_back(Vector3D::point(0, 0, h));
    vector<int> fn;
    for (int i = n-1; i >= 0; i--){
        fn.push_back(i);
    }
    faces.push_back(Face(fn));
    return {points, faces, color};
}

Figure createCylinder(const int n, const double h, const img::Color& color){
    vector<Vector3D> points;
    vector<Face> faces;
    for (int i = 0; i < n; i++){
        points.push_back(Vector3D::point(cos((2*i*pi)/n), sin((2*i*pi)/n), 0));
    }
    for (int i = 0; i < n; i++){
        points.push_back(Vector3D::point(cos((2*i*pi)/n), sin((2*i*pi)/n), h));
    }
    for (int i = 0; i < n; i++){
        faces.push_back(Face({i, (i+1)%n, (i+1)%n+n, i+n}));
    }
    vector<int> f1;
    vector<int> f2;
    for (int i = n-1; i >= 0; i--){
        f1.push_back(i);
    }
    faces.push_back(Face(f1));
    for (int i = 0; i < n; i++){
        f2.push_back(i+n);
    }
    faces.push_back(Face(f2));
    return {points, faces, color};
}

Figure createDonut(const double r, const double R, const int n, const int m, const img::Color& color){
    vector<Vector3D> points;
    vector<Face> faces;
    for (int i = 0; i != n; i++){
        for (int j = 0; j != m; j++){
            double u = (2*i*pi)/n;
            double v = (2*j*pi)/m;
            points.push_back(Vector3D::point((R+r*cos(v))*cos(u), (R+r*cos(v))*sin(u), r*sin(v)));
        }
    }
    for (int i = 0; i != n; i++){
        for (int j = 0; j != m; j++){
            faces.push_back(Face({i*m+j, (((i+1))%n)*m+j, (((i+1))%n)*m+((j+1)%m), i*m+((j+1)%m)}));
        }
    }
    return {points, faces, color};
}

Figure generate3DLsystem(const string& fileName, const img::Color& color){
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
    for (int iteration = 0; iteration != LSystem.get_nr_iterations(); iteration++){
        string newDrawString;
        for (auto oldchar : drawString){
            if (LSystem.get_alphabet().find(oldchar) != LSystem.get_alphabet().end()){
                newDrawString += LSystem.get_replacement(oldchar);
            }
            else{
                newDrawString += oldchar;
            }
        }
        drawString = newDrawString;
    }
    for (auto character : drawString){
        Vector3D tempH = H;
        Vector3D tempL = L;
        Vector3D tempU = U;
        if (character == '+'){
            H = tempH*cos((angle*pi)/180) + tempL*sin((angle*pi)/180);
            L = -tempH*sin((angle*pi)/180) + tempL*cos((angle*pi)/180);
        }
        else if (character == '-'){
            H = tempH*cos((-angle*pi)/180) + tempL*sin((-angle*pi)/180);
            L = -tempH*sin((-angle*pi)/180) + tempL*cos((-angle*pi)/180);
        }
        else if (character == '^'){
            H = tempH*cos((angle*pi)/180) + tempU*sin((angle*pi)/180);
            U = -tempH*sin((angle*pi)/180) + tempU*cos((angle*pi)/180);
        }
        else if (character == '&'){
            H = tempH*cos((-angle*pi)/180) + tempU*sin((-angle*pi)/180);
            U = -tempH*sin((-angle*pi)/180) + tempU*cos((-angle*pi)/180);
        }
        else if (character == '\\'){
            L = tempL*cos((angle*pi)/180) - tempU*sin((angle*pi)/180);
            U = tempL*sin((angle*pi)/180) + tempU*cos((angle*pi)/180);
        }
        else if (character == '/'){
            L = tempL*cos((-angle*pi)/180) - tempU*sin((-angle*pi)/180);
            U = tempL*sin((-angle*pi)/180) + tempU*cos((-angle*pi)/180);
        }
        else if (character == '|'){
            H = -H;
            L = -L;
        }
        else if (character == '('){
            Lsystem3DElement element(p2, H, L, U);
            stack.push_back(element);
        }
        else if (character == ')'){
            p1 = stack[stack.size()-1].position;
            p2 = stack[stack.size()-1].position;
            H = stack[stack.size()-1].H;
            L = stack[stack.size()-1].L;
            U = stack[stack.size()-1].U;
            stack.pop_back();
        }
        else {
            p2 = Vector3D::point(p2.x+H.x, p2.y+H.y, p2.z+H.z);
            if (LSystem.draw(character)) {
                points.push_back(p1);
                points.push_back(p2);
                faces.push_back(Face({(int) points.size()-2, (int) points.size()-1}));
            }
            p1 = p2;
        }
    }
    return {points, faces, color};
}

vector<Face> triangulate(const Face& face){
    vector<Face> faces;
    for (int i = 1; i < face.point_indexes.size()-1; i++){
        vector<int> indexes;
        indexes.push_back(face.point_indexes[0]);
        indexes.push_back(face.point_indexes[i]);
        indexes.push_back(face.point_indexes[i+1]);
        faces.push_back(Face(indexes));
    }
    return faces;
}

img::EasyImage generate_image(const ini::Configuration &configuration) {
    Matrix m;
    string type = configuration["General"]["type"].as_string_or_die();
    int size = configuration["General"]["size"].as_int_or_die();
    vector<double> bgColorTuple = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    img::Color bgColor(bgColorTuple[0]*255, bgColorTuple[1]*255, bgColorTuple[2]*255);
    if (type == "2DLSystem"){
        string LSystem = configuration["2DLSystem"]["inputfile"].as_string_or_die();
        vector<double> lineColorTuple = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
        img::Color lineColor(lineColorTuple[0]*255, lineColorTuple[1]*255, lineColorTuple[2]*255);
        return generateLsystem(LSystem, size, bgColor, lineColor);
    }
    else {
        Matrix transform;
        vector<double> eye = configuration["General"]["eye"].as_double_tuple_or_die();
        Vector3D eyeVector = Vector3D::point(eye[0], eye[1], eye[2]);
        double theta, phi, r;
        toPolar(eyeVector, theta, phi, r);
        Matrix eyeMat = eyeMatrix(theta, phi, r);
        int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        Figures3D figures;
        for (int i = 0; i != nrFigures; i++){
            Figure fig;
            string figure = "Figure";
            figure += to_string(i);
            vector<double> color = configuration[figure]["color"].as_double_tuple_or_die();
            img::Color figColor(color[0] * 255, color[1] * 255, color[2] * 255);
            double rotate_x_angle = configuration[figure]["rotateX"].as_double_or_die();
            double rotate_y_angle = configuration[figure]["rotateY"].as_double_or_die();
            double rotate_z_angle = configuration[figure]["rotateZ"].as_double_or_die();
            vector<double> center = configuration[figure]["center"].as_double_tuple_or_die();
            double scale = configuration[figure]["scale"].as_double_or_die();
            string figType = configuration[figure]["type"].as_string_or_die();
            Vector3D translationVector = Vector3D::point(center[0], center[1], center[2]);
            Matrix xRot = rotateX(rotate_x_angle);
            Matrix yRot = rotateY(rotate_y_angle);
            Matrix zRot = rotateZ(rotate_z_angle);
            Matrix translation = translate(translationVector);
            Matrix scaleMatrix = scaleFigure(scale);
            transform = scaleMatrix * xRot * yRot * zRot * translation;
            fig.color = figColor;
            if (figType == "LineDrawing") {
                int nrPoints = configuration[figure]["nrPoints"].as_int_or_die();
                for (int j = 0; j != nrPoints; j++) {
                    string point = "point";
                    point += to_string(j);
                    vector<double> coPoint = configuration[figure][point].as_double_tuple_or_die();
                    Vector3D p = Vector3D::point(coPoint[0], coPoint[1], coPoint[2]);
                    fig.points.push_back(p);
                }
                int nrLines = configuration[figure]["nrLines"].as_int_or_die();
                for (int j = 0; j != nrLines; j++) {
                    string line = "line";
                    line += to_string(j);
                    vector<int> indexLine = configuration[figure][line].as_int_tuple_or_die();
                    Face face(indexLine);
                    fig.faces.push_back(face);
                }
            }
            else if (figType == "Cube"){
                fig = createCube((figColor));
            }
            else if (figType == "Tetrahedron"){
                fig = createTetrahedron(figColor);
            }
            else if (figType == "Octahedron"){
                fig = createOctahedron(figColor);
            }
            else if (figType == "Icosahedron"){
                fig = createIcosahedron(figColor);
            }
            else if (figType == "Dodecahedron"){
                fig = createDodecahedron(figColor);
            }
            else if (figType == "Cylinder" or figType == "Cone"){
                int n = configuration[figure]["n"].as_int_or_die();
                double h = configuration[figure]["height"].as_double_or_die();
                if (figType == "Cylinder"){
                    fig = createCylinder(n, h, figColor);
                }
                else{
                    fig = createCone(n, h, figColor);
                }
            }
            else if (figType == "Sphere"){
                int n = configuration[figure]["n"].as_int_or_die();
                fig = createSphere(69, n, figColor);
            }
            else if (figType == "Torus"){
                int n = configuration[figure]["n"].as_int_or_die();
                double R = configuration[figure]["R"].as_double_or_die();
                double radius = configuration[figure]["r"].as_double_or_die();
                int temp = configuration[figure]["m"].as_int_or_die();
                fig = createDonut(radius, R, n, temp, figColor);
            }
            else if (figType == "3DLSystem"){
                string filename = configuration[figure]["inputfile"].as_string_or_die();
                fig = generate3DLsystem(filename, figColor);
            }
            if (type == "ZBuffering" and figType != "LineDrawing" and figType != "3DLSystem"){
                vector<Face> faces;
                for (const auto& face : fig.faces){
                    for (const auto& f : triangulate(face)){
                        faces.push_back(f);
                    }
                }
                fig.faces = faces;
            }
            applyTransformation(fig, transform);
            figures.push_back(fig);
        }
        applyTransformation(figures, eyeMat);
        Lines2D lines = doProjection(figures);
        return draw2DLines(figures, lines, size, type, bgColor);
    }
    return {};
}

int main(int argc, char const* argv[])
{
        int retVal = 0;
        try
        {
                std::vector<std::string> args = std::vector<std::string>(argv+1, argv+argc);
                if (args.empty()) {
                        std::ifstream fileIn("filelist");
                        std::string filelistName;
                        while (std::getline(fileIn, filelistName)) {
                                args.push_back(filelistName);
                        }
                }
                for(std::string fileName : args)
                {
                        ini::Configuration conf;
                        try
                        {
                                std::ifstream fin(fileName);
                                fin >> conf;
                                fin.close();
                        }
                        catch(ini::ParseException& ex)
                        {
                                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                                retVal = 1;
                                continue;
                        }

                        img::EasyImage image = generate_image(conf);
                        if(image.get_height() > 0 && image.get_width() > 0)
                        {
                                std::string::size_type pos = fileName.rfind('.');
                                if(pos == std::string::npos)
                                {
                                        //filename does not contain a '.' --> append a '.bmp' suffix
                                        fileName += ".bmp";
                                }
                                else
                                {
                                        fileName = fileName.substr(0,pos) + ".bmp";
                                }
                                try
                                {
                                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;

                                }
                                catch(std::exception& ex)
                                {
                                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                                        retVal = 1;
                                }
                        }
                        else
                        {
                                std::cout << "Could not generate image for " << fileName << std::endl;
                        }
                }
        }
        catch(const std::bad_alloc &exception)
        {
    		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
    		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
    		//(Unless of course you are already consuming the maximum allowed amount of memory)
    		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
		//mark the test as failed while in reality it just needed a bit more memory
                std::cerr << "Error: insufficient memory" << std::endl;
                retVal = 100;
        }
        return retVal;
}
