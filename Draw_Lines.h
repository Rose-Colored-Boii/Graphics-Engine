//
// Created by liam on 5/11/22.
//

#ifndef ENGINE_DRAW_LINES_H
#define ENGINE_DRAW_LINES_H

#endif //ENGINE_DRAW_LINES_H

#include "ZBuffer.h"
#include "easy_image.h"
#include "cmath"
#include "Light.h"
#include "Line2D.h"
#include "Figure.h"

using Figures3D = vector<Figure>;
using Lights3D = vector<Light>;
using Lines2D = vector<Line2D>;
constexpr double pi = M_PI;

void
draw_zbuf_line(ZBuffer &zBuffer, img::EasyImage &image, unsigned int x0, unsigned int y0, double z0, unsigned int x1,
               unsigned int y1, double z1, const img::Color &color) {
    int a = 0;
    if (x0 == x1 and y0 == y1) {
        double z = 1 / z0;
        if (z < zBuffer.matrix[x0][y0]) {
            zBuffer.matrix[x0][y0] = z;
            image(x0, y0) = color;
        }
    } else if (x0 == x1) {
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++) {
            a++;
        }
        //special case for x0 == x1
        int j = 0;
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++) {
            double z = 0;
            if (y0 < y1) {
                z = (j / (double) (a - 1)) / z1 + (1 - (j / (double) (a - 1))) / z0;
            } else {
                z = (j / (double) (a - 1)) / z0 + (1 - (j / (double) (a - 1))) / z1;
            }
            if (z < zBuffer.matrix[x0][i]) {
                zBuffer.matrix[x0][i] = z;
                image(x0, i) = color;
            }
            if (j > a - 1) {
                cerr << "pain" << endl;
                exit(1);
            }
            j++;
        }
    } else if (y0 == y1) {
        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++) {
            a++;
        }
        //special case for y0 == y1
        int j = 0;
        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++) {
            double z = 0;
            if (x0 < x1) {
                z = (j / (double) (a - 1)) / z1 + (1 - (j / (double) (a - 1))) / z0;
            } else {
                z = (j / (double) (a - 1)) / z0 + (1 - (j / (double) (a - 1))) / z1;
            }
            if (z < zBuffer.matrix[i][y0]) {
                zBuffer.matrix[i][y0] = z;
                image(i, y0) = color;
            }
            j++;
        }
    } else {
        if (x0 > x1) {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(z0, z1);
            std::swap(y0, y1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0) {
            for (unsigned int i = 0; i <= (x1 - x0); i++) {
                a++;
            }
            int j = 0;
            for (unsigned int i = 0; i <= (x1 - x0); i++) {
                double z = (j / (double) (a - 1)) / z1 + (1 - (j / (double) (a - 1))) / z0;
                if (z < zBuffer.matrix[x0 + i][lround(y0 + m * i)]) {
                    zBuffer.matrix[x0 + i][lround(y0 + m * i)] = z;
                    image(x0 + i, (unsigned int) round(y0 + m * i)) = color;
                }
                j++;
            }
        } else if (m > 1.0) {
            for (unsigned int i = 0; i <= (y1 - y0); i++) {
                a++;
            }
            int j = 0;
            for (unsigned int i = 0; i <= (y1 - y0); i++) {
                double z = (j / (double) (a - 1)) / z1 + (1 - (j / (double) (a - 1))) / z0;
                if (z < zBuffer.matrix[lround(x0 + (i / m))][y0 + i]) {
                    zBuffer.matrix[lround(x0 + (i / m))][y0 + i] = z;
                    image((unsigned int) round(x0 + (i / m)), y0 + i) = color;
                }
                j++;
            }
        } else if (m < -1.0) {
            for (unsigned int i = 0; i <= (y0 - y1); i++) {
                a++;
            }
            int j = 0;
            for (unsigned int i = 0; i <= (y0 - y1); i++) {
                double z = (j / (double) (a - 1)) / z1 + (1 - (j / (double) (a - 1))) / z0;
                if (z < zBuffer.matrix[lround(x0 - (i / m))][y0 - i]) {
                    zBuffer.matrix[lround(x0 - (i / m))][y0 - i] = z;
                    image((unsigned int) round(x0 - (i / m)), y0 - i) = color;
                }
                j++;
            }
        }
    }
}

void draw_zbuf_line(ZBuffer &zBuffer, img::EasyImage &image, const unsigned int x0, const unsigned int x1,
                    const unsigned int y, const double xG, const double yG, const double zG,
                    const double dzdx, const double dzdy, const img::Color &color, Lights3D &lights, Vector3D &n,
                    const double d, const vector<double> &diffuseReflection, const vector<double> &specularReflection,
                    const double reflectionCoeff, const double dx, const double dy, Matrix &eyeMat,
                    const bool shadowing) {
    for (unsigned int i = x0; i <= x1; i++) {
        double z = zG + (i - xG) * dzdx + (y - yG) * dzdy;
        if (z < zBuffer.matrix[i][y]) {
            double diffuseRed = 0;
            double diffuseGreen = 0;
            double diffuseBlue = 0;
            double specularRed = 0;
            double specularGreen = 0;
            double specularBlue = 0;
            for (auto &light: lights) {
                double zCo = 1 / z;
                double xCo = ((i - dx) * -zCo) / d;
                double yCo = ((y - dy) * -zCo) / d;
                Vector3D point = Vector3D::point(xCo, yCo, zCo);
                double lightZ, lightZ_1;
                if (shadowing) {
                    Vector3D originalPoint = point * Matrix::inv(eyeMat);
                    Vector3D lightPoint = originalPoint * light.eye;
                    double lightX = (light.d * lightPoint.x) / (-lightPoint.z) + light.dx;
                    double lightY = (light.d * lightPoint.y) / (-lightPoint.z) + light.dy;
                    lightZ = 1 / lightPoint.z;
                    double alphaX = lightX - floor(lightX);
                    double alphaY = lightY - floor(lightY);
                    double zE_1 = (1 - alphaX) / (1 / light.shadowMask.matrix[floor(lightX)][ceil(lightY)]) +
                                  alphaX / (1 / light.shadowMask.matrix[ceil(lightX)][ceil(lightY)]);
                    double zF_1 = (1 - alphaX) / (1 / light.shadowMask.matrix[floor(lightX)][floor(lightY)]) +
                                  alphaX / (1 / light.shadowMask.matrix[ceil(lightX)][floor(lightY)]);
                    lightZ_1 = alphaY * zE_1 + (1 - alphaY) * zF_1;
                }
                if (lightZ < lightZ_1 or !shadowing or
                    (lightZ >= (lightZ_1 - pow(10, -4)) and lightZ <= (lightZ_1 + pow(10, -4)))) {
                    if (!light.infinity) {
                        Vector3D l = Vector3D::normalise(light.location - point);
                        double cosAngle = (n.x * l.x + n.y * l.y + n.z * l.z);
                        if (cosAngle > cos(light.angle * pi / 180) and cosAngle > 0) {
                            diffuseRed += (double) light.diffuseLight[0] * diffuseReflection[0] *
                                          (1 - (1 - cosAngle) / (1 - cos(light.angle * pi / 180)));
                            diffuseGreen += (double) light.diffuseLight[1] * diffuseReflection[1] *
                                            (1 - (1 - cosAngle) / (1 - cos(light.angle * pi / 180)));
                            diffuseBlue += (double) light.diffuseLight[2] * diffuseReflection[2] *
                                           (1 - (1 - cosAngle) / (1 - cos(light.angle * pi / 180)));
                        }
                    }
                    Vector3D l;
                    if (!light.infinity) {
                        l = Vector3D::normalise(light.location - point);
                    } else {
                        l = -Vector3D::normalise(light.direction);
                    }
                    double cosAngle = (n.x * l.x + n.y * l.y + n.z * l.z);
                    Vector3D r = 2 * cosAngle * n - l;
                    Vector3D normalisedPoint = Vector3D::normalise(point);
                    double cosBeta = -normalisedPoint.x * r.x - normalisedPoint.y * r.y - normalisedPoint.z * r.z;
                    if (cosBeta > 0) {
                        specularRed +=
                                (double) light.specularLight[0] * specularReflection[0] * pow(cosBeta, reflectionCoeff);
                        specularGreen +=
                                (double) light.specularLight[1] * specularReflection[1] * pow(cosBeta, reflectionCoeff);
                        specularBlue +=
                                (double) light.specularLight[2] * specularReflection[2] * pow(cosBeta, reflectionCoeff);
                    }
                }
            }
            double red = (diffuseRed + specularRed) * 255 + color.red;
            double green = (diffuseGreen + specularGreen) * 255 + color.green;
            double blue = (diffuseBlue + specularBlue) * 255 + color.blue;
            if (red > 255) {
                red = 255;
            }
            if (green > 255) {
                green = 255;
            }
            if (blue > 255) {
                blue = 255;
            }
            img::Color newColor = img::Color(round(red), round(green),
                                             round(blue));
            zBuffer.matrix[i][y] = z;
            image(i, y) = newColor;
        }
    }
}

void draw_zbuf_triag(ZBuffer &zBuffer, img::EasyImage &image, const Vector3D &A, const Vector3D &B, const Vector3D &C,
                     double d, double dx, double dy, vector<double> &ambientReflection,
                     vector<double> &diffuseReflection,
                     vector<double> &specularReflection, double reflectionCoeff, Lights3D &lights, const string &type,
                     const bool clipping, Matrix &eyeMat, const bool shadowing, const bool textureMapping) {
    Vector3D u = B - A;
    Vector3D v = C - A;
    double w1 = u.y * v.z - u.z * v.y;
    double w2 = u.z * v.x - u.x * v.z;
    double w3 = u.x * v.y - u.y * v.x;
    double k = w1 * A.x + w2 * A.y + w3 * A.z;
    Vector3D n = Vector3D::normalise(Vector3D::cross(u, v));
    if (type == "LightedZBuffering" and k > 0) {
        n = -n;
    }
    double ambientRed = 0;
    double ambientGreen = 0;
    double ambientBlue = 0;
    double diffuseRed = 0;
    double diffuseGreen = 0;
    double diffuseBlue = 0;
    for (auto &light: lights) {
        ambientRed += (double) light.ambientLight[0] * ambientReflection[0];
        ambientGreen += (double) light.ambientLight[1] * ambientReflection[1];
        ambientBlue += (double) light.ambientLight[2] * ambientReflection[2];
        if (light.infinity) {
            Vector3D l = -Vector3D::normalise(light.direction);
            double cosAngle = n.x * l.x + n.y * l.y + n.z * l.z;
            if (cosAngle > 0) {
                diffuseRed += (double) light.diffuseLight[0] * diffuseReflection[0] * cosAngle;
                diffuseGreen +=
                        (double) light.diffuseLight[1] * diffuseReflection[1] * cosAngle;
                diffuseBlue +=
                        (double) light.diffuseLight[2] * diffuseReflection[2] * cosAngle;
            }
        }
    }
    double red = (ambientRed + diffuseRed) * 255;
    double green = (ambientGreen + diffuseGreen) * 255;
    double blue = (ambientBlue + diffuseBlue) * 255;
    if (red > 255) {
        red = 255;
    }
    if (green > 255) {
        green = 255;
    }
    if (blue > 255) {
        blue = 255;
    }
    img::Color color = img::Color(round(red), round(green), round(blue));
    if ((!clipping and k < 0) or clipping) {
        Point2D newA((d * A.x) / (-A.z) + dx, (d * A.y) / (-A.z) + dy);
        Point2D newB((d * B.x) / (-B.z) + dx, (d * B.y) / (-B.z) + dy);
        Point2D newC((d * C.x) / (-C.z) + dx, (d * C.y) / (-C.z) + dy);
        Point2D G((newA.getX() + newB.getX() + newC.getX()) / 3, (newA.getY() + newB.getY() + newC.getY()) / 3);
        double zG = 1 / (3 * A.z) + 1 / (3 * B.z) + 1 / (3 * C.z);
        double dzdx = w1 / (-d * k);
        double dzdy = w2 / (-d * k);
        int yMin = lround(min(newA.getY(), min(newB.getY(), newC.getY())) + 0.5);
        int yMax = lround(max(newA.getY(), max(newB.getY(), newC.getY())) - 0.5);
        vector<Point2D> points = {newA, newB, newC};
        for (int i = yMin; i <= yMax; i++) {
            double x1 = numeric_limits<double>::infinity();
            double x2 = numeric_limits<double>::infinity();
            double x3 = numeric_limits<double>::infinity();
            double x4 = -numeric_limits<double>::infinity();
            double x5 = -numeric_limits<double>::infinity();
            double x6 = -numeric_limits<double>::infinity();
            if ((i - points[0].getY()) * (i - points[1].getY()) <= 0) {
                x1 = x4 = points[1].getX() + (points[0].getX() - points[1].getX()) *
                                             ((i - points[1].getY()) / (points[0].getY() - points[1].getY()));
            }
            if ((i - points[1].getY()) * (i - points[2].getY()) <= 0) {
                x2 = x5 = points[2].getX() + (points[1].getX() - points[2].getX()) *
                                             ((i - points[2].getY()) / (points[1].getY() - points[2].getY()));
            }
            if ((i - points[0].getY()) * (i - points[2].getY()) <= 0) {
                x3 = x6 = points[2].getX() + (points[0].getX() - points[2].getX()) *
                                             ((i - points[2].getY()) / (points[0].getY() - points[2].getY()));
            }
            int xL = lround(min(x1, min(x2, x3)) + 0.5);
            int xR = lround(max(x4, max(x5, x6)) - 0.5);
            draw_zbuf_line(zBuffer, image, xL, xR, i, G.getX(), G.getY(), zG, dzdx, dzdy, color, lights, n, d,
                           diffuseReflection, specularReflection, reflectionCoeff, dx, dy, eyeMat, shadowing);
        }
    }
}

img::EasyImage
draw2DLines(const Figures3D &figs, const Lines2D &lines, const int size, const string &type, Lights3D &lights,
            const bool clipping, Matrix &eyeMat, const bool shadowing, const bool textureMapping,
            const img::Color &bgColor = img::Color()) {
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
    img::EasyImage image(lround(imageX), lround(imageY), bgColor);
    ZBuffer zBuffer((int) lround(imageX), (int) lround(imageY));
    double d = 0.95 * (imageX) / xRange;
    double DCx = d * (xMin + xMax) / 2;
    double DCy = d * (yMin + yMax) / 2;
    double dx = (imageX / 2) - DCx;
    double dy = (imageY / 2) - DCy;
    if (type == "ZBuffering" or type == "LightedZBuffering") {
        for (auto fig: figs) {
            for (auto face: fig.faces) {
                draw_zbuf_triag(zBuffer, image, fig.points[face.point_indexes[0]], fig.points[face.point_indexes[1]],
                                fig.points[face.point_indexes[2]], d, dx, dy, fig.ambientReflection,
                                fig.diffuseReflection, fig.specularReflection, fig.reflectionCoefficient, lights, type,
                                clipping, eyeMat, shadowing, textureMapping);
            }
        }
    } else {
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