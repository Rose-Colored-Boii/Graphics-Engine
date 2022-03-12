#include "easy_image.h"
#include "ini_configuration.h"
#include "Line2D.h"
#include "l_parser.h"
#include "LSystemElement.h"
#include "vector3d.h"
#include "Figure.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include "math.h"
#define _USE_MATH_DEFINES
using namespace std;
using Lines2D = vector<Line2D>;
using Figures3D = vector<Figure>;
constexpr double pi = 3.14159265358979323846264338327950288;

img::EasyImage draw2DLines(const Lines2D& lines, const int size, const img::Color bgColor = img::Color()){
    double xMin = min(lines[0].getP1().getX(), lines[0].getP2().getX());
    double xMax = max(lines[0].getP1().getX(), lines[0].getP2().getX());
    double yMin = min(lines[0].getP1().getY(), lines[0].getP2().getY());
    double yMax = max(lines[0].getP1().getY(), lines[0].getP2().getY());
    for (auto line : lines){
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
    img::EasyImage image(imageX, imageY, bgColor);
    double d = 0.95*(imageX)/xRange;
    double DCx = d*(xMin + xMax)/2;
    double DCy = d*(yMin + yMax)/2;
    double dx = (imageX/2) - DCx;
    double dy = (imageY/2) - DCy;
    for (auto line : lines){
        img::Color color(line.getColor().red, line.getColor().green, line.getColor().blue);
        image.draw_line(round(line.getP1().getX()*d + dx), round((line.getP1().getY()*d + dy)), round((line.getP2().getX()*d + dx)), round((line.getP2().getY()*d + dy)), color);
    }
    return image;
}

img::EasyImage generateLsystem(const string& fileName, const int size, const img::Color bgColor, const img::Color lineColor){
    ifstream input(fileName);
    LParser::LSystem2D LSystem(input);
    input.close();
    Lines2D lines;
    double pi = 3.14159265358979323846264338327950288;
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
    return draw2DLines(lines, size, bgColor);
}

Matrix scaleFigure(const double scaleFactor){
    Matrix scalingMatrix;
    for (int i = 1; i != 5; i++){
        for (int j = 1; j != 5; j++){
            scalingMatrix(i, j) = 0;
            if (j == i and j != 4){
                scalingMatrix(i, j) = scaleFactor;
            }
            else if (j == i and j == 4){
                scalingMatrix(i, j) = 1;
            }
        }
    }
    return scalingMatrix;
}

Matrix rotateX(const double angle){
    Matrix rotationMatrix;
    for (int i = 1; i != 5; i++){
        for (int j = 1; j != 5; j++){
            rotationMatrix(i, j) = 0;
            if (i == j and (i == 1 or i == 4)){
                rotationMatrix(i, j) = 1;
            }
            else if (i == j and (i == 2 or i == 3)){
                rotationMatrix(i, j) = cos(angle * pi/180);
            }
        }
    }
    rotationMatrix(3, 1) = -sin(angle * pi/180);
    rotationMatrix(2, 3) = sin(angle * pi/180);
    return rotationMatrix;
}

Matrix rotateY(const double angle){
    Matrix rotationMatrix;
    for (int i = 1; i != 5; i++){
        for (int j = 1; j != 5; j++){
            rotationMatrix(i, j) = 0;
            if (i == j and (i == 2 or i == 4)){
                rotationMatrix(i, j) = 1;
            }
            else if (i == j and (i == 1 or i == 3)){
                rotationMatrix(i, j) = cos(angle * pi/180);
            }
        }
    }
    rotationMatrix(1, 3) = -sin(angle * pi/180);
    rotationMatrix(3, 1) = sin(angle * pi/180);
    return rotationMatrix;
}

Matrix rotateZ(const double angle){
    Matrix rotationMatrix;
    for (int i = 1; i != 5; i++){
        for (int j = 1; j != 5; j++){
            rotationMatrix(i, j) = 0;
            if (i == j and (i == 3 or i == 4)){
                rotationMatrix(i, j) = 1;
            }
            else if (i == j and (i == 1 or i == 2)){
                rotationMatrix(i, j) = cos(angle * pi/180);
            }
        }
    }
    rotationMatrix(2, 1) = -sin(angle * pi/180);
    rotationMatrix(1, 2) = sin(angle * pi/180);
    return rotationMatrix;
}

Matrix translate(const Vector3D& vector){
    Matrix translationMatrix;
    for (int i = 1; i != 5; i++){
        for (int j = 1; j != 5; j++){
            translationMatrix(i, j) = 0;
            if (i == j){
                translationMatrix(i, j) = 1;
            }
        }
    }
    translationMatrix(4, 1) = vector.x;
    translationMatrix(4, 2) = vector.y;
    translationMatrix(4, 3) = vector.z;
    return translationMatrix;
}

void applyTransformation(Figure& fig, const Matrix &m){
    for (auto& point : fig.points){
        point = point * m;
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
    for (int i = 1; i != 5; i++){
        for (int j = 1; j != 5; j++){
            transformationMatrix(i, j) = 0;
        }
    }
    transformationMatrix(1, 1) = -sin(theta);
    transformationMatrix(2, 1) = cos(theta);
    transformationMatrix(1, 2) = -cos(theta) * cos(phi);
    transformationMatrix(1, 3) = cos(theta) * sin(phi);
    transformationMatrix(2, 2) = -sin(theta) * cos(phi);
    transformationMatrix(2, 3) = sin(theta) * sin(phi);
    transformationMatrix(3, 2) = sin(phi);
    transformationMatrix(3, 3) = cos(phi);
    transformationMatrix(4, 3) = -r;
    transformationMatrix(4, 4) = 1;
    return transformationMatrix;
}

Lines2D doProjection(const Figures3D& figures, const double d = 1){
    Lines2D lines;
    for (auto figure : figures){
        for (const auto& face : figure.faces){
            vector<Point2D> points;
            for (auto index : face.point_indexes){
                double x = (d * figure.points[index].x)/(-figure.points[index].z);
                double y = (d * figure.points[index].y)/(-figure.points[index].z);
                Point2D point(x, y);
                points.push_back(point);
            }
            for (Point2D p1 : points){
                for (Point2D p2 : points){
                    Line2D line(p1, p2, figure.color);
                    lines.push_back(line);
                }
            }
        }
    }
    return lines;
}

img::EasyImage generate_image(const ini::Configuration &configuration) {
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
    if (type == "Wireframe"){
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
            string figType = configuration[figure]["type"].as_string_or_die();
            double rotate_x_angle = configuration[figure]["rotateX"].as_double_or_die();
            double rotate_y_angle = configuration[figure]["rotateY"].as_double_or_die();
            double rotate_z_angle = configuration[figure]["rotateZ"].as_double_or_die();
            vector<double> center = configuration[figure]["center"].as_double_tuple_or_die();
            Vector3D translationVector = Vector3D::point(center[0], center[1], center[2]);
            double scale = configuration[figure]["scale"].as_double_or_die();
            Matrix xRot = rotateX(rotate_x_angle);
            Matrix yRot = rotateY(rotate_y_angle);
            Matrix zRot = rotateZ(rotate_z_angle);
            Matrix translation = translate(translationVector);
            Matrix scaleMatrix = scaleFigure(scale);
            transform = scaleMatrix * xRot * yRot * zRot * translation;
            vector<double> color = configuration[figure]["color"].as_double_tuple_or_die();
            img::Color figColor(color[0]*255, color[1]*255, color[2]*255);
            fig.color = figColor;
            int nrPoints = configuration[figure]["nrPoints"].as_int_or_die();
            for (int j = 0; j != nrPoints; j++) {
                string point = "point";
                point += to_string(j);
                vector<double> coPoint = configuration[figure][point].as_double_tuple_or_die();
                fig.points.push_back(Vector3D::point(coPoint[0], coPoint[1], coPoint[2])*transform);
                int nrLines = configuration[figure]["nrLines"].as_int_or_die();
                for (int j = 0; j != nrLines; j++) {
                    string line = "line";
                    line += to_string(j);
                    vector<int> indexLine = configuration[figure][line].as_int_tuple_or_die();
                    Face face(indexLine);
                    fig.faces.push_back(face);
                }
            }
            figures.push_back(fig);
        }
        applyTransformation(figures, eyeMat);
        Lines2D lines = doProjection(figures);
        return draw2DLines(lines, size, bgColor);
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
