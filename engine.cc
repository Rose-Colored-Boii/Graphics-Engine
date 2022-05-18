#include "easy_image.h"
#include "Light.h"
#include "ini_configuration.h"
#include "vector3d.h"
#include "Figure.h"
#include <iostream>
#include <string>
#include <algorithm>
#include "ZBuffer.h"
#include "Create_Figures.h"

img::EasyImage generate_image(const ini::Configuration &configuration) {
    Matrix m;
    string type = configuration["General"]["type"].as_string_or_die();
    int size = configuration["General"]["size"].as_int_or_die();
    vector<double> bgColorTuple = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    img::Color bgColor(bgColorTuple[0] * 255, bgColorTuple[1] * 255, bgColorTuple[2] * 255);
    if (type == "2DLSystem") {
        string LSystem = configuration["2DLSystem"]["inputfile"].as_string_or_die();
        vector<double> lineColorTuple = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
        img::Color lineColor(lineColorTuple[0] * 255, lineColorTuple[1] * 255, lineColorTuple[2] * 255);
        return generateLsystem(LSystem, size, bgColor, lineColor);
    } else {
        Lights3D lights;
        if (type != "LightedZBuffering") {
            lights.push_back({});
        }
        if (type == "LightedZBuffering") {
            int nrLights = configuration["General"]["nrLights"].as_int_or_die();
            for (int i = 0; i < nrLights; i++) {
                vector<double> def = {0, 0, 0};
                string light = "Light" + to_string(i);
                bool infinity = configuration[light]["infinity"].as_bool_or_default(false);
                vector<double> ambientLight = configuration[light]["ambientLight"].as_double_tuple_or_default(def);
                vector<double> diffuseLight = configuration[light]["diffuseLight"].as_double_tuple_or_default(def);
                vector<double> specularLight = configuration[light]["specularLight"].as_double_tuple_or_default(def);
                Light l;
                if (infinity) {
                    l.ambientLight = ambientLight;
                    l.diffuseLight = diffuseLight;
                    l.specularLight = specularLight;
                    vector<double> direction = configuration[light]["direction"].as_double_tuple_or_die();
                    l.direction = Vector3D::vector(direction[0], direction[1], direction[2]);
                    l.infinity = true;
                } else {
                    l.ambientLight = ambientLight;
                    l.diffuseLight = diffuseLight;
                    l.specularLight = specularLight;
                    vector<double> location = configuration[light]["location"].as_double_tuple_or_default(def);
                    l.location = Vector3D::point(location[0], location[1], location[2]);
                    double angle = configuration[light]["spotAngle"].as_double_or_default(90);
                    l.angle = angle;
                }
                lights.push_back(l);
            }
        }
        Matrix transform;
        Matrix eyeMat;
        double theta, phi, r, a, b, c;
        vector<double> eye = configuration["General"]["eye"].as_double_tuple_or_die();
        Vector3D eyeVector = Vector3D::point(eye[0], eye[1], eye[2]);
        toPolar(eyeVector, theta, phi, r);
        bool clipping = configuration["General"]["clipping"].as_bool_or_default(false);
        bool textureMapping = configuration["General"]["textureMapping"].as_bool_or_default(false);
        if (clipping) {
            vector<double> viewDirection = configuration["General"]["viewDirection"].as_double_tuple_or_die();
            Vector3D viewVector = Vector3D::vector(-viewDirection[0], -viewDirection[1], -viewDirection[2]);
            toPolar(viewVector, a, b, c);
            eyeMat = eyeMatrix(a, b, r);
        } else {
            eyeMat = eyeMatrix(theta, phi, r);
        }
        int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        vector<string> textureNames(nrFigures, "");
        Figures3D figures;
        for (int i = 0; i != nrFigures; i++) {
            Figure fig;
            Figures3D figs;
            string figure = "Figure";
            figure += to_string(i);
            double z = configuration[figure]["radius"].as_double_or_default(0);
            int x = configuration[figure]["n"].as_int_or_default(0);
            int y = configuration[figure]["m"].as_int_or_default(0);
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
            vector<double> color;
            img::Color figColor;
            if (type != "LightedZBuffering") {
                color = configuration[figure]["color"].as_double_tuple_or_die();
                figColor = img::Color(color[0] * 255, color[1] * 255, color[2] * 255);
                fig.color = figColor;
            }
            if (figType == "LineDrawing" or figType == "ThickLineDrawing") {
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
                if (figType == "ThickLineDrawing") {
                    generateThiccccccFigure(fig, figs, z, x, y, figColor);
                }
            } else if (figType == "Cube") {
                fig = createCube((figColor));
            } else if (figType == "Tetrahedron") {
                fig = createTetrahedron(figColor);
            } else if (figType == "Octahedron") {
                fig = createOctahedron(figColor);
            } else if (figType == "Icosahedron") {
                fig = createIcosahedron(figColor);
            } else if (figType == "Dodecahedron") {
                fig = createDodecahedron(figColor);
            } else if (figType == "Cylinder" or figType == "Cone") {
                int n = configuration[figure]["n"].as_int_or_die();
                double h = configuration[figure]["height"].as_double_or_die();
                if (figType == "Cylinder") {
                    fig = createCylinder(n, h, figColor);
                } else {
                    fig = createCone(n, h, figColor);
                }
            } else if (figType == "Sphere") {
                int n = configuration[figure]["n"].as_int_or_die();
                fig = createSphere(69, n, figColor);
            } else if (figType == "Torus") {
                int n = configuration[figure]["n"].as_int_or_die();
                double R = configuration[figure]["R"].as_double_or_die();
                double radius = configuration[figure]["r"].as_double_or_die();
                int temp = configuration[figure]["m"].as_int_or_die();
                fig = createDonut(radius, R, n, temp, figColor);
            } else if (figType == "3DLSystem") {
                string filename = configuration[figure]["inputfile"].as_string_or_die();
                fig = generate3DLsystem(filename, figColor);
            } else if (figType == "FractalCube") {
                int n = configuration[figure]["nrIterations"].as_int_or_die();
                double s = configuration[figure]["fractalScale"].as_double_or_die();
                fig = createCube(figColor);
                for (const auto &f: generateFractal(fig, n, s)) {
                    figs.push_back(f);
                }
            } else if (figType == "FractalDodecahedron") {
                int n = configuration[figure]["nrIterations"].as_int_or_die();
                double s = configuration[figure]["fractalScale"].as_double_or_die();
                fig = createDodecahedron(figColor);
                for (const auto &f: generateFractal(fig, n, s)) {
                    figs.push_back(f);
                }
            } else if (figType == "FractalIcosahedron") {
                int n = configuration[figure]["nrIterations"].as_int_or_die();
                double s = configuration[figure]["fractalScale"].as_double_or_die();
                fig = createIcosahedron(figColor);
                for (const auto &f: generateFractal(fig, n, s)) {
                    figs.push_back(f);
                }
            } else if (figType == "FractalOctahedron") {
                int n = configuration[figure]["nrIterations"].as_int_or_die();
                double s = configuration[figure]["fractalScale"].as_double_or_die();
                fig = createOctahedron(figColor);
                for (const auto &f: generateFractal(fig, n, s)) {
                    figs.push_back(f);
                }
            } else if (figType == "FractalTetrahedron") {
                int n = configuration[figure]["nrIterations"].as_int_or_die();
                double s = configuration[figure]["fractalScale"].as_double_or_die();
                fig = createTetrahedron(figColor);
                for (const auto &f: generateFractal(fig, n, s)) {
                    figs.push_back(f);
                }
            } else if (figType == "BuckyBall" or figType == "FractalBuckyBall" or figType == "MengerSponge") {
                fig = createIcosahedron(figColor);
            } else if (figType == "ThickCube") {
                fig = createCube(figColor);
                generateThiccccccFigure(fig, figs, z, x, y, figColor);
            } else if (figType == "ThickDodecahedron") {
                fig = createDodecahedron(figColor);
                generateThiccccccFigure(fig, figs, z, x, y, figColor);
            } else if (figType == "ThickIcosahedron" or figType == "ThickBuckyBall") {
                fig = createIcosahedron(figColor);
                generateThiccccccFigure(fig, figs, z, x, y, figColor);
            } else if (figType == "ThickOctahedron") {
                fig = createOctahedron(figColor);
                generateThiccccccFigure(fig, figs, z, x, y, figColor);
            } else if (figType == "ThickTetrahedron") {
                fig = createTetrahedron(figColor);
                generateThiccccccFigure(fig, figs, z, x, y, figColor);
            } else if (figType == "Thick3DLSystem") {
                string filename = configuration[figure]["inputfile"].as_string_or_die();
                fig = generate3DLsystem(filename, figColor);
                generateThiccccccFigure(fig, figs, z, x, y, figColor);
            }
            if (textureMapping) {
                Texture texture;
                texture.fileName = configuration[figure]["texture"].as_string_or_default("");
                vector<double> vectorA = configuration[figure]["a"].as_double_tuple_or_default({0, 0, 0});
                vector<double> vectorB = configuration[figure]["b"].as_double_tuple_or_default({0, 0, 0});
                vector<double> vectorP = configuration[figure]["p"].as_double_tuple_or_default({0, 0, 0});
                texture.a = Vector3D::vector(vectorA[0], vectorA[1], vectorA[2]);
                texture.b = Vector3D::vector(vectorB[0], vectorB[1], vectorB[2]);
                texture.p = Vector3D::vector(vectorP[0], vectorP[1], vectorP[2]);
                fig.texture = texture;
            }
            if (type == "LightedZBuffering") {
                vector<double> def = {0, 0, 0};
                vector<double> ambientReflection = configuration[figure]["ambientReflection"].as_double_tuple_or_default(
                        def);
                vector<double> diffuseReflection = configuration[figure]["diffuseReflection"].as_double_tuple_or_default(
                        def);
                vector<double> specularReflection = configuration[figure]["specularReflection"].as_double_tuple_or_default(
                        def);
                double reflectionCoefficient = configuration[figure]["reflectionCoefficient"].as_double_or_default(0);
                if (figs.empty()) {
                    fig.ambientReflection = ambientReflection;
                    fig.diffuseReflection = diffuseReflection;
                    fig.specularReflection = specularReflection;
                    fig.reflectionCoefficient = reflectionCoefficient;
                } else {
                    for (auto &f: figs) {
                        f.ambientReflection = ambientReflection;
                        f.diffuseReflection = diffuseReflection;
                        f.specularReflection = specularReflection;
                        f.reflectionCoefficient = reflectionCoefficient;
                    }
                }
            } else {
                if (figs.empty()) {
                    fig.ambientReflection = color;
                } else {
                    for (auto &f: figs) {
                        f.ambientReflection = color;
                    }
                }
            }
            if ((type == "ZBuffering" or type == "LightedZBuffering") and figType != "LineDrawing" and
                figType != "3DLSystem") {
                if (figs.empty()) {
                    vector<Face> faces;
                    for (const auto &face: fig.faces) {
                        for (const auto &f: triangulate(face)) {
                            faces.push_back(f);
                        }
                    }
                    fig.faces = faces;
                } else {
                    for (auto &temp: figs) {
                        vector<Face> faces;
                        for (const auto &face: temp.faces) {
                            for (const auto &f: triangulate(face)) {
                                faces.push_back(f);
                            }
                        }
                        temp.faces = faces;
                    }
                }
            }
            if (figs.empty()) {
                applyTransformation(fig, transform);
                figures.push_back(fig);
            } else {
                applyTransformation(figs, transform);
                for (const auto &f: figs) {
                    figures.push_back(f);
                }
                figs = {};
            }
        }
        bool shadowing = configuration["General"]["shadowEnabled"].as_bool_or_default(false);
        if (type == "LightedZBuffering" and shadowing) {
            for (auto &l: lights) {
                if (!l.infinity) {
                    double x;
                    double y;
                    double z;
                    toPolar(l.location, x, y, z);
                    Matrix lightEye = eyeMatrix(x, y, z);
                    int shadowMask = configuration["General"]["shadowMask"].as_int_or_die();
                    l.eye = lightEye;
                    Figures3D temp = figures;
                    generateShadowMask(l, temp, shadowMask, clipping);
                }
            }
        }
        applyTransformation(figures, eyeMat);
        for (auto &l: lights) {
            l.direction *= eyeMat;
            l.location *= eyeMat;
        }
        if (clipping) {
            double dNear = configuration["General"]["dNear"].as_double_or_die();
            double dFar = configuration["General"]["dFar"].as_double_or_die();
            double hfov = configuration["General"]["hfov"].as_double_or_die();
            double aspectRatio = configuration["General"]["aspectRatio"].as_double_or_die();
            clip(figures, dNear, dFar, hfov, aspectRatio);
        }
        Lines2D lines = doProjection(figures);
        return draw2DLines(figures, lines, size, type, lights, clipping, eyeMat, shadowing, bgColor);
    }
}

int main(int argc, char const *argv[]) {
    int retVal = 0;
    try {
        std::vector<std::string> args = std::vector<std::string>(argv + 1, argv + argc);
        if (args.empty()) {
            std::ifstream fileIn("filelist");
            std::string filelistName;
            while (std::getline(fileIn, filelistName)) {
                args.push_back(filelistName);
            }
        }
        for (std::string fileName: args) {
            ini::Configuration conf;
            try {
                std::ifstream fin(fileName);
                fin >> conf;
                fin.close();
            }
            catch (ini::ParseException &ex) {
                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                retVal = 1;
                continue;
            }

            img::EasyImage image = generate_image(conf);
            if (image.get_height() > 0 && image.get_width() > 0) {
                std::string::size_type pos = fileName.rfind('.');
                if (pos == std::string::npos) {
                    //filename does not contain a '.' --> append a '.bmp' suffix
                    fileName += ".bmp";
                } else {
                    fileName = fileName.substr(0, pos) + ".bmp";
                }
                try {
                    std::ofstream f_out(fileName.c_str(), std::ios::trunc | std::ios::out | std::ios::binary);
                    f_out << image;

                }
                catch (std::exception &ex) {
                    std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                    retVal = 1;
                }
            } else {
                std::cout << "Could not generate image for " << fileName << std::endl;
            }
        }
    }
    catch (const std::bad_alloc &exception) {
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
