//
// Created by liam on 3/24/22.
//

#ifndef ENGINE_CC_ZBUFFER_H
#define ENGINE_CC_ZBUFFER_H

#include "vector"
#include "limits"

using namespace std;

class ZBuffer {
public:
    vector<vector<double>> matrix;

    ZBuffer(const int width, const int height) {
        for (int i = 0; i < width; i++) {
            matrix.push_back(vector<double>(height));
            for (int j = 0; j < height; j++) {
                matrix[i][j] = numeric_limits<double>::infinity();
            }
        }
    }

    ZBuffer() {
        matrix = {};
    }
};


#endif //ENGINE_CC_ZBUFFER_H
