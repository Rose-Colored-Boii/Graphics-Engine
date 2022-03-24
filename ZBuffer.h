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

    ZBuffer(const int width, const int height){
        for (int i = 0; i < height; i++){
            matrix.push_back(vector<double>(width));
            for (int j = 0; j < width; j++){
                matrix[i][j] = numeric_limits<double>::infinity();
            }
        }
    }
};


#endif //ENGINE_CC_ZBUFFER_H
