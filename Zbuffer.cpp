//
// Created by mohammed on 30.07.22.
//

#include "Zbuffer.h"

ZBuffer::ZBuffer(const int width, const int height) {
    double posInf = numeric_limits<double>::infinity();
    for (int i = 0; i <width ; ++i) {
        vector<double>h;;
        if (i == 244) {
            int a = 0;
        }
        for (int j = 0; j <height ; ++j) {
            h.push_back(posInf);

        }
        zbuf.push_back(h);
        h.clear();
    }
}

const vector<vector<double>> &ZBuffer::getZbuf() const {
    return zbuf;
}

void ZBuffer::setZbuf(int width, int height, double value) {
    ZBuffer::zbuf[width][height] = value;
}
