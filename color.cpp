//
// Created by student on 02.03.22.
//

#include "color.h"
#include "cmath"
using namespace std;
img::Color Color::getcolor() {
    img::Color r;
    r.red=lround(this->red*255);
    r.green=lround(this->green*255);
    r.blue=lround(this->blue*255);
    return r;
}

