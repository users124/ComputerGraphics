//
// Created by mohammed on 31.07.22.
//

#ifndef ENGINE_DRAWZBUFF_H
#define ENGINE_DRAWZBUFF_H


#include "Line2D.h"
#include <cmath>
#include <iostream>
#include <list>
using Lines2D=std::list<Line2D>;

img::EasyImage draw2DLines1(Lines2D &lines, const int size,img::Color);


#endif //ENGINE_DRAWZBUFF_H
