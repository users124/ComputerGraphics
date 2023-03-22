//
// Created by mohammed on 09.03.22.
//

#ifndef ENGINE_DRAW_H
#define ENGINE_DRAW_H
#include "Line2D.h"
#include <cmath>
#include <iostream>
#include <list>
using Lines2D=std::list<Line2D>;

img::EasyImage draw2DLines(Lines2D &lines, const int size,img::Color);



#endif //ENGINE_DRAW_H
