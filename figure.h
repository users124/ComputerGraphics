//
// Created by mohammed on 23.05.22.
//

#ifndef ENGINE_FIGURE_H
#define ENGINE_FIGURE_H
#include "vector3d.h"
#include "face.h"
#include "color.h"
#include "list"



class figure
{
public:
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    Color color;
};

typedef std::vector<figure> Figures3D;


#endif //ENGINE_FIGURE_H
