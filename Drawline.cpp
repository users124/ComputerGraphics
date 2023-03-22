//
// Created by student on 02.03.22.
//
#include "Line2D.h"
#include <cmath>
#include <iostream>
#include <list>
using Lines2D=std::list<Line2D>;

img::EasyImage draw2DLines( Lines2D &lines,const int size) {
    int ymax = lines.front().p1.getY();
    int xmax = lines.front().p1.getX();
    int ymin = lines.front().p1.getY();
    int xmin = lines.front().p1.getX();
    for (auto i: lines) {
        i.color.getcolor();
        if (i.p1.getY() < i.p2.getY() and i.p2.getY() > ymax) {
            ymax = i.p2.getY();
            if (i.p1.getY() < ymin) {
                ymin = i.p1.getY();
            }
        }
        if (i.p1.getY() > i.p2.getY() and i.p1.getY() > ymax) {
            ymax = i.p1.getY();
            if (i.p2.getY() < ymin) {
                ymin = i.p2.getY();
            }
        }
        if (i.p1.getX() < i.p2.getX() and i.p2.getX() > xmax) {
            xmax = i.p2.getX();
            if (i.p1.getX() < xmin) {
                xmin = i.p1.getY();
            }
        }
        if (i.p1.getX() > i.p2.getX() and i.p1.getX() > xmax) {
            xmax = i.p1.getX();
            if (i.p2.getX() < xmin) {
                xmin = i.p2.getX();
            }
        }
    }
    double xrange = xmax - xmin;
    double yrange = ymax - ymin;
    double imagex = xrange / (max(xrange, yrange));
    double imagey = xrange / (max(xrange, yrange));
    double d = 0.95 * (imagex / xrange);
    double DCx = d * ((xmax + xmax) / 2);
    double DCy = d * ((ymax + ymax) / 2);
    double dx = (imagex / 2) * DCx;
    double dy = (imagey / 2) * DCy;
    img::EasyImage image(imagex,imagey);
    for (auto i: lines) {
        int nx1 = lround((d*i.p1.x) + dx);
        int nx2 = lround((d*i.p2.x) + dx);
        int ny1 = lround((d*i.p1.y) + dy);
        int ny2 = lround((d*i.p2.y) + dy);
        image.draw_line(nx1,ny1,nx2,ny2,i.color.getcolor());
    }
    return image;
}