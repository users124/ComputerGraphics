//
// Created by mohammed on 09.03.22.
//

#include "draw.h"
//#include "easy_image.h"

img::EasyImage draw2DLines(Lines2D &lines, const int size,img::Color achtergrond) {
    double ymax = lines.front().p1.getY();
    double xmax = lines.front().p1.getX();
    double ymin = lines.front().p1.getY();
    double xmin = lines.front().p1.getX();

    for (auto i: lines) {
        int as;
        xmax=max(i.p1.x,xmax);
        ymax=max(i.p1.y,ymax);
        xmin=min(i.p1.x,xmin);
        ymin=min(i.p1.y,ymin);
        xmax=max(i.p2.x,xmax);
        ymax=max(i.p2.y,ymax);
        xmin=min(i.p2.x,xmin);
        ymin=min(i.p2.y,ymin);
        as+=1;
    }


    double xrange = xmax - xmin;
    double yrange = ymax - ymin;
    double imagex = size*(xrange / (max(xrange, yrange)));

    double imagey = size*(yrange / (max(xrange, yrange)));
    double d = 0.95 * (imagex / xrange);
    double DCx = d * ((xmin + xmax) / 2);
    double DCy = d * ((ymin + ymax) / 2);
    double dx = (imagex / 2) - DCx;
    double dy = (imagey / 2) - DCy;

    img::EasyImage image(imagex,imagey,achtergrond);
    for (auto i: lines) {
        int nx1 = lround((d*i.p1.x) + dx);
        int nx2 = lround((d*i.p2.x) + dx);
        int ny1 = lround((d*i.p1.y) + dy);
        int ny2 = lround((d*i.p2.y) + dy);
        img::Color r=i.color.getcolor();
        image.draw_line(nx1,ny1,nx2,ny2,i.color.getcolor());
    }
    return image;

}