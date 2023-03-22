////
//// Created by mohammed on 02.08.22.
////
//
//#include "drawzt.h"
//#include "Zbuffer.h"
//#include "Extra_functies.h"
//
//img::EasyImage draw2DLinesz(Lines2D &lines, const int size,img::Color achtergrond) {
//    double ymax = lines.front().p1.getY();
//    double xmax = lines.front().p1.getX();
//    double ymin = lines.front().p1.getY();
//    double xmin = lines.front().p1.getX();
//
//    for (auto i: lines) {
//        int as;
//        xmax=max(i.p1.x,xmax);
//        ymax=max(i.p1.y,ymax);
//        xmin=min(i.p1.x,xmin);
//        ymin=min(i.p1.y,ymin);
//        xmax=max(i.p2.x,xmax);
//        ymax=max(i.p2.y,ymax);
//        xmin=min(i.p2.x,xmin);
//        ymin=min(i.p2.y,ymin);
//        as+=1;
//    }
//
//
//
//    double xrange = xmax - xmin;
//    double yrange = ymax - ymin;
//    double imagex = size*(xrange / (max(xrange, yrange)));
//
//    double imagey = size*(yrange / (max(xrange, yrange)));
//    double d = 0.95 * (imagex / xrange);
//    double DCx = d * ((xmin + xmax) / 2);
//    double DCy = d * ((ymin + ymax) / 2);
//    double dx = (imagex / 2) - DCx;
//    double dy = (imagey / 2) - DCy;
//    ZBuffer zbuf(imagex,imagey) ;
//    img::EasyImage image(imagex,imagey,achtergrond);
//    int o = 0;
//    for (auto i: lines) {
//        int nx1 = lround((d*i.p1.x+0.5) + dx);
//        int nx2 = lround((d*i.p2.x+0.5) + dx);
//        int ny1 = lround((d*i.p1.y+0.5) + dy);
//        int ny2 = lround((d*i.p2.y+0.5) + dy);
//        img::Color r=i.color.getcolor();
//        draw_zline1(image,zbuf,nx1,ny1,nx2,ny2,r,i.z1,i.z2);
//    }
//    return image;
//
//}