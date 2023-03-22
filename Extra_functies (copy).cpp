//
// Created by mohammed on 17.05.22.
//

#include "Extra_functies.h"
#include "drawzbuff.h"
#include "assert.h"
#include "easy_image.h"


vector<Face> triangulate(const Face& face){
    vector<Face>faces;
    for (int i = 1; i <face.point_indexes.size()-1 ; ++i) {
        Face face1;
        face1.point_indexes.push_back(face.point_indexes[0]);
        face1.point_indexes.push_back(face.point_indexes[i]);
        face1.point_indexes.push_back(face.point_indexes[i+1]);
        faces.push_back(face1);

    }
    return faces;
}




void draw_zline(img::EasyImage &image,ZBuffer &Zbuffer,unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, img::Color color,double z1,double z2)
{
    assert(x0 < image.get_width() && y0 < image.get_height());
    assert(x1 < image.get_width() && y1 < image.get_height());
    double z=0;
    double a=0;
    double p=0;
    if (x0 == x1)
    {
        //special case for x0 == x1
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++)
        {
            if (y1<y0){
                swap(y1,y0);
                swap(z2,z1);
            }

             a= y1-(y0);
             p=(double) ((y1)-i)/a;
             z=(p/z1)+((1-p)/z2);
            if (z<Zbuffer.getZbuf()[x0][i]){
                Zbuffer.setZbuf(x0,i,z);
                (image)(x0, i) = color;

            }
        }
    }
    else if (y0 == y1)
    {
        //special case for y0 == y1
        if (x1<x0){
            swap(x1,x0);
            swap(z2,z1);
        }

        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++)
        {
            a= x1-x0;
            p=(double) (max(x0, x1)-i)/a;
            z=(p/z1)+((1-p)/z2);
            if(z<Zbuffer.getZbuf()[i][y0]){
                Zbuffer.setZbuf(i,y0,z);
                (image)(i, y0) = color;

            }
        }
    }
    else
    {
        if (x0 > x1)
        {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
            swap(z1,z2);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0)
        {
            for (unsigned int i = 0; i <= (x1 - x0); i++)
            {
                a=(x1 - x0);
                p= (double) i/a;
                z=(p/z1)+((1-p)/z2);
                if(z<Zbuffer.getZbuf()[x0 + i][(unsigned int) round(y0 + m * i)]){
                    Zbuffer.setZbuf(x0+i,(unsigned int) round(y0 + m * i),z);
                    (image)(x0 + i, (unsigned int) round(y0 + m * i)) = color;

                }
            }
        }
        else if (m > 1.0)
        {
            a=(y1 - y0)-1;
            for (unsigned int i = 0; i <= (y1 - y0); i++)
            {
                p=i/a;
                z=(p/z1)+((1-p)/z2);
                if (z<Zbuffer.getZbuf()[(unsigned int) round(x0 + (i / m))][y0 + i]){
                    Zbuffer.setZbuf((int) round(x0 + (i / m)),y0 + i,z);
                    (image)((unsigned int) round(x0 + (i / m)), y0 + i) = color;
              }
            }
        }
        else if (m < -1.0)
        {
            a=(y0 - y1)-1;
            for (unsigned int i = 0; i <= (y0 - y1); i++)
            {
                p=i/a;
                z=(p/z1)+((1-p)/z2);
                if (z<Zbuffer.getZbuf()[(unsigned int) round(x0 - (i / m))][ y0 - i]){
                    Zbuffer.setZbuf(( int) round(x0 - (i / m)),y0-i,z);
                    (image)((unsigned int) round(x0 - (i / m)), y0 - i) = color;

                }
            }
        }
    }
}


void draw_zline1(img::EasyImage& image, ZBuffer &Zbuffer,Vector3D &A, Vector3D& B, Vector3D C ,double d, img::Color color, double dx, double dy)
{
    double posInf = numeric_limits<double>::infinity();
    double neginf = numeric_limits<double>::infinity() * -1;
    Point2D A1;
    A1.x=((d*A.x)/(-A.z))+dx;
    A1.y=((d*A.y)/(-A.z))+dy;

    Point2D B1;
    B1.x=((d*B.x)/(-B.z))+dx;
    B1.y=((d*B.y)/(-B.z))+dy;

    Point2D C1;
    C1.x=((d*C.x)/(-C.z))+dx;
    C1.y=((d*C.y)/(-C.z))+dy;


    double xg=(A1.x+B1.x+C1.x)/3;
    double yg=(A1.y+B1.y+C1.y)/3;
    double zg=(1/(3*A.z))+(1/(3*B.z))+(1/(3*C.z));
    Vector3D u;
    u.x=B.x-A.x;
    u.y=B.y-A.y;
    u.z=B.z-A.z;

    Vector3D v;
    v.x=C.x-A.x;
    v.y=C.y-A.y;
    v.z=C.z-A.z;

    Vector3D w;
    w.x=(u.y*v.z)-(u.z*v.y);
    w.y=(u.z*v.x)-(u.x*v.z);
    w.z=(u.x*v.y)-(u.y*v.x);

    double k=((w.x*A.x)+(w.y*A.y)+((w.z*A.z)));
    double dzdx=((w.x)/(-d*k));
    double dzdy=((w.y)/(-d*k));


    int ymax= lround(max(max(A1.y, B1.y),C1.y)+0.5);
    int ymin= lround(min(min(A1.y, B1.y),C1.y)+0.5);
    int teller=0;

    for (int yi = ymin; yi <=ymax-1 ; ++yi) {
        double XabL=posInf;
        teller++;
        double XbcL=posInf;
        double XacL=posInf;

        double XabR=neginf;
        double XbcR=neginf;
        double XacR=neginf;
        if ((yi-A1.y)*(yi-B1.y)<=0 and A1.y!=B1.y){
            XabR=B1.x+((A1.x-B1.x)*((yi-B1.y)/(A1.y-B1.y)));
            XabL=B1.x+((A1.x-B1.x)*((yi-B1.y)/(A1.y-B1.y)));
        }
        if ((yi-B1.y)*(yi-C1.y)<=0 and B1.y!=C1.y){
            XbcR=C1.x+((B1.x-C1.x)*((yi-C1.y)/(B1.y-C1.y)));
            XbcL=C1.x+((B1.x-C1.x)*((yi-C1.y)/(B1.y-C1.y)));
        }
        if ((yi-A1.y)*(yi-C1.y)<=0 and A1.y!=C1.y){
            XacR=C1.x+((A1.x-C1.x)*((yi-C1.y)/(A1.y-C1.y)));
            XacL=C1.x+((A1.x-C1.x)*((yi-C1.y)/(A1.y-C1.y)));
        }
        int Xl= round(min(min(XabL,XbcL),XacL)+0.5);
        int XR= round(max(max(XabR,XbcR),XacR)+0.5);
//        if (Xl>XR){
//            swap(Xl,XR);
//        }

        for (unsigned int i =Xl; i <= XR; i++)
        {
            double a= XR-Xl;
            double z= (1.0001*zg)+((i-xg)*dzdx)+((yi-yg)*dzdy);
            if (yi < 0) {
//                continue;
            }
            if(z<Zbuffer.getZbuf()[i][yi]){
                Zbuffer.setZbuf(i,yi,z);
            (image)(i, yi) = color;

            }
        }

    }

}


string correct (string s ,LParser::LSystem2D& l_system,unsigned int m){
    string l;
    if(m==0){
        return s;
    }
    for (auto d: s) {
        if(l_system.get_alphabet().count(d)== true){
            l += l_system.get_replacement(d);
        }
        else{
            l+=d;}
    }
    m-=1;
    l = correct(l,l_system,m);
    return l;
}

string correct1 (string s ,LParser::LSystem3D& l_system,unsigned int m){
    string l;
    if(m==0){
        return s;
    }
    for (auto d: s) {
        if(l_system.get_alphabet().count(d)== true){
            l += l_system.get_replacement(d);
        }
        else{
            l+=d;}
    }
    m-=1;
    l = correct1(l,l_system,m);
    return l;
}


img::EasyImage twoD (vector<double>kleur,vector<double>achtergrondkleur,LParser::LSystem2D l_system,int size){
    img::EasyImage image;
    stack<tuple<double,double,double>>pos;
    Color color;
    img::Color background;
    color.red = (kleur[0]);
    color.green = (kleur[1]);
    color.blue = (kleur[2]);

    background.red = lround(achtergrondkleur[0]*255);
    background.green = lround(achtergrondkleur[1]*255);
    background.blue = lround(achtergrondkleur[2]*255);

    Lines2D lijnen;
    string Juist = correct(l_system.get_initiator(),l_system,l_system.get_nr_iterations());
    double nx = 0;
    double ny = 0;
    double angle=l_system.get_starting_angle()*(M_PI/180);
    for(auto z:Juist){
        if(z == '+'){
            angle = angle + (l_system.get_angle()*(M_PI/180));
            continue;
        }
        if(z == '-'){
            angle = angle - (l_system.get_angle()*(M_PI/180));
            continue;
        }
        if (z=='('){
            pos.push(tuple<double,double,double>(nx,ny,angle));
        }
        if (z==')'){
            nx = get<0>(pos.top());
            ny = get<1>(pos.top());
            angle= get<2>(pos.top());
            pos.pop();
        }
        if(l_system.get_alphabet().count(z) == true){
            Point2D p;
            p.x=nx;
            p.y=ny;
            nx = nx + cos(angle);
            ny = ny + (1* sin(angle));
            Point2D q;
            q.x=nx;
            q.y=ny;
            Line2D l;
            l.p1=p;
            l.p2=q;
            l.color=color;
            if(l_system.draw(z)){
                lijnen.push_back(l);}
        }

    }

    image = draw2DLines(lijnen,size,background);
    return image;
}


figure createTetrahedron(vector<double> kleur){
    Color color;
    color.red = (kleur[0]);
    color.green = (kleur[1]);
    color.blue = (kleur[2]);

    vector<Face>face;
    vector<Vector3D>k;
    Vector3D  p1;
    p1.z=-1;
    p1.y=-1;
    p1.x=1;
    k.push_back(p1);
    Vector3D p2;
    p2.z=-1;
    p2.y=1;
    p2.x=-1;
    k.push_back(p2);
    Vector3D p3;
    p3.z=1;
    p3.y=1;
    p3.x=1;
    k.push_back(p3);
    Vector3D p4;
    p4.z=1;
    p4.y=-1;
    p4.x=-1;
    k.push_back(p4);

    Face l1;
    l1.point_indexes.push_back(1);
    l1.point_indexes.push_back(2);
    l1.point_indexes.push_back(3);
    face.push_back(l1);

    Face l2;
    l2.point_indexes.push_back(2);
    l2.point_indexes.push_back(4);
    l2.point_indexes.push_back(3);
    face.push_back(l2);

    Face l3;
    l3.point_indexes.push_back(1);
    l3.point_indexes.push_back(4);
    l3.point_indexes.push_back(2);
    face.push_back(l3);

    Face l4;
    l4.point_indexes.push_back(1);
    l4.point_indexes.push_back(3);
    l4.point_indexes.push_back(4);
    face.push_back(l4);

    figure fig;
    fig.points=k;
    fig.faces=face;
    fig.color=color;

    return fig;

}

figure createSphere(const int n,vector<double>kleur){
    Color color;
    color.red = (kleur[0]);
    color.green = (kleur[1]);
    color.blue = (kleur[2]);
    int teller=0;
    vector<Face>face;
    vector<Vector3D>k;
    Face d;
    vector<Face>face1;


    figure fig= createIcosahedron(kleur);

    for (int i = 0; i <n ; ++i) {
        int l = fig.faces.size();
        for (int j = 0; j <l ; ++j) {
            Vector3D f;
            Face d;
            f.x=(fig.points[fig.faces[j].point_indexes[0]-1].x + fig.points[fig.faces[j].point_indexes[2]-1].x)/2;
            f.y=(fig.points[fig.faces[j].point_indexes[0]-1].y + fig.points[fig.faces[j].point_indexes[2]-1].y)/2;
            f.z=(fig.points[fig.faces[j].point_indexes[0]-1].z + fig.points[fig.faces[j].point_indexes[2]-1].z)/2;
            fig.points.push_back(f);
//            punt tussen 1 en 3
            f.x=(fig.points[fig.faces[j].point_indexes[0]-1].x + fig.points[fig.faces[j].point_indexes[1]-1].x)/2;
            f.y=(fig.points[fig.faces[j].point_indexes[0]-1].y + fig.points[fig.faces[j].point_indexes[1]-1].y)/2;
            f.z=(fig.points[fig.faces[j].point_indexes[0]-1].z + fig.points[fig.faces[j].point_indexes[1]-1].z)/2;
            fig.points.push_back(f);
            //            punt tussen 1 en 2
            f.x=(fig.points[fig.faces[j].point_indexes[1]-1].x + fig.points[fig.faces[j].point_indexes[2]-1].x)/2;
            f.y=(fig.points[fig.faces[j].point_indexes[1]-1].y + fig.points[fig.faces[j].point_indexes[2]-1].y)/2;
            f.z=(fig.points[fig.faces[j].point_indexes[1]-1].z + fig.points[fig.faces[j].point_indexes[2]-1].z)/2;
            fig.points.push_back(f);
//            punt tussen 2 en 3

            d.point_indexes.push_back(fig.faces[j].point_indexes[0]);
            d.point_indexes.push_back(fig.points.size()-1);
            d.point_indexes.push_back(fig.points.size()-2);
            face.push_back(d);
            d.point_indexes.clear();

            d.point_indexes.push_back(fig.faces[j].point_indexes[1]);
            d.point_indexes.push_back(fig.points.size());
            d.point_indexes.push_back(fig.points.size()-1);
            face.push_back(d);
            d.point_indexes.clear();

            d.point_indexes.push_back(fig.faces[j].point_indexes[2]);
            d.point_indexes.push_back(fig.points.size()-2);
            d.point_indexes.push_back(fig.points.size());
            face.push_back(d);
            d.point_indexes.clear();

            d.point_indexes.push_back(fig.points.size()-1);
            d.point_indexes.push_back(fig.points.size());
            d.point_indexes.push_back(fig.points.size()-2);
            face.push_back(d);
            d.point_indexes.clear();

        }
        fig.faces = face;
        face.clear();
    }
    for (auto & point : fig.points) {
        point.normalise();
    }
    return fig;


}

figure createCylinder(const int n,const double h,vector<double>kleur){
    Color color;
    color.red = (kleur[0]);
    color.green = (kleur[1]);
    color.blue = (kleur[2]);

    vector<Face>face;
    vector<Vector3D>k;

    for (int i = 0; i < n; ++i) {
        double t=(2*i*M_PI);
        double x = cos((2*i*M_PI)/n);
        Vector3D pi=Vector3D::point(x, sin((2*i*M_PI)/n),0);
        k.push_back(pi);
    }

    for (int i = 0; i < n; ++i) {
        double t=(2*i*M_PI);
        double x = cos((2*i*M_PI)/n);
        Vector3D pi=Vector3D::point(x, sin((2*i*M_PI)/n),h);
        k.push_back(pi);
    }


    for (int i = 1; i < (n)+1; ++i) {
        Face f;
        if ((i+1)%(n+1)==0) {
            f.point_indexes.push_back((1));
        }
        if ((i+1)%(n+1)!=0) {
            f.point_indexes.push_back(((i+1)%(n+1)));
        }
        f.point_indexes.push_back(i);

        f.point_indexes.push_back((n)+i);
        if ((i+1)%(n+1)==0) {
            f.point_indexes.push_back((n+1));
        }
        if ((i+1)%(n+1)!=0) {
            f.point_indexes.push_back((n)+i+1);
        }
        face.push_back(f);
    }

    Face f1;
    for (int j = n; j > 0; --j) {
        f1.point_indexes.push_back(j);
    }

    Face f2;
    for (int j = 2*n; j > n; --j) {
        f2.point_indexes.push_back(j);
    }

    face.push_back(f1);
    face.push_back(f2);

    figure fig;
    fig.points=k;
    fig.faces=face;
    fig.color=color;

    return fig;

}

figure createCone(const int n,const double h,vector<double>kleur){
    Color color;
    color.red = (kleur[0]);
    color.green = (kleur[1]);
    color.blue = (kleur[2]);

    vector<Face>face;
    vector<Vector3D>k;

    for (int i = 0; i < n; ++i) {
        double t=(2*i*M_PI);
        double x = cos((2*i*M_PI)/n);
        Vector3D pi=Vector3D::point(x, sin((2*i*M_PI)/n),0);
        k.push_back(pi);
    }

    Vector3D pn= Vector3D::point(0,0,h);
    k.push_back(pn);

    for (int i = n+1; i > 1; --i) {
        Face f;
//        f.point_indexes = {i,((i+1)%n),n+1};
//        face.push_back(f);
        f.point_indexes.push_back(i);
        if ((i+1)%(n+1)==0) {
        f.point_indexes.push_back((1));
        }
        if ((i+1)%(n+1)!=0) {
            f.point_indexes.push_back(((i+1)%(n+1)));
        }
        f.point_indexes.push_back(n+1);
        face.push_back(f);
    }
    Face f1;
    for (int j = n; j > 0; --j) {
        f1.point_indexes.push_back(j);
    }

    face.push_back(f1);

    figure fig;
    fig.points=k;
    fig.faces=face;
    fig.color=color;

    return fig;



}

figure createTorus(const double r,const double R, int n, int m,vector<double>kleur){
    Color color;
    color.red = (kleur[0]);
    color.green = (kleur[1]);
    color.blue = (kleur[2]);

    vector<Face>face;
    vector<Vector3D>k;
    vector<vector<double>>i1;
    vector<double>i2;
    int g=1;
    for (int i = 0; i <n ; ++i) {
        for (int j = 0; j < m; ++j) {
            double u=((2*i*M_PI)/n);
            double v=((2*j*M_PI)/m);
            double x=(R+(r* cos((v))))* cos(u);
            double y=(R+(r* cos((v))))* sin(u);
            double z=r* sin(v);
            Vector3D p1= Vector3D::point(x,y,z);
            k.push_back(p1);
            i2.push_back(g);
            g+=1;
        }
        i1.push_back(i2);
        i2.clear();
    }
//    m=m-1;
//    n=n-1;
    for (int i = 0; i <i1.size() ; ++i) {
        for (int j = 0; j < m; ++j) {
            Face f;
            f.point_indexes.push_back(i1[i][j]);
            f.point_indexes.push_back(i1[(i+1) % n][j]);
            f.point_indexes.push_back(i1[(i + 1) % n][(j + 1) % m]);
            f.point_indexes.push_back(i1[(i)][(j + 1) % m]);
            face.push_back(f);
//            g+=1;
//        if (g!=0){
//            f.point_indexes.push_back(i);
//            f.point_indexes.push_back((i+m)%m);
//            f.point_indexes.push_back(((i)%n) + 1);
//            f.point_indexes.push_back(((i+m)%n) + 1);
//            face.push_back(f);
//            g+=1;
//        }

        }
    }
//    for (auto f4:face){
//        for(auto f3 : f4.point_indexes){
//            f3+=1;
//        }
//    }

    figure fig;
    fig.points=k;
    fig.faces=face;
    fig.color=color;

    return fig;


}

figure createDodecahedron(vector<double>kleur){
    Color color;
    color.red = (kleur[0]);
    color.green = (kleur[1]);
    color.blue = (kleur[2]);

    vector<Face>face;
    vector<Vector3D>k;
    figure th= createIcosahedron(kleur);
    for (int i = 0; i <th.faces.size() ; ++i) {
        Vector3D p1=Vector3D::point((th.points[th.faces[i].point_indexes[0]-1].x +th.points[th.faces[i].point_indexes[1]-1].x +th.points[th.faces[i].point_indexes[2]-1].x)/3,(th.points[th.faces[i].point_indexes[0]-1].y +th.points[th.faces[i].point_indexes[1]-1].y +th.points[th.faces[i].point_indexes[2]-1].y)/3,(th.points[th.faces[i].point_indexes[0]-1].z +th.points[th.faces[i].point_indexes[1]-1].z +th.points[th.faces[i].point_indexes[2]-1].z)/3);
        k.push_back(p1);
    }
    Face l1;
    l1.point_indexes.push_back(1);
    l1.point_indexes.push_back(2);
    l1.point_indexes.push_back(3);
    l1.point_indexes.push_back(4);
    l1.point_indexes.push_back(5);
    face.push_back(l1);

    Face l2;
    l2.point_indexes.push_back(1);
    l2.point_indexes.push_back(6);
    l2.point_indexes.push_back(7);
    l2.point_indexes.push_back(8);
    l2.point_indexes.push_back(2);

    face.push_back(l2);

    Face l3;
    l3.point_indexes.push_back(2);
    l3.point_indexes.push_back(8);
    l3.point_indexes.push_back(9);
    l3.point_indexes.push_back(10);
    l3.point_indexes.push_back(3);
    face.push_back(l3);

    Face l4;
    l4.point_indexes.push_back(3);
    l4.point_indexes.push_back(10);
    l4.point_indexes.push_back(11);
    l4.point_indexes.push_back(12);
    l4.point_indexes.push_back(4);
    face.push_back(l4);


    Face l5;
    l5.point_indexes.push_back(4);
    l5.point_indexes.push_back(12);
    l5.point_indexes.push_back(13);
    l5.point_indexes.push_back(14);
    l5.point_indexes.push_back(5);
    face.push_back(l5);


    Face l6;
    l6.point_indexes.push_back(5);
    l6.point_indexes.push_back(14);
    l6.point_indexes.push_back(15);
    l6.point_indexes.push_back(6);
    l6.point_indexes.push_back(1);
    face.push_back(l6);

    Face l7;
    l7.point_indexes.push_back(20);
    l7.point_indexes.push_back(19);
    l7.point_indexes.push_back(18);
    l7.point_indexes.push_back(17);
    l7.point_indexes.push_back(16);
    face.push_back(l7);

    Face l8;
    l8.point_indexes.push_back(20);
    l8.point_indexes.push_back(15);
    l8.point_indexes.push_back(14);
    l8.point_indexes.push_back(13);
    l8.point_indexes.push_back(19);
    face.push_back(l8);

    Face l9;
    l9.point_indexes.push_back(19);
    l9.point_indexes.push_back(13);
    l9.point_indexes.push_back(12);
    l9.point_indexes.push_back(11);
    l9.point_indexes.push_back(18);
    face.push_back(l9);

    Face l10;
    l10.point_indexes.push_back(18);
    l10.point_indexes.push_back(11);
    l10.point_indexes.push_back(10);
    l10.point_indexes.push_back(9);
    l10.point_indexes.push_back(17);
    face.push_back(l10);

    Face l11;
    l11.point_indexes.push_back(17);
    l11.point_indexes.push_back(9);
    l11.point_indexes.push_back(8);
    l11.point_indexes.push_back(7);
    l11.point_indexes.push_back(16);
    face.push_back(l11);

    Face l12;
    l12.point_indexes.push_back(16);
    l12.point_indexes.push_back(7);
    l12.point_indexes.push_back(6);
    l12.point_indexes.push_back(15);
    l12.point_indexes.push_back(20);
    face.push_back(l12);

    figure fig;
    fig.points=k;
    fig.faces=face;
    fig.color=color;

    return fig;


}

figure createIcosahedron(vector<double>kleur){
    Color color;
    color.red = (kleur[0]);
    color.green = (kleur[1]);
    color.blue = (kleur[2]);

    vector<Face>face;
    vector<Vector3D>k;

    Vector3D  p1=Vector3D::point(0,0, sqrt(5)/2);
    k.push_back(p1);
    for (int i = 2; i <7 ; ++i) {
        Vector3D  p2=Vector3D::point(cos((i-2)*(2*M_PI/5)), sin((i-2)*(2*M_PI/5)), 0.5);
        k.push_back(p2);
    }
    for (int i = 7; i <12 ; ++i) {
        Vector3D  p2=Vector3D::point(cos(((M_PI/5) +(i-7)*(2*M_PI/5))), sin(((M_PI/5) +(i-7)*(2*M_PI/5))), -0.5);
        k.push_back(p2);
    }
    Vector3D  p3=Vector3D::point(0,0, -sqrt(5)/2);
    k.push_back(p3);

    Face l1;
    l1.point_indexes.push_back(1);
    l1.point_indexes.push_back(2);
    l1.point_indexes.push_back(3);
    face.push_back(l1);

    Face l2;
    l2.point_indexes.push_back(1);
    l2.point_indexes.push_back(3);
    l2.point_indexes.push_back(4);

    face.push_back(l2);

    Face l3;
    l3.point_indexes.push_back(1);
    l3.point_indexes.push_back(4);
    l3.point_indexes.push_back(5);
    face.push_back(l3);

    Face l4;
    l4.point_indexes.push_back(1);
    l4.point_indexes.push_back(5);
    l4.point_indexes.push_back(6);
    face.push_back(l4);


    Face l5;
    l5.point_indexes.push_back(1);
    l5.point_indexes.push_back(6);
    l5.point_indexes.push_back(2);
    face.push_back(l5);


    Face l6;
    l6.point_indexes.push_back(2);
    l6.point_indexes.push_back(7);
    l6.point_indexes.push_back(3);
    face.push_back(l6);

    Face l7;
    l7.point_indexes.push_back(3);
    l7.point_indexes.push_back(7);
    l7.point_indexes.push_back(8);
    face.push_back(l7);

    Face l8;
    l8.point_indexes.push_back(3);
    l8.point_indexes.push_back(8);
    l8.point_indexes.push_back(4);
    face.push_back(l8);

    Face l9;
    l9.point_indexes.push_back(4);
    l9.point_indexes.push_back(8);
    l9.point_indexes.push_back(9);
    face.push_back(l9);

    Face l10;
    l10.point_indexes.push_back(4);
    l10.point_indexes.push_back(9);
    l10.point_indexes.push_back(5);
    face.push_back(l10);

    Face l11;
    l11.point_indexes.push_back(5);
    l11.point_indexes.push_back(9);
    l11.point_indexes.push_back(10);
    face.push_back(l11);

    Face l12;
    l12.point_indexes.push_back(5);
    l12.point_indexes.push_back(10);
    l12.point_indexes.push_back(6);
    face.push_back(l12);

    Face l13;
    l13.point_indexes.push_back(6);
    l13.point_indexes.push_back(10);
    l13.point_indexes.push_back(11);
    face.push_back(l13);

    Face l14;
    l14.point_indexes.push_back(6);
    l14.point_indexes.push_back(11);
    l14.point_indexes.push_back(2);
    face.push_back(l14);

    Face l15;
    l15.point_indexes.push_back(2);
    l15.point_indexes.push_back(11);
    l15.point_indexes.push_back(7);
    face.push_back(l15);

    Face l16;
    l16.point_indexes.push_back(12);
    l16.point_indexes.push_back(8);
    l16.point_indexes.push_back(7);
    face.push_back(l16);

    Face l17;
    l17.point_indexes.push_back(12);
    l17.point_indexes.push_back(9);
    l17.point_indexes.push_back(8);
    face.push_back(l17);

    Face l18;
    l18.point_indexes.push_back(12);
    l18.point_indexes.push_back(10);
    l18.point_indexes.push_back(9);
    face.push_back(l18);

    Face l19;
    l19.point_indexes.push_back(12);
    l19.point_indexes.push_back(11);
    l19.point_indexes.push_back(10);
    face.push_back(l19);

    Face l20;
    l20.point_indexes.push_back(12);
    l20.point_indexes.push_back(7);
    l20.point_indexes.push_back(11);
    face.push_back(l20);

    figure fig;
    fig.points=k;
    fig.faces=face;
    fig.color=color;

    return fig;







}

figure createOctahedron(vector<double>kleur){
    Color color;
    color.red = (kleur[0]);
    color.green = (kleur[1]);
    color.blue = (kleur[2]);

    vector<Face>face;
    vector<Vector3D>k;
    Vector3D  p1;
    p1.x=1;
    p1.y=0;
    p1.z=0;
    k.push_back(p1);
    Vector3D  p2;
    p2.x=0;
    p2.y=1;
    p2.z=0;
    k.push_back(p2);
    Vector3D  p3;
    p3.x=-1;
    p3.y=0;
    p3.z=0;
    k.push_back(p3);

    Vector3D  p4;
    p4.x=0;
    p4.y=-1;
    p4.z=0;
    k.push_back(p4);

    Vector3D  p5;
    p5.x=0;
    p5.y=0;
    p5.z=-1;
    k.push_back(p5);

    Vector3D  p6;
    p6.x=0;
    p6.y=0;
    p6.z=1;
    k.push_back(p6);

    Face l1;
    l1.point_indexes.push_back(1);
    l1.point_indexes.push_back(2);
    l1.point_indexes.push_back(6);
    face.push_back(l1);

    Face l2;
    l2.point_indexes.push_back(2);
    l2.point_indexes.push_back(3);
    l2.point_indexes.push_back(6);

    face.push_back(l2);

    Face l3;
    l3.point_indexes.push_back(3);
    l3.point_indexes.push_back(4);
    l3.point_indexes.push_back(6);
    face.push_back(l3);

    Face l4;
    l4.point_indexes.push_back(4);
    l4.point_indexes.push_back(1);
    l4.point_indexes.push_back(6);
    face.push_back(l4);


    Face l5;
    l5.point_indexes.push_back(2);
    l5.point_indexes.push_back(1);
    l5.point_indexes.push_back(5);
    face.push_back(l5);


    Face l6;
    l6.point_indexes.push_back(3);
    l6.point_indexes.push_back(2);
    l6.point_indexes.push_back(5);
    face.push_back(l6);

    Face l7;
    l7.point_indexes.push_back(4);
    l7.point_indexes.push_back(3);
    l7.point_indexes.push_back(5);
    face.push_back(l7);

    Face l8;
    l8.point_indexes.push_back(1);
    l8.point_indexes.push_back(4);
    l8.point_indexes.push_back(5);
    face.push_back(l8);

    figure fig;
    fig.points=k;
    fig.faces=face;
    fig.color=color;

    return fig;

}

figure createCube(vector<double>kleur){
    Color color;
    color.red = (kleur[0]);
    color.green = (kleur[1]);
    color.blue = (kleur[2]);

    vector<Face>face;
    vector<Vector3D>k;
    Vector3D  p1;
    p1.z=-1;
    p1.y=-1;
    p1.x=1;
    k.push_back(p1);
    Vector3D  p2;
    p2.z=-1;
    p2.y=1;
    p2.x=-1;
    k.push_back(p2);
    Vector3D  p3;
    p3.z=1;
    p3.y=1;
    p3.x=1;
    k.push_back(p3);

    Vector3D  p4;
    p4.z=1;
    p4.y=-1;
    p4.x=-1;
    k.push_back(p4);

    Vector3D  p5;
    p5.z=-1;
    p5.y=1;
    p5.x=1;
    k.push_back(p5);

    Vector3D  p6;
    p6.z=-1;
    p6.y=-1;
    p6.x=-1;
    k.push_back(p6);

    Vector3D  p7;
    p7.z=1;
    p7.y=-1;
    p7.x=1;
    k.push_back(p7);

    Vector3D  p8;
    p8.z=1;
    p8.y=1;
    p8.x=-1;
    k.push_back(p8);

    Face l1;
    l1.point_indexes.push_back(1);
    l1.point_indexes.push_back(5);
    l1.point_indexes.push_back(3);
    l1.point_indexes.push_back(7);
    face.push_back(l1);

    Face l2;
    l2.point_indexes.push_back(5);
    l2.point_indexes.push_back(2);
    l2.point_indexes.push_back(8);
    l2.point_indexes.push_back(3);

    face.push_back(l2);

    Face l3;
    l3.point_indexes.push_back(2);
    l3.point_indexes.push_back(6);
    l3.point_indexes.push_back(4);
    l3.point_indexes.push_back(8);
    face.push_back(l3);

    Face l4;
    l4.point_indexes.push_back(6);
    l4.point_indexes.push_back(1);
    l4.point_indexes.push_back(7);
    l4.point_indexes.push_back(4);
    face.push_back(l4);


    Face l5;
    l5.point_indexes.push_back(7);
    l5.point_indexes.push_back(3);
    l5.point_indexes.push_back(8);
    l5.point_indexes.push_back(4);
    face.push_back(l5);


    Face l6;
    l6.point_indexes.push_back(1);
    l6.point_indexes.push_back(6);
    l6.point_indexes.push_back(2);
    l6.point_indexes.push_back(5);
    face.push_back(l6);

    figure fig;
    fig.points=k;
    fig.faces=face;
    fig.color=color;

    return fig;



}

Matrix scaleFigure(const double scale){
    Matrix scaleM;
    for (int i = 1; i <5 ; ++i) {
        scaleM(i,i)=scale;
    }
    scaleM(4,4)=1;
    return scaleM;
}

Matrix rotateX1(const double angle){
    Matrix rotX;
    rotX(2,2)= cos(angle);
    rotX(2,3)= sin(angle);
    rotX(3,2)=-sin(angle);
    rotX(3,3)= cos(angle);
    return rotX;
}

Matrix rotateY1(const double angle){
    Matrix rotY;
    rotY(1,1)=cos(angle);
    rotY(1,3)=-sin(angle);
    rotY(3,1)=sin(angle);
    rotY(3,3)= cos(angle);
    return rotY;

}

Matrix rotateZ1(const double angle){
    Matrix rotZ;
    rotZ(1,1)= cos(angle);
    rotZ(1,2)= sin(angle);
    rotZ(2,1)=-sin(angle);
    rotZ(2,2)= cos(angle);
    return rotZ;
}

Matrix translate(const Vector3D &vector){
    Matrix verschuif;
    verschuif(4,1)=vector.x;
    verschuif(4,2)=vector.y;
    verschuif(4,3)=vector.z;
    return verschuif;
}

void applyTransformation(figure &fig, const Matrix m){
    for (int i = 0; i < fig.points.size(); i++) {
        fig.points[i]=fig.points[i]*m;
    }
}

void toPolar(const Vector3D &point,double &r,double &phi, double &theta){
    r= sqrt(pow(point.x,2)+ pow(point.y,2)+ pow(point.z,2));
    phi= acos((point.z/r));
    theta= atan2(point.y,point.x);

}

Matrix eyePointTrans(const Vector3D &eyepoint){
    Matrix eye;
    double r;
    double theta;
    double phi;
    toPolar(eyepoint,r,phi,theta);
    eye(1,1)=-sin(theta);
    eye(1,2)=-(cos(theta) * cos(phi));
    eye(1,3)=cos(theta)*sin(phi);
    eye(2,1)=cos(theta);
    eye(2,2)=-(sin(theta)*cos(phi));
    eye(2,3)= sin(theta)* sin(phi);
    eye(3,2)=sin(phi);
    eye(3,3)= cos(phi);
    eye(4,3)=-r;
    return eye;
}

Point2D doProjection(const Vector3D& point, const double d){
    double x= (d*point.x)/(-point.z);
    double y=(d*point.y)/(-point.z);
    Point2D h;
    h.x=x;
    h.y=y;
    return h;

}

Lines2D doProjection1(const Figures3D & fig,double d){
    Lines2D r;
    for (int i = 0; i < fig.size() ; ++i) {
        for (int j = 0; j < fig[i].faces.size(); ++j) {
            Line2D f;
            f.p1=doProjection(fig[i].points[fig[i].faces[j].point_indexes[0]],d);
            f.p2=doProjection(fig[i].points[fig[i].faces[j].point_indexes[1]],d);
            f.color=fig[i].color;
            r.push_back(f);
        }
    }
    return r;
}

Lines2D doProjection2(const Figures3D & fig,double e){
    Lines2D r;
    int teller=0;
    for (int i = 0; i < fig.size() ; ++i) {
        for (int j = 0; j < fig[i].faces.size(); ++j) {
            for (int k = 0; k <fig[i].faces[j].point_indexes.size(); ++k) {
                if (teller==0){
                    Line2D f;
                    int aq=fig[i].faces[j].point_indexes[k]-1;
                    f.p1=doProjection(fig[i].points[aq],e);
                    int er=fig[i].faces[j].point_indexes[fig[i].faces[j].point_indexes.size()-1]-1;
                    f.p2=doProjection(fig[i].points[er],e);
                    Vector3D f1;
                    f1=fig[i].points[er];
                    f.color=fig[i].color;
                    f.z1=fig[i].points[aq].z;
                    f.z2=fig[i].points[er].z;
                    r.push_back(f);
                    teller+=1;
                    k=-1;
                    continue;
                }
                if (k!=fig[i].faces[j].point_indexes.size()-2){
                    Line2D f;
                    int d=(fig[i].faces[j].point_indexes[k]);
                    f.p1=doProjection(fig[i].points[d-1],e);
                    int d1=fig[i].faces[j].point_indexes[k+1];
                    f.p2=doProjection(fig[i].points[d1-1],e);
                    f.color=fig[i].color;
                    f.z1=fig[i].points[d-1].z;
                    f.z2=fig[i].points[d1-1].z;
                    r.push_back(f);
                    teller+=1;
                    continue;
                }
                if (k==fig[i].faces[j].point_indexes.size()-2){
                    Line2D f;
                    int d=(fig[i].faces[j].point_indexes[k]);
                    int d1=fig[i].faces[j].point_indexes[k+1];
                    f.p1=doProjection(fig[i].points[fig[i].faces[j].point_indexes[k]-1],e);
                    f.p2=doProjection(fig[i].points[fig[i].faces[j].point_indexes[k+1]-1],e);
                    f.color=fig[i].color;
                    f.z1=fig[i].points[fig[i].faces[j].point_indexes[k]-1].z;
                    f.z2=fig[i].points[fig[i].faces[j].point_indexes[k+1]-1].z;
                    r.push_back(f);
                    teller=0;
                    break;
                }
            }
        }
    }
    return r;
}

img::EasyImage  threeD(vector<vector<double>>kleur, vector<double>achtergrondkleur, vector<vector<vector<double>>>Points,int size,  vector<vector<vector<double>>>line,vector<vector<double>>center,vector<double>eye,   vector<double> rotateX,   vector<double> rotateY,   vector<double> rotateZ,   vector<double> scale,double aantalF,double d){
    img::EasyImage image;
    stack<tuple<double,double,double>>pos;
    img::Color background;

    background.red = lround(achtergrondkleur[0]*255);
    background.green = lround(achtergrondkleur[1]*255);
    background.blue = lround(achtergrondkleur[2]*255);

    vector<Face>face;
    vector<Vector3D>k;
    Figures3D f;
    Vector3D oog = Vector3D::vector(eye[0],eye[1],eye[2]);
    Vector3D cent;

    for (int i = 0; i <aantalF ; i++) {
        cent.x=center[i][0];
        cent.y=center[i][1];
        cent.z=center[i][2];
        for (int j = 0; j <Points[i].size() ; ++j) {
            Vector3D vec;
            vec.x=Points[i][j][0];
            vec.y=Points[i][j][1];
            vec.z=Points[i][j][2];
            k.push_back(vec);
        }
        for (int l = 0; l < line[i].size(); ++l) {
            Face e;
            e.point_indexes.push_back(line[i][l][0]);
            e.point_indexes.push_back(line[i][l][1]);
            face.push_back(e);
        }
        Color color;
        color.red = (kleur[i][0]);
        color.green = (kleur[i][1]);
        color.blue = (kleur[i][2]);
        figure fig;
        fig.points=k;
        fig.color=color;
        fig.faces=face;
        k.clear();
        face.clear();
        Matrix m;
        Matrix ar;
        Matrix rot=rotateX1(rotateX[i])* rotateY1(rotateY[i])* rotateZ1(rotateZ[i]);
        ar=eyePointTrans(oog);
        m =scaleFigure(scale[i])* rot* translate(cent);
        applyTransformation(fig,m);
        applyTransformation(fig, ar);
        f.push_back(fig);
    }
    Lines2D r;
    r= doProjection1(f,d);
    image=draw2DLines(r,size,background);
    return image;


}

img::EasyImage threeD1(vector<vector<double>> kleur, vector<vector<double>>center,vector<double>rotateX,vector<double>rotateY,vector<double>rotateZ,vector<double>scale,double aantalF,  vector<double>eye,  vector<double>achtergrondkleur,vector<string> type,double size,vector<int>  n,vector<double>h,vector<double>r1,vector<double>R,vector<double>m1,double d){
    img::EasyImage image;
    stack<tuple<double,double,double>>pos;
    img::Color background;

    background.red = lround(achtergrondkleur[0]*255);
    background.green = lround(achtergrondkleur[1]*255);
    background.blue = lround(achtergrondkleur[2]*255);

    vector<Face>face;
    vector<Vector3D>k;
    Figures3D f1;
    Vector3D oog = Vector3D::vector(eye[0],eye[1],eye[2]);
    Vector3D cent;


    for (int i = 0; i <aantalF ; i++) {
        figure fig1;
        cent.x=center[i][0];
        cent.y=center[i][1];
        cent.z=center[i][2];
        if(type[i]=="Tetrahedron"){
            fig1=(createTetrahedron(kleur[i]));
        }
//        if(type[i]=="icosahedron"){
//            fig1=(createSphere(kleur[i]));
//        }
        if(type[i]=="Cube"){
            fig1=(createCube(kleur[i]));
        }

        if(type[i]=="Octahedron"){
            fig1=(createOctahedron(kleur[i]));
        }

        if(type[i]=="Icosahedron"){
            fig1=(createIcosahedron(kleur[i]));
        }

        if(type[i]=="Dodecahedron"){
            fig1=(createDodecahedron(kleur[i]));
        }

        if(type[i]=="Sphere"){
            fig1=(createSphere(n[i],kleur[i]));
        }

        if(type[i]=="Cone"){
            fig1=(createCone(n[i],h[i],kleur[i]));
        }

        if(type[i]=="Cylinder"){
            fig1=(createCylinder(n[i],h[i],kleur[i]));
        }

        if(type[i]=="Torus"){
            fig1=(createTorus(r1[i],R[i],n[i],m1[i],kleur[i]));
        }

        Matrix m;
        Matrix ar;
        Matrix rot=rotateX1(rotateX[i])* rotateY1(rotateY[i])* rotateZ1(rotateZ[i]);
        ar=eyePointTrans(oog);
        m =scaleFigure(scale[i])* rot* translate(cent);
        applyTransformation(fig1,m);
        applyTransformation(fig1, ar);
        f1.push_back(fig1);

    }
    Lines2D r;
    r= doProjection2(f1,d);
//    image=draw2DLines(r,size,background);
    image=draw2DLines1(r,size,background);
    return image;


}

img::EasyImage threeDl ( vector<vector<double>>kleur,vector<double>achtergrondkleur,LParser::LSystem3D l_system,int size,   vector<double> X,   vector<double> Y,   vector<double> Z,vector<vector<double>>center,vector<double>eye,   vector<double> scale,double aantalF){
    img::EasyImage image;
    double angle=l_system.get_angle()*(M_PI/180);
    Figures3D fig;
    Vector3D cent;

    img::Color background;
    background.red = lround(achtergrondkleur[0]*255);
    background.green = lround(achtergrondkleur[1]*255);
    background.blue = lround(achtergrondkleur[2]*255);
    Vector3D oog = Vector3D::vector(eye[0],eye[1],eye[2]);
    int teller=1;
    int teller1=0;
//    stack<tuple<double,double,double>>pos;
    for (int i = 0; i <aantalF ; ++i) {

        figure figs;
        cent.x=center[i][0];
        cent.y=center[i][1];
        cent.z=center[i][2];
        stack<tuple<Vector3D,Vector3D,Vector3D,Vector3D,int>>posi;
        bool q= false;
        figs.points.push_back(Vector3D::point(0,0,0));
        Color color;
        img::Color background;
        color.red = (kleur[i][0]);
        color.green = (kleur[i][1]);
        color.blue = (kleur[i][2]);

        background.red = lround(achtergrondkleur[0]*255);
        background.green = lround(achtergrondkleur[1]*255);
        background.blue = lround(achtergrondkleur[2]*255);

        string Juist = correct1(l_system.get_initiator(),l_system,l_system.get_nr_iterations());
        Vector3D H=Vector3D ::vector(1,0,0);
        Vector3D L=Vector3D ::vector(0,1,0);
        Vector3D U=Vector3D ::vector(0,0,1);
        Vector3D pos;
        Vector3D OH;
        Vector3D OL;
        Vector3D OU;

        for(auto z:Juist){
            if(z == '+'){
                OH=H;
                H.x=((H.x)* cos(angle))+(L.x * sin(angle));
                H.y=((H.y)* cos(angle))+(L.y *sin(angle));
                H.z=((H.z)* cos(angle))+(L.z * sin(angle));

                L.x=((-OH.x)* sin(angle))+(L.x * cos(angle));
                L.y=((-OH.y)* sin(angle))+(L.y *cos(angle));
                L.z=((-OH.z)* sin(angle))+(L.z * cos(angle));
                continue;
            }
            if(z == '-'){
                OH=H;
                H.x=((H.x)* cos(-angle))+(L.x * sin(-angle));
                H.y=((H.y)* cos(-angle))+(L.y *sin(-angle));
                H.z=((H.z)* cos(-angle))+(L.z * sin(-angle));

                L.x=((-OH.x)* sin(-angle))+(L.x * cos(-angle));
                L.y=((-OH.y)* sin(-angle))+(L.y *cos(-angle));
                L.z=((-OH.z)* sin(-angle))+(L.z * cos(-angle));
                continue;
            }
            if (z=='&'){
                OH=H;
                H.x=((H.x)* cos(-angle))+((U.x) * sin(-angle));
                H.y=((H.y)* cos(-angle))+((U.y) *sin(-angle));
                H.z=((H.z)* cos(-angle))+((U.z) * sin(-angle));

                U.x=((-OH.x)* sin(-angle))+(U.x * cos(-angle));
                U.y=((-OH.y)* sin(-angle))+(U.y *cos(-angle));
                U.z=((-OH.z)* sin(-angle))+(U.z * cos(-angle));
                continue;

            }
            if (z=='^'){
                OH=H;
                H.x=((H.x)* cos(angle))+((U.x) * sin(angle));
                H.y=((H.y)* cos(angle))+((U.y) *sin(angle));
                H.z=((H.z)* cos(angle))+((U.z) * sin(angle));

//                U.x=(U.x* cos(angle))-(H.x* sin(angle));
                double s = U.x* cos(angle);
                double b = OH.x* sin(angle);
                double bvd = OH.x;
                double c = s - b;
                U.x = c;
                U.y=((-OH.y)* sin(angle))+(U.y *cos(angle));
                U.z=((-OH.z)* sin(angle))+(U.z * cos(angle));
                continue;

            }
            if (z=='/'){
                OL=L;
                L.x=((L.x)* cos(-angle))-((U.x) * sin(-angle));
                L.y=((L.y)* cos(-angle))-((U.y) *sin(-angle));
                L.z=((L.z)* cos(-angle))-((U.z) * sin(-angle));

                U.x=(OL.x)* sin(-angle)+(U.x * cos(-angle));
                U.y=(OL.y)* sin(-angle)+(U.y *cos(-angle));
                U.z=(OL.z)* sin(-angle)+(U.z * cos(-angle));
                continue;

            }
            if (z== '\\'){
                OL=L;
                L.x=((L.x)* cos(angle))-((U.x) * sin(angle));
                L.y=((L.y)* cos(angle))-((U.y) *sin(angle));
                L.z=((L.z)* cos(angle))-((U.z) * sin(angle));

                U.x=(OL.x)* sin(angle)+(U.x * cos(angle));
                U.y=(OL.y)* sin(angle)+(U.y *cos(angle));
                U.z=(OL.z)* sin(angle)+(U.z * cos(angle));
                continue;

            }

            if (z=='|'){
                H=-H;
                L=-L;
                continue;


            }
            if (z=='('){
                if (q){
                    posi.push(tuple<Vector3D,Vector3D,Vector3D,Vector3D,int>(pos,H,L,U,teller1));
                }
                if (!q){
                    posi.push(tuple<Vector3D,Vector3D,Vector3D,Vector3D,int>(pos,H,L,U,teller));

                }
            }
            if (z==')'){
                pos = get<0>(posi.top());
                H = get<1>(posi.top());
                L = get<2>(posi.top());
                U = get<3>(posi.top());
                teller1 = get<4>(posi.top());
                 q= true;
                posi.pop();
            }
            if(l_system.get_alphabet().count(z) == true){
                Vector3D pos1;
                pos1.x=pos.x+H.x;
                pos1.y=pos.y+H.y;
                pos1.z=pos.z+H.z;
                pos=pos1;
                if (!q){
                    teller1+=1;
                }
                teller+=1;
                figs.points.push_back(pos1);
                figs.color=color;
                if(l_system.draw(z)){
                    Face face;
                    face.point_indexes.push_back(teller1);
                    face.point_indexes.push_back(teller);
                    figs.faces.push_back(face);
                    teller1=teller-1;
                    q= false;
                }
            }

        }
        Matrix m;
        Matrix ar;
        Matrix rot=rotateX1(X[i])* rotateY1(Y[i])* rotateZ1(Z[i]);
        ar=eyePointTrans(oog);
        m =scaleFigure(scale[i])* rot* translate(cent);
        applyTransformation(figs,m);
        applyTransformation(figs, ar);
        fig.push_back(figs);
    }
    Lines2D r;
    r= doProjection2(fig,1);
    image=draw2DLines1(r,size,background);
    return image;
}
img::EasyImage threeD2(vector<vector<double>> kleur, vector<vector<double>>center,vector<double>rotateX,vector<double>rotateY,vector<double>rotateZ,vector<double>scale,double aantalF,  vector<double>eye,  vector<double>achtergrondkleur,vector<string> type,double size,vector<int>  n,vector<double>h,vector<double>r1,vector<double>R,vector<double>m1,double d){
    img::EasyImage image;
    stack<tuple<double,double,double>>pos;
    img::Color background;

    background.red = lround(achtergrondkleur[0]*255);
    background.green = lround(achtergrondkleur[1]*255);
    background.blue = lround(achtergrondkleur[2]*255);

    vector<Face>face;
    vector<Vector3D>k;
    Figures3D f1;
    Vector3D oog = Vector3D::vector(eye[0],eye[1],eye[2]);
    Vector3D cent;


    for (int i = 0; i <aantalF ; i++) {
        figure fig1;
        cent.x = center[i][0];
        cent.y = center[i][1];
        cent.z = center[i][2];
        if (type[i] == "Tetrahedron") {
            fig1 = (createTetrahedron(kleur[i]));
        }
//        if(type[i]=="icosahedron"){
//            fig1=(createSphere(kleur[i]));
//        }
        if (type[i] == "Cube") {
            fig1 = (createCube(kleur[i]));
        }

        if (type[i] == "Octahedron") {
            fig1 = (createOctahedron(kleur[i]));
        }

        if (type[i] == "Icosahedron") {
            fig1 = (createIcosahedron(kleur[i]));
        }

        if (type[i] == "Dodecahedron") {
            fig1 = (createDodecahedron(kleur[i]));
        }

        if (type[i] == "Sphere") {
            fig1 = (createSphere(n[i], kleur[i]));
        }

        if (type[i] == "Cone") {
            fig1 = (createCone(n[i], h[i], kleur[i]));
        }

        if (type[i] == "Cylinder") {
            fig1 = (createCylinder(n[i], h[i], kleur[i]));
        }

        if (type[i] == "Torus") {
            fig1 = (createTorus(r1[i], R[i], n[i], m1[i], kleur[i]));
        }

        Matrix m;
        Matrix ar;
        Matrix rot = rotateX1(rotateX[i]) * rotateY1(rotateY[i]) * rotateZ1(rotateZ[i]);
        ar = eyePointTrans(oog);
        m = scaleFigure(scale[i]) * rot * translate(cent);
        applyTransformation(fig1, m);
        applyTransformation(fig1, ar);
        f1.push_back(fig1);
    }
    Lines2D lines;
    lines = doProjection2(f1,1);
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
    int ag=0;
    int az=1;
    int ad=2;
    if (ag==1 xor ad==2  and az==1){
        int te=10;
    }

    double xrange = xmax - xmin;
    double yrange = ymax - ymin;
    double imagex = size*(xrange / (max(xrange, yrange)));

    double imagey = size*(yrange / (max(xrange, yrange)));
    double d1 = 0.95 * (imagex / xrange);
    double DCx = d1 * ((xmin + xmax) / 2);
    double DCy = d1 * ((ymin + ymax) / 2);
    double dx = (imagex / 2) - DCx;
    double dy = (imagey / 2) - DCy;
    ZBuffer zbuf(imagex,imagey) ;
    image=img::EasyImage (imagex,imagey,background);
    vector<Face>zfaces;
    for (int i = 0; i <f1.size() ; ++i) {
        for(const auto& a:f1[i].faces){
            for (auto b: triangulate(a)) {
                zfaces.push_back(b);
                img::Color color;
                color.red = (kleur[i][0]*255);
                color.green = (kleur[i][1]*255);
                color.blue = (kleur[i][2]*255);
                draw_zline1(image,zbuf,f1[i].points[b.point_indexes[0]-1],f1[i].points[b.point_indexes[1]-1],f1[i].points[b.point_indexes[2]-1],d1,color,dx,dy);
            }
        }
    }
    return image;

}
//    Lines2D r;
//    r= doProjection2(f1,d);
////    image=draw2DLines(r,size,background);
//    image=draw2DLines1(r,size,background);
//    return image;


