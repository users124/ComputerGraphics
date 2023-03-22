//
// Created by mohammed on 17.05.22.
//

#ifndef ENGINE_EXTRA_FUNCTIES_H
#define ENGINE_EXTRA_FUNCTIES_H
#include "ini_configuration.h"
#include <fstream>
#include "l_parser.h"
#include "draw.h"
#include <iostream>
#include <stdexcept>
#include <string>
#include "cmath"
#include "stack"
#include "figure.h"
#include "vector3d.h"
#include "Zbuffer.h"
using namespace std;


vector<Face> triangulate(const Face& face);
void draw_zline1(img::EasyImage& image, ZBuffer &Zbuffer,Vector3D &A, Vector3D& B, Vector3D C ,double d, img::Color color, double dx, double dy);
void draw_zline(img::EasyImage& image, ZBuffer &Zbuffer, unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, img::Color color, double , double);
string correct (string s ,LParser::LSystem2D& l_system,unsigned int m);
string correct1 (string s ,LParser::LSystem3D& l_system,unsigned int m);
img::EasyImage threeD (vector<vector<double>>kleur,vector<double>achtergrondkleur, vector<vector<vector<double>>>points,int size, vector<vector<vector<double>>>,vector<vector<double>>,vector<double>,   vector<double> ,   vector<double> ,   vector<double> ,   vector<double>,double,double d);
img::EasyImage twoD ( vector<double>, vector<double>,LParser::LSystem2D,int);
Matrix scaleFigure(const double scale);
Matrix rotateX1(const double angle);
Matrix rotateY1(const double angle);
Matrix rotateZ1(const double angle);
Matrix translate(const Vector3D &vector);
void applyTransformation(figure &fig, const Matrix m);
Matrix eyePointTrans(const Vector3D &eyepoint);
void toPolar(const Figures3D&,double&,double &,double &);
Point2D doProjection(const Vector3D &,const double );
Lines2D doProjection1(const Figures3D &,double d);
Lines2D doProjection2(const Figures3D & fig,double d);
img::EasyImage threeD1(vector<vector<double>>, vector<vector<double>>,vector<double>,vector<double>,vector<double>,vector<double>,double, vector<double>, vector<double>,vector<string>,double,vector<int>,vector<double>,vector<double>,vector<double>,vector<double>,double d);
figure createTetrahedron(const vector<double>kleur);
figure createCube(vector<double>);
figure createSphere(const vector<int> n,vector<double>);
figure createOctahedron(vector<double>kleur);
figure createIcosahedron(vector<double>kleur);
figure createDodecahedron(vector<double>);
figure createCylinder(const int n,const double ,vector<double>kleur);
figure createCone(const int n,const double h,vector<double>kleur);
figure createTorus(const double r,const double R, int n, int m,vector<double>);
img::EasyImage threeDl ( vector<vector<double>>kleur,vector<double>achtergrondkleur,LParser::LSystem3D l_system,int size,   vector<double> X,   vector<double> Y,   vector<double> Z,vector<vector<double>>center,vector<double>eye,   vector<double> scale,double aantalF);
img::EasyImage threeD2(vector<vector<double>>, vector<vector<double>>,vector<double>,vector<double>,vector<double>,vector<double>,double, vector<double>, vector<double>,vector<string>,double,vector<int>,vector<double>,vector<double>,vector<double>,vector<double>,double d);

#endif //ENGINE_EXTRA_FUNCTIES_H
