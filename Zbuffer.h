//
// Created by mohammed on 30.07.22.
//

#ifndef ENGINE_ZBUFFER_H
#define ENGINE_ZBUFFER_H
#include "vector"
#include "limits"

using namespace std;
class ZBuffer: public vector<vector<double> >
{
public:
private:
    vector<vector<double> >zbuf;
public:
    const vector<vector<double>> &getZbuf() const;

    void setZbuf(int , int, double);

    ZBuffer(const int width, const int height);
};

#endif //ENGINE_ZBUFFER_H
