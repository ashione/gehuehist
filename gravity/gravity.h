#ifndef GRAVITY_H
#define GRAVITY_H

#include <vector>
#include <cmath>

using namespace std;

struct Point {
    double x;
    double y;
    Point(double a, double b) : x(a), y(b){};
};

extern void gravity(vector<double>&,vector<double>&);

#endif
