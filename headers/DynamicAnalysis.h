#pragma once

#include <iostream>
#include "vec2D.h"

#ifdef CAPI
#define CAPI extern "C" __declspec(dllexport)
#else
#define CAPI extern "C" __declspec(dllimport)
#endif

struct Force
{
    double k, Qu;
    double Ep, up;
};

CAPI double forceFunc(Force *f, double y);

struct SDOF
{
    double *yg;
    double *a, *v, *y;
    double dt, C, m;
    double *Q, *u;
    Force *forces[10];
    size_t fCount, count;
    

    SDOF(double m, double C)
        : yg(nullptr), a(nullptr), v(nullptr), y(nullptr), dt(0), C(C), fCount(0), m(m), count(0)
    {
    }
};

CAPI SDOF *createSDOF(double m, double c);

CAPI int addForce(SDOF *s, Force *f);

CAPI double calcForceSDOF(SDOF *s, double y);

CAPI void setGroundDis(SDOF *s, double *ag, size_t count, double dt);
CAPI void solve(SDOF *sys);

double interPolate(double t, double dt, double *s);

CAPI double changeTimeStep(double dt_old, double *s, size_t count, double *ag, double *vg, double *yg, size_t newCount);
