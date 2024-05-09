#pragma once
#include "vec.h"
#include <vector>

typedef Tensor<double> Ten;

typedef unsigned char UCHAR;

typedef struct structure
{
    Ten y, p; // displacment , External Forces
    Ten K;    // stiffness      K =    [k11    k12]
              //                  [K21    K22]
    Tensor<bool> isSupport;
    size_t supports;

public:
    //------------STRUCTURE-----------------------
    structure(size_t DF) : y(DF, 1), p(DF, 1), K(DF, DF), isSupport(DF, 1), supports(0)
    {
        y.ZERO();
        p.ZERO();
        K.ZERO();
        isSupport.ZERO();
    }

} structure;


void StructureAddStiffMat(structure* pstr, size_t i, size_t j, Ten &Aii, Ten &Aij, Ten &Aji, Ten &Ajj);
void StructureAddForce(structure* pstr,size_t i, double f);
void StructureAddDis(structure* pstr,size_t i, double dis);

double StructureGetDis(structure* pstr,size_t i);

void print(structure* pstr);

void solve(structure* pstr);