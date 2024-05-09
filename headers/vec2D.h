#pragma once
#include <iostream>

class vec2D
{
private:
    double v[2]{};

public:
    vec2D()
    {
    }
    ~vec2D()
    {
    }

public:
    void operator=(vec2D &vec2);
    vec2D operator+(vec2D &vec2);
    vec2D operator*(double scaler);
    vec2D operator-(vec2D &vec2);
    double operator*(vec2D &vec2);

    double &operator[](unsigned char i);

    void print();
};

class Matrix2D
{
private:
    double mat[2][2];

public:
    Matrix2D();
    ~Matrix2D();

public:
    void operator=(Matrix2D &mat2);

    Matrix2D operator+(Matrix2D &mat2);

    Matrix2D operator*(double scaler);

    Matrix2D operator-(Matrix2D &mat2);

    Matrix2D operator*(Matrix2D &mat2);
    vec2D operator*(vec2D &vec);

    double *operator[](unsigned char i);
    double Det();
    Matrix2D Inverse();
    void print();
};

Matrix2D DF(vec2D (*func)(vec2D), vec2D e);

vec2D solve(vec2D (*func)(vec2D), vec2D a, vec2D x0);

Matrix2D DF(void *obj, vec2D (*func)(void *obj, vec2D), vec2D e);

vec2D solve(void *obj, vec2D (*func)(void *obj, vec2D), vec2D a, vec2D x0);


vec2D rungeKutta_4(double t, vec2D &x, vec2D  (* vectorFeild)(double , vec2D&), double dt );