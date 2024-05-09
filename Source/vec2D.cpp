#include "../headers/vec2D.h"


void vec2D::operator=(vec2D &vec2)
{
    v[0] = vec2[0];
    v[1] = vec2[1];
}

vec2D vec2D::operator+(vec2D &vec2)
{
    vec2D result;
    result.v[0] = v[0] + vec2.v[0];
    result.v[1] = v[1] + vec2.v[1];
    return result;
}

vec2D vec2D::operator*(double scaler)
{
    vec2D result;
    result.v[0] = v[0] * scaler;
    result.v[1] = v[1] * scaler;
    return result;
}

vec2D vec2D::operator-(vec2D &vec2)
{
    vec2D result;
    result.v[0] = v[0] - vec2.v[0];
    result.v[1] = v[1] - vec2.v[1];
    return result;
}

double vec2D::operator*(vec2D &vec2)
{
    return v[0] * vec2.v[0] + v[1] * vec2.v[1];
}

double &vec2D::operator[](unsigned char i)
{
    return v[i];
}

void vec2D::print()
{
    printf("<%f %f>\n", v[0], v[1]);
}

Matrix2D::Matrix2D() {}
Matrix2D::~Matrix2D() {}

Matrix2D Matrix2D::operator+(Matrix2D &mat2)
{
    Matrix2D result;
    result.mat[0][0] = mat[0][0] + mat2.mat[0][0];
    result.mat[0][1] = mat[0][1] + mat2.mat[0][1];
    result.mat[1][0] = mat[1][0] + mat2.mat[1][0];
    result.mat[1][1] = mat[1][1] + mat2.mat[1][1];

    return result;
}

Matrix2D Matrix2D::operator*(double scaler)
{
    Matrix2D result;
    result.mat[0][0] = mat[0][0] * scaler;
    result.mat[0][1] = mat[0][1] * scaler;
    result.mat[1][0] = mat[1][0] * scaler;
    result.mat[1][1] = mat[1][1] * scaler;

    return result;
}

Matrix2D Matrix2D::operator-(Matrix2D &mat2)
{
    Matrix2D result;
    result.mat[0][0] = mat[0][0] - mat2.mat[0][0];
    result.mat[0][1] = mat[0][1] - mat2.mat[0][1];
    result.mat[1][0] = mat[1][0] - mat2.mat[1][0];
    result.mat[1][1] = mat[1][1] - mat2.mat[1][1];

    return result;
}

Matrix2D Matrix2D::operator*(Matrix2D &mat2)
{
    Matrix2D result;
    result.mat[0][0] = mat[0][0] * mat2.mat[0][0] + mat[0][1] * mat2.mat[1][0];
    result.mat[0][1] = mat[0][0] * mat2.mat[0][1] + mat[0][1] * mat2.mat[1][1];
    result.mat[1][0] = mat[1][0] * mat2.mat[0][0] + mat[1][1] * mat2.mat[1][0];
    result.mat[1][1] = mat[1][0] * mat2.mat[0][1] + mat[1][1] * mat2.mat[1][1];

    return result;
}

vec2D Matrix2D::operator*(vec2D &vec)
{
    vec2D result;
    result[0] = mat[0][0] * vec[0] + mat[0][1] * vec[1];
    result[1] = mat[1][0] * vec[0] + mat[1][1] * vec[1];
    return result;
}

double *Matrix2D::operator[](unsigned char i)
{
    return mat[i];
}

double Matrix2D::Det()
{
    return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
}

Matrix2D Matrix2D::Inverse()
{
    double det = Det();

    Matrix2D result;
    result.mat[0][0] = mat[1][1] / det;
    result.mat[0][1] = -mat[0][1] / det;
    result.mat[1][0] = -mat[1][0] / det;
    result.mat[1][1] = mat[0][0] / det;

    return result;
}

void Matrix2D::operator=(Matrix2D &mat2)
{
    mat[0][0] = mat2[0][0];
    mat[0][1] = mat2[0][1];
    mat[1][0] = mat2[1][0];
    mat[1][1] = mat2[1][1];
}

void Matrix2D::print()
{
    printf("[%lf %lf]\n", mat[0][0], mat[0][1]);
    printf("[%lf %lf]\n", mat[1][0], mat[1][1]);
}

Matrix2D DF(vec2D (*func)(vec2D), vec2D e)
{
    Matrix2D result;
    vec2D f0, df, en;
    double de = 0.00001;
    f0 = func(e);
    // change x
    en = e;
    en[0] += de;
    df = func(en) - func(e);

    result[0][0] = df[0] / de;
    result[1][0] = df[1] / de;
    // change y

    en = e;
    en[1] += de;
    df = func(en) - func(e);

    result[0][1] = df[0] / de;
    result[1][1] = df[1] / de;

    return result;
}

vec2D solve(vec2D (*func)(vec2D), vec2D a, vec2D x0)
{
    vec2D x = x0, y;

    size_t i = 0;
    do
    {
        y = (func(x) - a);
        x = x - DF(func, x).Inverse() * y;
        // x.print();
        i++;
    } while (((y * y) > 1e-20 * (a * a)) && (i < 100));

    return x;
}

Matrix2D DF(void *obj, vec2D (*func)(void *obj, vec2D), vec2D e)
{
    Matrix2D result;
    vec2D f0, df, en;
    double de = -0.00001;

    f0 = func(obj, e);
    // change x
    en = e;
    en[0] += de;
    df = func(obj, en) - func(obj, e);

    result[0][0] = df[0] / de;
    result[1][0] = df[1] / de;
    // change y

    en = e;
    en[1] += de;
    df = func(obj, en) - func(obj, e);

    result[0][1] = df[0] / de;
    result[1][1] = df[1] / de;

    return result;
}

vec2D solve(void *obj, vec2D (*func)(void *obj, vec2D), vec2D a, vec2D x0)
{
    vec2D x = x0, y;

    size_t i = 0;
    do
    {
        y = (func(obj, x) - a);
        x = x - DF(obj, func, x).Inverse() * y;

        i++;
    } while (((y * y) > 1e-20 * (a * a)) && (i < 100));

    return x;
}

vec2D rungeKutta_4(double t, vec2D &x, vec2D  (* vectorFeild)(double , vec2D&), double dt )
{
    vec2D v1, v2, v3, v4, v;
    
    v1 = vectorFeild(t, x);
    v2 = vectorFeild(t + 0.5 * dt, x + v1 * dt * 0.5);
    v3 = vectorFeild(t + 0.5 * dt, x + v2 * dt * 0.5);
    v4 = vectorFeild(t + dt, x + v3 * dt);
    // printf("------------------------------------\n");

    v = (v1 + v2 * 2 + v3 * 2 + v4) * (1.0 / 6.0);

    return v;
}