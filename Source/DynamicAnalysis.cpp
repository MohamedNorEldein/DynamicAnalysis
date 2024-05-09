#define CAPI
#define STRUCTUR_DEBUG
#include "../headers/DynamicAnalysis.h"

double interPolate(double t, double dt, double *s)
{
    size_t i = (size_t)(t / dt);
    double p = t / dt - i;

    return s[i] * (1.0 - p) + s[i + 1] * (p);
}

CAPI double changeTimeStep(double dt_old, double *s, size_t count, double *ag, double *vg, double *yg, size_t newCount)
{
    double dt_new = dt_old * count / newCount;
    ag[0] = 0;

    for (size_t i = 1; i < newCount; i++)
    {
        ag[i] = interPolate(dt_new * i, dt_old, s);

        vg[i] = vg[i - 1] + (ag[i] + ag[i - 1]) / 2 * dt_new;
        yg[i] = yg[i - 1] + (vg[i] + vg[i - 1]) / 2 * dt_new; // (ag[i] + ag[i - 1]) / 2 * dt_new * dt_new / 3;
    }
    return dt_new;
}

CAPI double forceFunc(Force *f, double y)
{
    double q1 = f->k * (y - f->up);

    if (q1 > f->Qu)
    {
        return f->Qu;
    }
    if (q1 < -f->Qu)
    {
        return -f->Qu;
    }
    return q1;
}

double updatePlasticDisFunc(Force *f, double y, double dy)
{
    double q1 = f->k * (y - f->up);

    if ((q1 * dy > 0) && (abs(q1) >= f->Qu))
    {
        f->up += dy;
        q1 = f->Qu * q1 / abs(q1);
        f->Ep += abs(q1 * dy);
    }
    return q1;
}

CAPI int addForce(SDOF *s, Force *f)
{
    if (s->fCount == 10)
        return -1;
    // printf("force k = %f  q = %f up = %f\n", f->k, f->Qu, f->up);
    s->forces[s->fCount] = f;

    s->fCount++;
    return 0;
}

CAPI SDOF *createSDOF(double m, double c)
{
    SDOF *sdof = new SDOF(m, c);
    // printf("%x\n", sdof);
    return sdof;
}

CAPI double calcForceSDOF(SDOF *s, double y)
{
    double f = 0;

    for (size_t i = 0; i < s->fCount; i++)
    {
        f += forceFunc(s->forces[i], y);
    }
    return f;
}

void updatePlasticDisSDOF(SDOF *s, double y, double dy, size_t index)
{
    double u = y - s->yg[index - 1];
    double du = dy - (s->yg[index] - s->yg[index - 1]);
    double Q = 0;
    for (size_t i = 0; i < s->fCount; i++)
    {
        Q += updatePlasticDisFunc(s->forces[i],u, du );
    }
    s->Q[index] = Q;
    return;
}

CAPI void setGroundDis(SDOF *s, double *ag, size_t count, double dt)
{
    s->yg = ag;
    s->dt = dt;
    s->count = count;

    free(s->y);
    free(s->v);
    free(s->a);
    free(s->Q);
    free(s->u);

    s->y = (double *)malloc(count * sizeof(double));
    s->v = (double *)malloc(count * sizeof(double));
    s->a = (double *)malloc(count * sizeof(double));
    s->Q = (double *)malloc(count * sizeof(double));
    s->u = (double *)malloc(count * sizeof(double));

    s->y[0] = 0;
    s->v[0] = 0;
    s->a[0] = 0;
}

CAPI void setInitialCondition(SDOF *sys, double y0, double v0)
{

    sys->y[0] = y0;
    sys->v[0] = v0;
}

vec2D Feild(SDOF *sys, vec2D &x, size_t i)
{
    vec2D dx;
    dx[0] = x[1];
    dx[1] = -calcForceSDOF(sys, x[0] - sys->yg[i]) / sys->m - (sys->C / sys->m * x[1]);
    return dx;
}

vec2D solveModeffiedEuler(SDOF *sys, vec2D &x, size_t i)
{
    double dt = sys->dt;
    vec2D v1, v2, v;

    v1 = Feild(sys, x, i);

    v2 = Feild(sys, x + v1 * dt, i + 1);

    v = (v1 + v2) * (1.0 / 2.0);
    return v;
}

CAPI void solve(SDOF *sys)
{
    vec2D x, y;
    
    for (size_t i = 1; i < sys->count; i++)
    {
        x[0] = sys->y[i - 1];
        x[1] = sys->v[i - 1];

        y = solveModeffiedEuler(sys, x, i);
        x = x + y * sys->dt;
        updatePlasticDisSDOF(sys, x[0], x[1] * sys->dt, i);

        sys->u[i] = x[0] - sys->yg[i];

        sys->y[i] = x[0];
        sys->v[i] = x[1];
        sys->a[i] = (sys->v[i] - sys->v[i - 1]) / sys->dt;
    }
}
