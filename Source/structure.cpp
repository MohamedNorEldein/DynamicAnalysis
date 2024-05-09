
#include "../headers/structure.h"


void StructureAddStiffMat(structure* pstr, size_t i, size_t j, Ten &Aii, Ten &Aij, Ten &Aji, Ten &Ajj)
{
    pstr->K.addSubSpace(Aii, i, i);
    pstr->K.addSubSpace(Aij, i, j);
    pstr->K.addSubSpace(Aji, j, i);
    pstr->K.addSubSpace(Ajj, j, j);
}

void StructureAddForce(structure* pstr,size_t i, double f)
{
    pstr->p[i][0] += f;
}

void StructureAddDis(structure* pstr,size_t i, double dis)
{
    pstr->y[i][0] += dis;
    pstr->isSupport[i][0] = 1;
    pstr->supports++;
}

double StructureGetDis(structure* pstr,size_t i)
{
    return pstr->y[i][0];
}

void print(structure* pstr)
{
    printf(pstr->K);
}

void solve(structure* pstr)
{
    size_t DF = pstr->K.getColumnNum() - pstr->supports;

    std::vector<size_t> map_a, map_b;
    Ten K11(DF, DF), K12(DF, pstr->supports);
    Ten y1(DF, 1), y2(pstr->supports, 1);
    Ten p1(DF, 1), p2(pstr->supports, 1);

    size_t i = 0, j = 0, s = 0;

    // decompose
    for (size_t i = 0; i < pstr->K.getColumnNum(); i++)
    {
        if (pstr->isSupport[i][0] == 0)
        {
            map_a.push_back(i);
            j++;
        }
        else
        {
            map_b.push_back(i);
            s++;
        }
    }

    for (i = 0; i < DF; i++)
    {
        for (j = 0; j < DF; j++)
        {
            K11[i][j] = pstr->K[map_a[i]][map_a[j]];
        }
        y1[i][0] = pstr->y[map_a[i]][0];
        p1[i][0] = pstr->p[map_a[i]][0];
    }

    for (i = 0; i < DF; i++)
    {
        for (j = 0; j < pstr->supports; j++)
        {
            K12[i][j] = pstr->K[map_a[i]][map_b[j]];
        }
    }

    for (i = 0; i < pstr->supports; i++)
    {
        y2[i][0] = pstr->y[map_b[i]][0];
        p2[i][0] = pstr->p[map_b[i]][0];
    }

    // solve

    y1 = K11.inverse() * (p1 - K12 * y2);

    // out put
    for (i = 0; i < DF; i++)
    {
        pstr->y[map_a[i]][0] = y1[i][0];
    }
    pstr->p = pstr->K * pstr->y;

#ifdef STRUCTUR_DEBUG
    printf("free motions\n");
    for (auto a : map_a)
    {
        printf(" %d\n", a);
    }

    printf("supports \n");
    for (auto a : map_b)
    {
        printf(" %d\n", a);
    }

    printf(K);
    printf(y);

    printf("/***********************************/\n");
    printf("K11\n");
    printf(K11);
    printf("y1\n");
    printf(y1);

    printf("/***********************************/\n");
    printf("K12\n");
    printf(K12);
    printf("y2\n");
    printf(y2);

    printf("/***********************************/\n");
    printf("y\n");
    printf(y);
    printf("P\n");
    printf(p);
#endif
}
