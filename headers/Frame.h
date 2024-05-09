#pragma once

#include "structure.h"

#define DIMMINSION 3

struct Node
{
    double pos[DIMMINSION];           //  position tensor
    double dis[DIMMINSION]{0};        //  displacment
    double nodalForce[DIMMINSION]{0}; //  nodal Force
    bool support[DIMMINSION]{0};

public:
    Node(double x, double y);
    void setDis(UCHAR i, double value);
    void setNodalForce(UCHAR i, double value);
    void print();
};

struct Member
{
    size_t iNode, jNode;
    double IendForces[DIMMINSION]{0}, JendForces[DIMMINSION]{0};
    double IFixedForces[DIMMINSION]{0}, JFixedForces[DIMMINSION]{0};

    Ten K11, K12, K21, K22, TransformationMatrix;
    
    double l, E, A, I;

    Member(size_t iNode, size_t jNode, double E, double A, double I);
    Member(Member &) = default;
    Member(Member &&) = default;

    ~Member();

public:
    Ten &getGlopalK11();
    Ten &getGlopalK12();
    Ten &getGlopalK21();
    Ten &getGlopalK22();

public:
    size_t getIndex_I();
    size_t getIndex_J();

   

public:
    void print();
};

class Frame
{

private:
    std::vector<Node> nodes;
    std::vector<Member> members;

public:
    Frame();

public:
    size_t addNode(double x, double y);

    void addnodeForce(size_t nodeIndex, size_t forceIndex, double force);
    void addSupport(size_t nodeIndex, size_t disIndex, double dis);

    size_t addMember(size_t i, size_t j, double E, double A, double I);

private:
    void buildMember(size_t index);

public:
    void solve();
    void print();
};
