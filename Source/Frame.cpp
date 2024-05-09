

#include "..\headers\Frame.h"

#define DIMMINSION 3

Node::Node(double x, double y)
{
    pos[0] = x;
    pos[1] = y;
}
void Node::setDis(UCHAR i, double value)
{
    dis[i] = value;
    support[i] = 1;
}

void Node::setNodalForce(UCHAR i, double value) { nodalForce[i] = value; }
void Node::print()
{
    printf("Node\n"

           "\tposition \tx = %f \ty = %f\n"
           "\tnodal force \tpx = %f \tpy = %f \tM = %f\n"
           "\tdisplacment \tdx = %f \tdy = %f \tdR = %f\n",

           pos[0], pos[1],
           nodalForce[0], nodalForce[1], nodalForce[2],
           dis[0], dis[1], dis[2]

    );
}

Member::Member(size_t iNode, size_t jNode, double E, double A, double I)
    : iNode(iNode), jNode(jNode),
      E(E), A(A), I(I),
      K11(DIMMINSION, DIMMINSION),
      K12(DIMMINSION, DIMMINSION),
      K21(DIMMINSION, DIMMINSION),
      K22(DIMMINSION, DIMMINSION),
      TransformationMatrix(DIMMINSION, DIMMINSION)
{
    K11.ZERO();
    K12.ZERO();
    K21.ZERO();
    K22.ZERO();


    TransformationMatrix.ZERO();
}

Member::~Member() {}

Ten &Member::getGlopalK11()
{
    return TransformationMatrix.transpose() * K11 * TransformationMatrix;
}
Ten &Member::getGlopalK12()
{
    return TransformationMatrix.transpose() * K12 * TransformationMatrix;
}
Ten &Member::getGlopalK21()
{
    return TransformationMatrix.transpose() * K21 * TransformationMatrix;
}
Ten &Member::getGlopalK22()
{
    return TransformationMatrix.transpose() * K22 * TransformationMatrix;
}

size_t Member::getIndex_I()
{
    return iNode * DIMMINSION;
}
size_t Member::getIndex_J()
{
    return jNode * DIMMINSION;
}

void Member::print()
{
    printf("member\n"
           "\tstart node %d , end node %d \n"
           "\tE = %f \tA=%f \tI=%f\n"
           "\tstart force in member \tN = %f \tQ = %f \tM = %f\n"
           "\tend force in member \tN = %f \tQ = %f \tM = %f\n",

           iNode, jNode, E, A, I,
           IendForces[0], IendForces[1], IendForces[2],
           JendForces[0], JendForces[1], JendForces[2]

    );
}

Frame::Frame()
{
}

size_t Frame::addNode(double x, double y)
{
    nodes.push_back(Node(x, y));
    return nodes.size() - 1;
}

void Frame::addnodeForce(size_t nodeIndex, size_t forceIndex, double force)
{
    nodes[nodeIndex].nodalForce[forceIndex] = force;
}
void Frame::addSupport(size_t nodeIndex, size_t disIndex, double dis)
{
    nodes[nodeIndex].dis[disIndex] = dis;
    nodes[nodeIndex].support[disIndex] = 1;
}

size_t Frame::addMember(size_t i, size_t j, double E, double A, double I)
{
    members.push_back(Member(i, j, E, A, I));
    return members.size() - 1;
}


void Frame::buildMember(size_t index)
{
    Member &member = members[index];

    size_t i = member.iNode;
    size_t j = member.jNode;
    // printf("%u %u \n", i, j);

    Node &INode = nodes[i];
    Node &JNode = nodes[j];

    double dx, dy, ex, ey;
    dx = INode.pos[0] - JNode.pos[0];
    dy = INode.pos[1] - JNode.pos[1];

    member.l = sqrt(dx * dx + dy * dy);
    ex = dx / member.l;
    ey = dy / member.l;

    member.TransformationMatrix.ZERO();

    member.TransformationMatrix[0][0] = ex;
    member.TransformationMatrix[0][1] = -ey;
    member.TransformationMatrix[1][0] = ey;
    member.TransformationMatrix[1][1] = ex;
    member.TransformationMatrix[2][2] = 1;

    member.K11[0][0] = member.E * member.A / member.l;
    member.K11[1][1] = 12 * member.E * member.I / member.l / member.l / member.l;
    member.K11[1][2] = 6 * member.E * member.I / member.l / member.l;
    member.K11[2][1] = 6 * member.E * member.I / member.l / member.l;
    member.K11[2][2] = 4 * member.E * member.I / member.l;

    member.K12[0][0] = -member.E * member.A / member.l;
    member.K12[1][1] = -12 * member.E * member.I / member.l / member.l / member.l;
    member.K12[1][2] = 6 * member.E * member.I / member.l / member.l;
    member.K12[2][1] = -6 * member.E * member.I / member.l / member.l;
    member.K12[2][2] = 2 * member.E * member.I / member.l;

    member.K21[0][0] = -member.E * member.A / member.l;
    member.K21[1][1] = -12 * member.E * member.I / member.l / member.l / member.l;
    member.K21[1][2] = -6 * member.E * member.I / member.l / member.l;
    member.K21[2][1] = 6 * member.E * member.I / member.l / member.l;
    member.K21[2][2] = 2 * member.E * member.I / member.l;

    member.K22[0][0] = member.E * member.A / member.l;
    member.K22[1][1] = 12 * member.E * member.I / member.l / member.l / member.l;
    member.K22[1][2] = -6 * member.E * member.I / member.l / member.l;
    member.K22[2][1] = -6 * member.E * member.I / member.l / member.l;
    member.K22[2][2] = 4 * member.E * member.I / member.l;

    
#ifdef STRUCTUR_DEBUG
    printf(member.TransformationMatrix);
    printf(member.K11);
    printf(member.K12);
    printf(member.K21);
    printf(member.K22);
#endif
}

void Frame::solve()
{

    structure system(nodes.size() * DIMMINSION);

    for (size_t i = 0; i < members.size(); i++)
    {
        Member &m = members[i];
        buildMember(i);

        StructureAddStiffMat(&system, m.iNode * DIMMINSION, m.jNode * DIMMINSION, m.getGlopalK11(), m.getGlopalK12(), m.getGlopalK21(), m.getGlopalK22());
    }

    for (size_t i = 0; i < nodes.size(); i++)
    {
        for (size_t j = 0; j < DIMMINSION; j++)
        {
            if (nodes[i].support[j] == 1)
                StructureAddDis(&system,i * DIMMINSION + j, nodes[i].dis[j]);

            StructureAddForce(&system,i * DIMMINSION + j, nodes[i].nodalForce[j]);
        }
    }

    ::solve(&system);
    for (size_t i = 0; i < nodes.size(); i++)
    {
        for (size_t j = 0; j < DIMMINSION; j++)
        {

            nodes[i].dis[j] = StructureGetDis(&system,i * DIMMINSION + j);
        }
    }
    Ten y1(DIMMINSION, 1), y2(DIMMINSION, 1);
    Ten p1(DIMMINSION, 1), p2(DIMMINSION, 1);

    for (size_t i = 0; i < members.size(); i++)
    {
        Member &m = members[i];
        memcpy(y1.getData(), nodes[m.iNode].dis, DIMMINSION * sizeof(m.JendForces[0]));
        memcpy(y2.getData(), nodes[m.jNode].dis, DIMMINSION * sizeof(m.IendForces[0]));

        p1 = m.K11 * m.TransformationMatrix * y1 + m.K12 * m.TransformationMatrix * y2;
        p2 = (m.K21 * m.TransformationMatrix * y1 + m.K22 * m.TransformationMatrix * y2) * (-1);

        memcpy(m.IendForces, p1.getData(), DIMMINSION * sizeof(m.JendForces[0]));
        memcpy(m.JendForces, p2.getData(), DIMMINSION * sizeof(m.IendForces[0]));
    }
}
void Frame::print()
{
    printf("-------------Nodes----------------------\n");

    for (auto &n : nodes)
    {
        n.print();
    }
    printf("-------------Members----------------------\n");
    for (auto &m : members)
    {
        m.print();
    }
}
