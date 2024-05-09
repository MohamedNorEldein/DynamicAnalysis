#include "../headers/Frame.h"





int main(){
    printf("start\n");
    Frame fr;
    size_t n0 =  fr.addNode(0,0);
    size_t n1 =fr.addNode(0,1);
    size_t n2 =fr.addNode(0,2);

    size_t n3 =fr.addNode(0,3);
    size_t n4 =fr.addNode(1,3);
    size_t n5 =fr.addNode(2,3);
    size_t n6 =fr.addNode(3,3);
    
    size_t n7 =fr.addNode(3,2);
    size_t n8 =fr.addNode(3,1);
    size_t n9 =fr.addNode(3,0);

    printf("nodes are defined\n");

    fr.addMember(n0,n1,1,1000,1);
    fr.addMember(n1,n2,1,1000,1);
    fr.addMember(n2,n3,1,1000,1);
    fr.addMember(n3,n4,1,1000,1);
    fr.addMember(n4,n5,1,1000,1);
    fr.addMember(n5,n6,1,1000,1);
    fr.addMember(n6,n7,1,1000,1);
    fr.addMember(n7,n8,1,1000,1);
    fr.addMember(n8,n9,1,1000,1);
    
    printf("members are defined\n");

    fr.addSupport(n0,0,0);
    fr.addSupport(n0,1,0);
    fr.addSupport(n9,0,0);
   // fr.addSupport(n9,1,0);

    printf("supports are defined\n");

    fr.addnodeForce(n3,1,-1);
    fr.addnodeForce(n4,1,-1);
    fr.addnodeForce(n5,1,-1);
    fr.addnodeForce(n6,1,-1);

    printf("nodal forces are defined\n");

    fr.solve();

    printf("system is solved \n");


    fr.print();
    return 0;
}