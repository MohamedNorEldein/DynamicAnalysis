#include "../headers/CAPI.h"

#define pTensor Tensor<float> *

CAPI void *createTensorFloat(size_t rowNum, size_t columnNum)
{
    pTensor tf = new Tensor<float>(rowNum, columnNum);
    // generate(*tf);
    return tf;
}

CAPI void *initTensorFloat(size_t rowNum, size_t columnNum, float *data)
{
    pTensor tf = new Tensor<float>(rowNum, columnNum, data);

    return tf;
}

CAPI void deleteTensorFloat(void *tf)
{
    delete ((pTensor)tf);
}

CAPI void *addTensorFloat(void *tf1, void *tf2)
{
    pTensor result = new Tensor<float>(((pTensor)tf1)->getRowsNum(), ((pTensor)tf1)->getColumnNum());
    Tensor<float>::add((pTensor)tf1, (pTensor)tf2, result);
#ifdef DEBUG_CAPI
    printf("%x, %x, %x\n", tf1, tf2, result);
#endif

    return result;
}

CAPI void *subTensorFloat(void *tf1, void *tf2)
{
    pTensor result = new Tensor<float>(((pTensor)tf1)->getRowsNum(), ((pTensor)tf1)->getColumnNum());
    Tensor<float>::sub((pTensor)tf1, (pTensor)tf2, result);
    return result;
}

CAPI void *mullTensorFloat(void *tf1, void *tf2)
{
    pTensor result = new Tensor<float>(((pTensor)tf1)->getRowsNum(), ((pTensor)tf2)->getColumnNum());
    Tensor<float>::matmull((pTensor)tf1, (pTensor)tf2, result);
    return result;
}

CAPI void *scalarmullTensorFloat(void *tf1, float a)
{
    pTensor result = new Tensor<float>(((pTensor)tf1)->getRowsNum(), ((pTensor)tf1)->getColumnNum());
    Tensor<float>::mull(a, (pTensor)tf1, result);
    return result;
}

CAPI void *transposeTensorFloat(void *tf1)
{
    pTensor result = new Tensor<float>(((pTensor)tf1)->getRowsNum(), ((pTensor)tf1)->getColumnNum());
    Tensor<float>::transpose((pTensor)tf1, result);
    return result;
}

CAPI void *inverseTensorFloat(void *tf1)
{
    pTensor result = new Tensor<float>(((pTensor)tf1)->inverse());

    return result;
}

CAPI void printTensorFloat(void *tf)
{
    printf(*(pTensor)tf);
}

CAPI void setElementTensorFloat(void *tf, size_t i, size_t j, float value)
{
    ((pTensor)tf)->setElement(i, j, value);
    // printf("%f \n", ((pTensor)tf)->getElement(i,j));
}

CAPI float getElementTensorFloat(void *tf, size_t i, size_t j)
{
    // #ifdef DEBUG_CAPI
    //  printf("%u, %u, %f\n", i,j,((pTensor)tf)->getElement(i,j));
    // #endif
    return ((pTensor)tf)->getElement(i, j);
}

CAPI void *getCorelationMatrix(void *a)
{

    Tensor<float> *b = new Tensor<float>(correlationMatrix(*(pTensor)a));
    return b;
}

CAPI size_t getRowsNum(void *a)
{
    return ((pTensor)a)->getRowsNum();
}

CAPI size_t getColumnNum(void *a)
{
    return ((pTensor)a)->getColumnNum();
}


CAPI void generateTensorFloat(void *a,float b)
{
     generate(*((pTensor)a),b);
}
