#pragma once

#include <iostream>
#include "vec.h"


template <typename Type>
Tensor<Type> correlationMatrix(Tensor<Type> &a)
{
    Tensor<Type> b = a.transpose();

    return ( b*a).inverse()*b;
}

