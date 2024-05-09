#pragma once
#include <math.h>
#include <stdlib.h>

typedef unsigned char uint8;

template <typename Type, uint8 DIMENSION>
struct vec
{
    /* data */
    Type data[DIMENSION]{0};

    /* methods */
    Type &operator[](uint8 i)
    {
        if (i < DIMENSION)
        {
            return data[i];
        }
        throw 0;
    }

    Type dot(const vec<Type, DIMENSION> &b)
    {
        Type c = 0;
        for (uint8 i = 0; i < DIMENSION; i++)
        {
            c += data[i] * b.data[i];
        }
        return c;
    }

    vec<Type, DIMENSION> &operator+=(const vec<Type, DIMENSION> &b)
    {

        for (uint8 i = 0; i < DIMENSION; i++)
        {
            data[i] += b.data[i];
        }
        return *this;
    }

    vec<Type, DIMENSION> &operator-=(const vec<Type, DIMENSION> &b)
    {

        for (uint8 i = 0; i < DIMENSION; i++)
        {
            data[i] -= b.data[i];
        }
        return *this;
    }

    vec<Type, DIMENSION> &operator*=(Type b)
    {

        for (uint8 i = 0; i < DIMENSION; i++)
        {
            data[i] = b * data[i];
        }
        return *this;
    }

    vec<Type, DIMENSION> operator+(const vec<Type, DIMENSION> &b)
    {
        vec<Type, DIMENSION> c;
        for (uint8 i = 0; i < DIMENSION; i++)
        {
            c.data[i] = data[i] + b.data[i];
        }
        return c;
    }

    vec<Type, DIMENSION> operator-(const vec<Type, DIMENSION> &b)
    {
        vec<Type, DIMENSION> c;
        for (uint8 i = 0; i < DIMENSION; i++)
        {
            c.data[i] = data[i] - b.data[i];
        }
        return c;
    }

    vec<Type, DIMENSION> operator*(Type b)
    {
        vec<Type, DIMENSION> c;
        for (uint8 i = 0; i < DIMENSION; i++)
        {
            c.data[i] = data[i] * b;
        }
        return c;
    }
};

/*************************************************/
/*          Matrix<Type,DIMENSION>          */

template <typename Type, uint8 DIMENSION>
struct Matrix
{
    /* data */
    Type data[DIMENSION][DIMENSION];

    /* methods */

    Type *operator[](uint8 i)
    {
        return data[i];
    }

    Matrix<Type, DIMENSION> &operator+=(Matrix<Type, DIMENSION> &b)
    {
        for (uint8 i = 0; i < DIMENSION; i++)
        {
            for (uint8 j = 0; j < DIMENSION; j++)
            {
                data[i][j] += b.data[i][j];
            }
        }
        return *this;
    }

    Matrix<Type, DIMENSION> &operator-=(Matrix<Type, DIMENSION> &b)
    {
        for (uint8 i = 0; i < DIMENSION; i++)
        {
            for (uint8 j = 0; j < DIMENSION; j++)
            {
                data[i][j] -= b.data[i][j];
            }
        }
        return *this;
    }

    Matrix<Type, DIMENSION> &operator*=(Type b)
    {
        for (uint8 i = 0; i < DIMENSION; i++)
        {
            for (uint8 j = 0; j < DIMENSION; j++)
            {
                data[i][j] *= b;
            }
        }
        return *this;
    }

    Matrix<Type, DIMENSION> operator-(const Matrix<Type, DIMENSION> &b)
    {
        Matrix<Type, DIMENSION> c{0};
        for (uint8 i = 0; i < DIMENSION; i++)
        {
            for (uint8 j = 0; j < DIMENSION; j++)
            {
                c.data[i][j] = data[i][j] - b.data[i][j];
            }
        }
        return c;
    }

    Matrix<Type, DIMENSION> operator+(const Matrix<Type, DIMENSION> &b)
    {
        Matrix<Type, DIMENSION> c{0};
        for (uint8 i = 0; i < DIMENSION; i++)
        {
            for (uint8 j = 0; j < DIMENSION; j++)
            {
                c.data[i][j] = data[i][j] + b.data[i][j];
            }
        }
        return c;
    }

    Matrix<Type, DIMENSION> operator*(Type b)
    {
        Matrix<Type, DIMENSION> c{0};
        for (uint8 i = 0; i < DIMENSION; i++)
        {
            for (uint8 j = 0; j < DIMENSION; j++)
            {
                c.data[i][j] = data[i][j] * b;
            }
        }
        return c;
    }

    Matrix<Type, DIMENSION> operator*(const Matrix<Type, DIMENSION> &b)
    {
        Matrix<Type, DIMENSION> c{0};

        for (uint8 i = 0; i < DIMENSION; i++)
        {
            for (uint8 j = 0; j < DIMENSION; j++)
            {
                c.data[i][j] = 0;
                for (uint8 k = 0; k < DIMENSION; k++)
                {
                    c.data[i][j] += data[i][k] * b.data[k][j];
                }
            }
        }

        return c;
    }

    vec<Type, DIMENSION> operator*(const vec<Type, DIMENSION> &b)
    {
        vec<Type, DIMENSION> c{0};
        for (uint8 i = 0; i < DIMENSION; i++)
        {
            c.data[i] = 0;
            for (uint8 j = 0; j < DIMENSION; j++)
            {
                c.data[i] += data[i][j] * b.data[j];
            }
        }
        return c;
    }

    vec<Type, DIMENSION> solve(const vec<Type, DIMENSION> &b)
    {
        Type v = 0;
        int k = 0;
        vec<Type, DIMENSION> c;
        while (k < DIMENSION * DIMENSION)
        {

            for (uint8 i = 0; i < DIMENSION; i++)
            {
                v = 0;
                for (uint8 j = 0; j < DIMENSION; j++)
                {
                    v -= data[i][j] * c.data[j];
                }
                v += data[i][i] * c.data[i];
                c.data[i] = (b.data[i] + v) / data[i][i];
            }
            k++;
        }
        return c;
    }

    vec<Type, DIMENSION> solve(const vec<Type, DIMENSION> &b, bool *restrains, Type *restrainValue)
    {
        Type v = 0;
        int k = 0;
        vec<Type, DIMENSION> c;
        while (k < DIMENSION * DIMENSION)
        {

            for (uint8 i = 0; i < DIMENSION; i++)
            {
                if (restrains[i] == true)
                    c.data[i] = restrainValue[i];
                else
                {
                    v = 0;
                    for (uint8 j = 0; j < DIMENSION; j++)
                    {

                        v -= data[i][j] * c.data[j];
                    }
                    v += data[i][i] * c.data[i];
                    c.data[i] = (b.data[i] + v) / data[i][i];
                }
            }
            k++;
        }
        return c;
    }

    Matrix<Type, DIMENSION> Transpose()
    {
        Matrix<Type, DIMENSION> c;
        for (uint8 i = 0; i < DIMENSION; i++)
        {
            for (uint8 j = 0; j < DIMENSION; j++)
            {
                c.data[i][j] = data[j][i];
            }
        }
        return c;
    }

private:
    static void swapRows(Type row1[], Type row2[], int size)
    {
        for (int i = 0; i < size; i++)
        {
            Type temp = row1[i];
            row1[i] = row2[i];
            row2[i] = temp;
        }
    }

    static void scaleRow(Type row[], Type scalar, int size)
    {
        for (int i = 0; i < size; i++)
        {
            row[i] *= scalar;
        }
    }

    static void addRows(Type row1[], Type row2[], Type scalar, int size)
    {
        for (int i = 0; i < size; i++)
        {
            row1[i] += (row2[i] * scalar);
        }
    }

public:
    static int inverseMatrix(Matrix<Type, DIMENSION> *mat, Matrix<Type, DIMENSION> *result)
    {
        Matrix<Type, DIMENSION> identity;
        int size = DIMENSION;

        // Initialize the identity matrix
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                if (i == j)
                {
                    identity.data[i][j] = 1.0;
                }
                else
                {
                    identity.data[i][j] = 0.0;
                }
            }
        }

        // Apply Gauss-Jordan elimination
        for (int i = 0; i < size; i++)
        {
            if (mat->data[i][i] == 0)
            {
                int found = 0;
                for (int j = i + 1; j < size; j++)
                {
                    if (mat->data[j][i] != 0)
                    {
                        swapRows(mat->data[i], mat->data[j], size);
                        swapRows(identity.data[i], identity.data[j], size);
                        found = 1;
                        break;
                    }
                }

                if (!found)
                {
                    printf("Matrix is not invertible.\n");
                    return 0;
                }
            }

            Type pivot = mat->data[i][i];
            scaleRow(mat->data[i], 1.0 / pivot, size);
            scaleRow(identity.data[i], 1.0 / pivot, size);

            for (int j = 0; j < size; j++)
            {
                if (j != i)
                {
                    Type factor = -mat->data[j][i];
                    addRows(mat->data[j], mat->data[i], factor, size);
                    addRows(identity.data[j], identity.data[i], factor, size);
                }
            }
        }

        // Copy the inverse matrix to the result
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                result->data[i][j] = identity.data[i][j];
            }
        }

        return 1;
    }

    Matrix<Type, DIMENSION> inverse()
    {
        Matrix<Type, DIMENSION> result;
        inverseMatrix(this, &result);
        return result;
    }
};

/****************tensor****************/

#ifdef DEBUG_VECTORS
static size_t counter;
#endif

template <typename Type>
class Tensor
{
private:
    Type *data;
    size_t DIMENSION1; // no of rows
    size_t DIMENSION2; // no of coulmns = no of elements in a row

private:
    void copy(Type *otherTensor, size_t startI = 0, Type startV = 0)
    {

#ifdef DEBUG_VECTORS
        printf("****************************copy\n");
#endif

        if (otherTensor == 0)
        {
            for (size_t i = 0; i < DIMENSION1 * DIMENSION2; i++)
            {
                data[i] = startV;
            }
        }
        else
        {
            for (size_t i = 0; i < startI; i++)
            {
                data[i] = startV;
            }

            for (size_t i = startI; i < DIMENSION1 * DIMENSION2; i++)
            {
                data[i] = otherTensor[i];
            }
        }
    }

public:
    static void add(Tensor *tensor1, Tensor *tensor2, Tensor *result)
    {
        size_t DIMENSION1 = tensor1->getRowsNum();
        size_t DIMENSION2 = tensor1->getColumnNum();

        for (size_t i = 0; i < DIMENSION1; i++)
        {
            for (size_t j = 0; j < DIMENSION2; j++)
            {
                result->setElement(i, j, tensor1->getElement(i, j) + tensor2->getElement(i, j));
            }
        }
    }

    static void sub(Tensor *tensor1, Tensor *tensor2, Tensor *result)
    {
        size_t DIMENSION1 = tensor1->getRowsNum();
        size_t DIMENSION2 = tensor1->getColumnNum();

        for (size_t i = 0; i < DIMENSION1; i++)
        {
            for (size_t j = 0; j < DIMENSION2; j++)
            {
                result->setElement(i, j, tensor1->getElement(i, j) - tensor2->getElement(i, j));
            }
        }
    }

    static void _elementMul(Tensor *tensor1, Tensor *tensor2, Tensor *result)
    {
        size_t DIMENSION1 = tensor1->getRowsNum();
        size_t DIMENSION2 = tensor1->getColumnNum();

        for (size_t i = 0; i < DIMENSION1; i++)
        {
            for (size_t j = 0; j < DIMENSION2; j++)
            {
                result->setElement(i, j, tensor1->getElement(i, j) * tensor2->getElement(i, j));
            }
        }
    }

    static void mull(Type scaler, Tensor *tensor, Tensor *result)
    {
        size_t DIMENSION1 = tensor->getRowsNum();
        size_t DIMENSION2 = tensor->getColumnNum();

        for (size_t i = 0; i < DIMENSION1; i++)
        {
            for (size_t j = 0; j < DIMENSION2; j++)
            {
                result->setElement(i, j, tensor->getElement(i, j) * scaler);
            }
        }
    }

    static void transpose(Tensor *tensor1, Tensor *result)
    {
        size_t DIMENSION1 = result->getRowsNum();
        size_t DIMENSION2 = result->getColumnNum();

        for (size_t i = 0; i < DIMENSION1; i++)
        {
            for (size_t j = 0; j < DIMENSION2; j++)
            {

                result->setElement(i, j, tensor1->getElement(j, i));
            }
        }
    }

    static void matmull(Tensor *tensor1, Tensor *tensor2, Tensor *result)
    {
        size_t DIMENSION1 = tensor1->getRowsNum();
        size_t DIMENSION2 = tensor1->getColumnNum();
        size_t DIMENSION3 = tensor2->getColumnNum();

#ifdef DEBUG_VECTORS
        printf("(%u x %u) * (%u x %u)\n", tensor1->getRowsNum(), tensor1->getColumnNum(), tensor2->getRowsNum(), tensor2->getColumnNum());
#endif

        if (tensor1->getColumnNum() != tensor2->getRowsNum())
        {
            printf("matmull dimminsion error\n");
            throw 0;
        }

        for (size_t i = 0; i < DIMENSION1; i++)
        {
            for (size_t j = 0; j < DIMENSION3; j++)
            {
                Type sum = 0.0;
                for (size_t k = 0; k < DIMENSION2; k++)
                {
                    sum += tensor1->getElement(i, k) * tensor2->getElement(k, j);
                }
                result->setElement(i, j, sum);
            }
        }
    }

    static Type regresionMull_inner(Tensor<Type> *matrix, Tensor<Type> *vector, size_t index, size_t index2)
    {
        Type a = matrix->getElement(index, 0);

        size_t count = vector->getRowsNum();

        for (size_t i = 0; i < count; i++)
        {
            a += matrix->getElement(index, i + 1) * vector->getElement(i, index2);
        }

        return a;
    }

    static int stat_regresionMull(Tensor<Type> *matrix, Tensor<Type> *vector, Tensor<Type> *result)
    {
        if (matrix->getColumnNum() == vector->getRowsNum() + 1)
        {
            size_t count1 = matrix->getRowsNum(), count2 = vector->getColumnNum();

            for (size_t i = 0; i < count1; i++)
            {
                for (size_t j = 0; j < count2; j++)
                {
                    result->setElement(i, 0, regresionMull_inner(matrix, vector, i, j));
                }
            }
            return 0;
        }

        return 1;
    }

private:
    void swapRows(Type *matrix, int dimension, int row1, int row2)
    {
        for (int j = 0; j < dimension; j++)
        {
            Type temp = matrix[row1 * dimension + j];
            matrix[row1 * dimension + j] = matrix[row2 * dimension + j];
            matrix[row2 * dimension + j] = temp;
        }
    }

    void scaleRow(Type *matrix, int dimension, int row, Type scalar)
    {
        for (int j = 0; j < dimension; j++)
        {
            matrix[row * dimension + j] *= scalar;
        }
    }

    void addScaledRow(Type *matrix, int dimension, int row1, int row2, Type scalar)
    {
        for (int j = 0; j < dimension; j++)
        {
            matrix[row1 * dimension + j] += scalar * matrix[row2 * dimension + j];
        }
    }

    void invertMatrix(const Type *matrix, int dimension, Type *result)
    {
        // Create an augmented matrix [matrix | identity matrix]
        int augmentedDimension = 2 * dimension;
        Type *augmentedMatrix = new Type[augmentedDimension * dimension];

        for (int i = 0; i < dimension; i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                augmentedMatrix[i * augmentedDimension + j] = matrix[i * dimension + j];
                augmentedMatrix[i * augmentedDimension + j + dimension] = (i == j) ? 1.0f : 0.0f;
            }
        }

        // Apply Gauss-Jordan elimination
        for (int i = 0; i < dimension; i++)
        {
            // Find pivot row
            int pivotRow = -1;
            for (int j = i; j < dimension; j++)
            {
                if (augmentedMatrix[j * augmentedDimension + i] != 0.0f)
                {
                    pivotRow = j;
                    break;
                }
            }

            if (pivotRow == -1)
            {
                // Matrix is singular, inverse does not exist
                delete[] augmentedMatrix;
                return;
            }

            // Swap pivot row with current row
            swapRows(augmentedMatrix, augmentedDimension, i, pivotRow);

            // Scale pivot row
            Type pivot = augmentedMatrix[i * augmentedDimension + i];
            scaleRow(augmentedMatrix, augmentedDimension, i, 1.0f / pivot);

            // Eliminate other rows
            for (int j = 0; j < dimension; j++)
            {
                if (j != i)
                {
                    Type factor = augmentedMatrix[j * augmentedDimension + i];
                    addScaledRow(augmentedMatrix, augmentedDimension, j, i, -factor);
                }
            }
#ifdef DEBUG_VECTORS
            printf("inverted\n");
#endif
        }

        // Extract the inverse matrix from the augmented matrix
        for (int i = 0; i < dimension; i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                result[i * dimension + j] = augmentedMatrix[i * augmentedDimension + j + dimension];
            }
        }

        delete[] augmentedMatrix;
    }

public:
    Tensor(size_t DIMENSION1, size_t DIMENSION2) : DIMENSION1(DIMENSION1), DIMENSION2(DIMENSION2)
    {

        data = (Type *)malloc(DIMENSION1 * DIMENSION2 * sizeof(Type));

#ifdef DEBUG_VECTORS
        printf("created\n");
        counter++;
#endif
    }
    Tensor(void) : DIMENSION1(0), DIMENSION2(0)
    {
        data = 0;
#ifdef DEBUG_VECTORS
        printf("created\n");
        counter++;
#endif
    }

    Tensor(size_t DIMENSION1, size_t DIMENSION2, float *fdata, size_t startI = 0, Type startV = 0) : DIMENSION1(DIMENSION1), DIMENSION2(DIMENSION2)
    {

        data = (Type *)malloc(DIMENSION1 * DIMENSION2 * sizeof(Type));
        copy(fdata, startI, startV);

#ifdef DEBUG_VECTORS
        printf("created\n");
        counter++;
#endif
    }

    Tensor(Tensor &otherTensor) : DIMENSION1(otherTensor.DIMENSION1), DIMENSION2(otherTensor.DIMENSION2)
    {

        data = (Type *)malloc(DIMENSION1 * DIMENSION2 * sizeof(Type));
        copy(otherTensor.data);

#ifdef DEBUG_VECTORS
        printf("created\n");
        counter++;
#endif
    }

    Tensor(Tensor &&otherTensor) : DIMENSION1(otherTensor.DIMENSION1), DIMENSION2(otherTensor.DIMENSION2)
    {
        data = otherTensor.data;
        otherTensor.data = 0;

#ifdef DEBUG_VECTORS
        printf("created\n");
        counter++;
#endif
    }

    ~Tensor()
    {
        free(data);

#ifdef DEBUG_VECTORS
        printf("Tensor deleted\n");
        counter--;
#endif
    }

public:
    Type *operator[](size_t i) const
    {
        return data + DIMENSION2 * i;
    }

    Type& getElement(size_t i, size_t j)
    {
        return data[i * DIMENSION2 + j];
    }

    void setElement(size_t i, size_t j, Type value)
    {
        data[i * DIMENSION2 + j] = value;
    }

    size_t getRowsNum() const
    {
        return DIMENSION1;
    }

    size_t getColumnNum() const
    {
        return DIMENSION2;
    }

public:
    Tensor operator+(const Tensor &other) const
    {
        if (DIMENSION1 != other.getRowsNum())
            throw 0;
        if (DIMENSION2 != other.getColumnNum())
            throw 0;

        Tensor result(DIMENSION1, DIMENSION2);
        add((Tensor *)this, (Tensor *)(&other), &result);
        return result;
    }

    const Tensor &operator+=(const Tensor &other) const
    {
        if (DIMENSION1 != other.getRowsNum())
            throw 0;
        if (DIMENSION2 != other.getColumnNum())
            throw 0;

        add((Tensor *)this, (Tensor *)(&other), (Tensor *)this);
        return *this;
    }

    void operator=(const Tensor &other)
    {
        if (DIMENSION1 != other.getRowsNum())
            throw 0;
        if (DIMENSION2 != other.getColumnNum())
            throw 0;

        for (size_t i = 0; i < DIMENSION1 * DIMENSION2; i++)
        {
            data[i] = other.data[i];
        }
    }

    Tensor operator-(const Tensor &other) const
    {
        if (DIMENSION1 != other.getRowsNum())
            throw 0;
        if (DIMENSION2 != other.getColumnNum())
            throw 0;

        Tensor result(DIMENSION1, DIMENSION2);
        sub((Tensor *)this, (Tensor *)(&other), &result);
        return result;
    }

    Tensor operator*(Type scaler) const
    {
        Tensor result(DIMENSION1, DIMENSION2);
        mull(scaler, (Tensor *)this, &result);
        return result;
    }

    void operator*=(Type scaler) const
    {

        mull(scaler, (Tensor *)this, (Tensor *)this);
    }

    Tensor operator*(const Tensor &other) const
    {

        if (DIMENSION2 != other.getRowsNum())
            throw 0;
        Tensor result(DIMENSION1, other.DIMENSION2);
        matmull((Tensor *)this, (Tensor *)(&other), &result);
        return result;
    }

    Tensor elementMul(const Tensor &other) const
    {
        if (DIMENSION1 != other.DIMENSION1)
            throw 0;

        if (DIMENSION2 != other.DIMENSION2)
            throw 0;

        Tensor result(DIMENSION1, DIMENSION2);
        _elementMul((Tensor *)this, (Tensor *)(&other), &result);
        return result;
    }

    void set(const Tensor &other)
    {
        delete this->data;

        DIMENSION1 = other.DIMENSION1;
        DIMENSION2 = other.DIMENSION2;

        this->data = (Type *)calloc(DIMENSION1 * DIMENSION2, sizeof(Type));

        //        print(other.data,DIMENSION1 * DIMENSION2);
        copy(other.data);
        // printf("%x %x\n",data, other.data);
    }

    void set(Tensor &&other)
    {

        DIMENSION1 = other.DIMENSION1;
        DIMENSION2 = other.DIMENSION2;

        data = other.data;
        other.data = 0;
    }

public:
    Tensor inverse()
    {
        size_t DIMENSION1 = this->getRowsNum();
        size_t DIMENSION2 = this->getColumnNum();

        if (DIMENSION2 != DIMENSION1)
            return Tensor(0, 0);

        Tensor result(DIMENSION1, DIMENSION2);

        invertMatrix(((Tensor *)this)->data, ((Tensor *)this)->DIMENSION1, result.data);

        return result;
    }

    Tensor transpose() const
    {
        Tensor result(DIMENSION2, DIMENSION1);
        transpose((Tensor *)this, &result);
        return result;
    }

private:
    void copyRow(Type *srcdata, Type *desdata, size_t oldStart_index, size_t newStart_index)
    {

        for (size_t i = 0; i < DIMENSION2; i++)
        {
            desdata[newStart_index * DIMENSION2 + i] = srcdata[oldStart_index * DIMENSION2 + i];
        }
    }

    void copyColumn(Type *srcdata, Type *desdata, size_t oldStart_index, size_t newStart_index)
    {
        for (size_t i = 0; i < DIMENSION1; i++)
        {
            desdata[newStart_index + (DIMENSION2 + 1) * i] = srcdata[oldStart_index + DIMENSION2 * i];
        }
    }

public:
    void addSubSpace(Tensor &otherTensor, size_t startI, size_t startJ)
    {

#ifdef DEBUG_VECTORS
        printf("****************************copy\n");
        printf(otherTensor);
        printf("%d , %d\n",otherTensor.DIMENSION1, otherTensor.DIMENSION2);

#endif

        for (size_t i = 0; i <  otherTensor.DIMENSION1; i++)
        {
            for (size_t j = 0; j <  otherTensor.DIMENSION2; j++)
            {
                this->getElement(startI+ i,startJ+ j) += otherTensor[i][j];
            }
        }
    }

    void setRow(size_t index, Type *dValues)
    {
        for (size_t i = 0; i < DIMENSION2; i++)
        {
            data[index * DIMENSION2 + i] = dValues[i];
        }
    }

    void setColumn(size_t index, Type *dValues)
    {
        for (size_t i = 0; i < DIMENSION1; i++)
        {
            data[index + (DIMENSION2)*i] = dValues[i];
        }
    }

    void insertColumn(size_t index, Type constantValue)
    {
        if (index > DIMENSION2)
            return;

        Type *newData = new Type[DIMENSION1 * (DIMENSION2 + 1)];

        for (size_t i = 0; i < index; i++)
        {
            copyColumn(data, newData, i, i);
        }

        for (size_t i = 0; i < DIMENSION1; i++)
        {
            newData[index + (DIMENSION2 + 1) * i] = constantValue;
        }

        for (size_t i = index + 1; i < DIMENSION2 + 1; i++)
        {
            copyColumn(data, newData, i - 1, i);
        }

        DIMENSION2++;
        delete data;
        data = newData;
    }

    void insertRow(size_t index, Type constantValue)
    {
        if (index > DIMENSION1)
            return;

        Type *newData = new Type[DIMENSION2 * (DIMENSION1 + 1)];

        for (size_t i = 0; i < index; i++)
        {
            copyRow(data, newData, i, i);
        }

        for (size_t i = 0; i < DIMENSION2; i++)
        {
            newData[index * DIMENSION2 + i] = constantValue;
        }

        for (size_t i = index + 1; i < DIMENSION1; i++)
        {
            copyRow(data, newData, i - 1, i);
        }

        DIMENSION1++;
        delete data;
        data = newData;
    }

    void copyData(Type *fdata, size_t startIndex = 0, Type startvalue = 0)
    {
        copy(fdata, startIndex, startvalue);
    }

    void setData(Type *fdata, size_t n, size_t m)
    {
        DIMENSION1 = n;
        DIMENSION2 = m;
        data = fdata;
    }

    void ZERO()
    {
        copy(0);
    }

    void ONES()
    {
        copy(0, 0, 1);
    }

    Type *getData()
    {
        return data;
    }

    Type norm() const
    {
        Type n = 0;
        for (size_t i = 0; i < DIMENSION1 * DIMENSION2; i++)
        {
            n += data[i] * data[i];
        }
        return sqrt(n);
    }

public:
    Tensor<Type> applyFunction(Type (*func)(Type)) const
    {
        Tensor result(DIMENSION1, DIMENSION2);

        for (size_t i = 0; i < DIMENSION1 * DIMENSION2; i++)
        {
            result.data[i] = func(data[i]);
        }
        return result;
    }

    static Tensor<Type> regresionMull(const Tensor<Type> &matrix, const Tensor<Type> &vector)
    {
        Tensor<Type> result(matrix.getRowsNum(), vector.getColumnNum());
        stat_regresionMull((Tensor<Type> *)&matrix, (Tensor<Type> *)&vector, &result);
        return result;
    }
};

/***************DEBUG_VECTORS*****************/

#include <iostream>
template <typename Type, uint8 DIMENSION>
void printf(Matrix<Type, DIMENSION> &mat)
{
    printf("{\n");
    for (short i = 0; i < DIMENSION; i++)
    {
        for (short j = 0; j < DIMENSION; j++)
        {
            printf(" %f ", mat.data[i][j]);
        }
        printf("\n");
    }
    printf("}\n");
}

template <typename Type, uint8 DIMENSION>
void printf(vec<Type, DIMENSION> &vec)
{
    std::cout << ("<");
    for (short i = 0; i < DIMENSION; i++)
    {
        printf(" %f ", vec[i]);
    }
    printf(">\n");
}

template <typename Type>
void printf(const Tensor<Type> &tensor)
{

    size_t DIMENSION1 = tensor.getRowsNum();
    size_t DIMENSION2 = tensor.getColumnNum();
    printf("{\n");
    // Display the result
    for (size_t i = 0; i < DIMENSION1; i++)
    {
        std::cout << " " << (tensor[i][0]);
        for (size_t j = 1; j < DIMENSION2; j++)
        {
            std::cout << ", " << tensor[i][j];
        }
        std::cout << ",\n";
    }
    printf("}\n");
}

template <typename Type, uint8 DIMENSION>
void generate(vec<Type, DIMENSION> &vec)
{
    for (short i = 0; i < DIMENSION; i++)
    {
        vec->data[i] = Type((rand() % 1000)) / 10;
    }
}

template <typename Type, uint8 DIMENSION>
void generate(Matrix<Type, DIMENSION> &mat)
{
    for (int i = 0; i < DIMENSION; i++)
    {
        for (int j = 0; j < DIMENSION; j++)
        {
            if (i == j)
                mat->data[i][j] = 10;
            else
                mat->data[i][j] = Type((rand() % 1000)) / 100;
        }
    }
}

template <typename Type>
void generate(Tensor<Type> &tensor, float a = 100)
{
    size_t DIMENSION1 = tensor.getRowsNum();
    size_t DIMENSION2 = tensor.getColumnNum();
    // Set tensor elements
    for (size_t i = 0; i < DIMENSION1; i++)
    {
        for (size_t j = 0; j < DIMENSION2; j++)
        {
            tensor.setElement(i, j, Type(rand() % 100) / a);
        }
    }
}

template <typename Type>
Tensor<Type> abs(Tensor<Type> &tensor)
{
    return tensor.applyFunction(abs);
}

template <typename Type>
void generate(Type *data, size_t num)
{
    for (size_t i = 0; i < num; i++)
    {
        data[i] = Type(rand() % 100);
    }
}

template <typename Type>
void print(Type *data, size_t num)
{
    for (size_t i = 0; i < num; i++)
    {
        printf("%f ", data[i]);
    }
    printf("\n");
}

template <typename Type>
Type derivative(Type x, Type (*func)(Type), Type dx = 0.00001)
{
    return (func(x + dx) - func(x)) / dx;
}

template <typename Type>
Type solve(Type x, Type (*func)(Type), Type ys = 0, Type dx = 0.00001, size_t count = 10)
{
    float df = 0, f = 0;
    for (size_t i = 0; i < count; i++)
    {
        df = derivative(x, func, dx);
        f = ys - func(x);
        if (df == 0 || f == 0)
            return x;

        x = x - f / df;
    }
    return x;
}
