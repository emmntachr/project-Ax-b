#include <stdio.h>
#include <stdlib.h>
#include "Algos.h"

void displayColumn(const int m, const int n, const int j, float** matrix)
{
    printf("( ");
    for (int i = 0; i < m-1; i++)
    {
        printf("%.2f, ", matrix[i][j]);
    }
    printf(" %.2f )",matrix[m-1][j]);
}

void giveValuesToMatrices(int m, int n, float** matrix)
{
    if (matrix == NULL)
    {
        printf("Error: Invalid matrix!\n");
        return 1;
    }

    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("Enter a value for the matrix at position (%d, %d): ", i, j);
            scanf_s("%f", &matrix[i][j]);
        }
    }
}

void freeMatrix(const int m, float*** matrix) {
    int i;
    for (i = 0; i < m; i++) {
        free((*matrix)[i]);
    }
    free(*matrix);
    *matrix = NULL;
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------

// It displays a 2D matrix
void displayMatrix(int m, int n, float** matrix)
{
    int i, j;
    for (i = 0; i < m; i++)
    {
        printf("| ");
        for (j = 0; j < n; j++)
        {
            printf("%.2f\t", matrix[i][j]);
        }
        printf("|\n");
    }
}

// This function add to matrices matrix1 and matrix2
void matrix_add(int m, int n, const float** matrix1, const float** matrix2, float** matrix3)
{
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            matrix3[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
}

//This function multiplies tho matrices
void matrix_mult(int m, int n, const float** matrix1, const float** matrix2, float** matrix3)
{
    int i, j, k;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            matrix3[i][j] = 0.0;
            for (k = 0; k < m; k++)
            {

                matrix3[i][j] = matrix3[i][j] + matrix1[i][k] * matrix2[k][j];
            }
        }
    }
}

/*
    This function given a value n it creates an identical matrix
    whose diagonal elements are 1 and the other elements are 0
*/
void identical_matrix(int n, float** identical)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j)
            {
                identical[i][j] = 1;
            }
            else
            {
                identical[i][j] = 0;
            }
        }
    }
}

// The code given a column j search for the pivot of a matrix with dimensions m x n
float findPivot(const int m, const int n, const int j, float** matrix)
{
    int i;
    float pivot=0.0;
    for (i = j; i < m; i++) 
    {
        if (matrix[i][j] != 0)
        {
            pivot = matrix[i][j];
            return pivot;
        }
    }
    return pivot; // Return a default value if no non zero pivot is found
}

//This function changes the rows i and j in a matrix
void changeRow(const int m, const int n, const int i, const int j, float** matrix)
{
    //matrix is copied to copymatrix
    float temp;
    int k;
    for (k = 0; k < n; k++)
    {
        temp = matrix[i][k];
        matrix[i][k] = matrix[j][k];
        matrix[j][k] = temp;
    }
}

// The Gauss algorithm takes as input a matrix and returns an upper triangular matrix
void gaussianElimination(const int m, const int n, float** matrix)
{
    int i, j, k;
    float pivot;

    for (j = 0; j < n - 1; j++) 
    {
        pivot = findPivot(m, n, j, matrix); // Find the pivot of a column j

        if (pivot == 0.0) 
            continue;

        // It switches the row
        for (i = j; i < m; i++)
        {
            if (matrix[i][j] == pivot)
            {
                changeRow(m, n, i, j, matrix);
                break;
            }
        }

        // It does linear transformations
        for (k = j + 1; k < m; k++)
        {
            float factor = matrix[k][j] / pivot;
            for (i = j; i < n; i++)// We start from the j and not the 0 value because we fixed the previous elements in previous steps.
                matrix[k][i] = matrix[k][i]- factor * matrix[j][i];
        }
    }
    // Now the matrix is upper triangular
}







// This function creates an augmented matrix give the matrix matrix and the vector b.
void augmentedMatrix(const int m, const int n, float** matrix, float* b, float** result)
{
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n + 1; j++)
        {
            if (j < n)
            {
                result[i][j] = matrix[i][j];
            }
            else
            {
                result[i][j] = b[i];

            }
        }
    }
}
// This function checks if a matrix has a row full of zeroes
int checkZeroRow(const int m, const int n, const int k, float** matrix)
{

    int zeroRow = 1;// i  has value 1 (true) if a column i is full of 0 else the zeroRow becomes 0 (false)
    int j;
    for (j = 0; j < n; j++)
    {
        if (matrix[k][j] != 0)
        {
            zeroRow = 0;
            break;
        }
    }

    return zeroRow;
}
//This function switches a row full of 0 in the bottom rows.
void putZeroRowsAtEnd(const int m, const int n, float** matrix)
{
    // Initialize two pointers, one pointing to the first row and the other to the last row
    int firstNonZeroRow = 0;
    int lastRow = m - 1;

    // Iterate through the rows of the matrix
    while (firstNonZeroRow < lastRow)
    {
        // Check if the current row is a zero row
        if (checkZeroRow(m, n, firstNonZeroRow, matrix))
        {
            // If it's a zero row, search for the next non-zero row from the bottom
            while (checkZeroRow(m, n, lastRow, matrix) && lastRow > firstNonZeroRow)
            {
                lastRow--;
            }

            // Swap the zero row with the last non-zero row
            if (lastRow > firstNonZeroRow)
            {
                changeRow(m, n, firstNonZeroRow, lastRow, matrix);
                lastRow--;
            }
        }

        // Move to the next row
        firstNonZeroRow++;
    }
}

//This function produces the row echelon form of a given system
void rref(const int m, const int n, float** matrix, float* b, float** resultMatrix)
{
    augmentedMatrix(m, n, matrix, b, resultMatrix);
    gaussianElimination(m, n + 1, resultMatrix);

    int i, j;

    // Do row operations for row-reduced echelon form (rref)
    for (i = m - 1; i >= 0; i--)
    {
        // Find pivot column
        j = 0;
        while (j < n && resultMatrix[i][j] == 0)
        {
            j++;
        }

        if (j == n) // Skip if row is all zeros
        {
            continue;
        }

        // Do pivot element 1
        float pivot = resultMatrix[i][j];
        for (int k = j; k <= n; k++)
        {
            resultMatrix[i][k] /= pivot;
        }

        // Eliminate non zero elements above the pivot
        for (int a = i - 1; a >= 0; a--)
        {
            float factor = resultMatrix[a][j];
            for (int k = j; k <= n; k++)
            {
                resultMatrix[a][k] -= factor * resultMatrix[i][k];
            }
        }
    }
    putZeroRowsAtEnd(m, n + 1, resultMatrix);

}

//This function multiplies a matrix A with a vector x as a dot product 
float* multipleAx(const int m, const int n, float** A, float* x) 
{
    float* result = (float*)malloc(m * sizeof(float));
    if (result == NULL) 
    {
        printf("Memory allocation failed. Exiting...\n");
        return NULL; // Return NULL to indicate an error
    }

    // Perform matrix-vector multiplication
    int i, j;
    for (i = 0; i < m; i++) 
    {
        result[i] = 0.0;
        for (j = 0; j < n; j++) 
        {
            result[i] += A[i][j] * x[j];
        }
    }

    return result;
}
// This function returns the rank of a given matrix 
int matrixRank(const int m, const int n, float** matrix)
{
    int i, j,k;
    gaussianElimination(m, n, matrix);

    // Perform row operations for row-reduced echelon form (rref)
    for (i = m - 1; i >= 0; i--)
    {
        // Find pivot column
        j = 0;
        while (j < n && matrix[i][j] == 0)
        {
            j++;
        }

        if (j == n) // Skip if row is all zeros
        {
            continue;
        }

        // Make pivot element 1
        float pivot = matrix[i][j];
        for (int k = j; k <= n; k++)
        {
            matrix[i][k] /= pivot;
        }

        // Eliminate non-zero elements above the pivot
        for (int a = i - 1; a >= 0; a--)
        {
            float factor = matrix[a][j];
            for (int k = j; k <= n; k++)
            {
                matrix[a][k] -= factor * matrix[i][k];
            }
        }
    }
    // Now the matrix is in rref form
    int count;
    count = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (matrix[i][j] == 1)
            {
                count++;
                break;
            }
        }
    }
    return count;

}
// This function performs backSubtitution to a matrix and find a solution fillling the matrix x
void backSubstitution(const int m, const int r, float** A, float* x, float* vectorb) 
{
    int i, j;
    float sum;

    for (i = r - 1; i >= 0; i--) {
        // Check if the row is not full of zeros
        int nonZeroRow = 0;
        for (j = 0; j < m; j++) {
            if (A[i][j] != 0) {
                nonZeroRow = 1;
                break;
            }
        }

        if (nonZeroRow) {
            // Do back substitution for non-zero row
            sum = 0.0;
            for (j = i + 1; j < m; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (vectorb[i] - sum) / A[i][i];
        }
        else {
            // Row is full of zeros, set x[i] to 0
            x[i] = 0.0;
        }
    }
}



















