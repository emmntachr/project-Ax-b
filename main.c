#include <stdio.h>
#include<stdlib.h>
#include "Algos.h"



int main()
{
    printf("Welcome to the programm that solves the lineas system of the form Ax=b.\nIn order the system to be solved you must enter \nthe dimensions m,n,  the matrix A and the vector b");
    printf("\n--------------------------------------------------------------\n\n");
    // We ask the user to give the dimensions of the matrix
    int m, n;

    printf("Enter the value for m: ");
    scanf_s("%d", &m);

    printf("Enter the value for n: ");
    scanf_s("%d", &n);

    // We define the matrices we are going to use
    float** matrix = (float**)malloc(m * sizeof(float*));
    if (matrix == NULL)
    {
        printf("Memory allocation failed. Exiting...\n");
        exit(1);
    }

    for (int i = 0; i < m; i++)
    {
        matrix[i] = (float*)malloc(n * sizeof(float));
        if (matrix[i] == NULL)
        {
            printf("Memory allocation failed. Exiting...\n");
            exit(1);
        }
    }
    float* x = (float*)malloc(n * sizeof(float));
    if (x == NULL)
    {
        printf("Memory allocation failed. Exiting...\n");
        exit(1);
    }

    float* b = (float*)malloc(m * sizeof(float));
    if (b == NULL)
    {
        printf("Memory allocation failed. Exiting...\n");
        exit(1);
    }
    float** resultMatrix = (float**)malloc(m * sizeof(float*));
    if (resultMatrix == NULL)
    {
        printf("Memory allocation failed. Exiting...\n");
        exit(1);
    }
    for (int i = 0; i < m; i++)
    {
        resultMatrix[i] = (float*)malloc((n + 1) * sizeof(float));
        if (resultMatrix[i] == NULL)
        {
            printf("Memory allocation failed. Exiting...\n");
            exit(1);
        }
    }
    float* z = (float*)malloc(n * sizeof(float));

    if (z == NULL) {
        printf("Memory allocation failed. Exiting...\n");
        return 1; // Exit with an error code
    }

    // We ask user to insert values to the matrix
    printf("Enter the elements of the matrix:\n");
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("Enter element at position (%d, %d): ", i + 1, j + 1);
            scanf_s("%f", &matrix[i][j]);
        }
    }

    // Using malloc we define the matrix matrixU with dimensions mxn and the vector, vectorC with dimension mx1.
    float** matrixU;
    float* vectorC;
    matrixU = (float**)malloc(m * sizeof(float*));
    for (int i = 0; i < m; i++)
    {
        matrixU[i] = (float*)malloc(n * sizeof(float));
    }

    vectorC = (float*)malloc(m * sizeof(float));

    if (matrixU == NULL || vectorC == NULL)
    {
        printf("Memory allocation failed. Exiting...\n");
        exit(1);
    }
    // We ask the  user to insert values to the vector b
    printf("Enter the elements of vector b:\n");
    for (int i = 0; i < m; i++) {
        printf("Enter element at position %d: ", i + 1);
        scanf_s("%f", &b[i]);
    }
    /*
        We store in the matrix resultMatrix the augmented system and then we bring the system to row reduced echelon form, 
        hence the system goes from [A|b] to [U|c].
        In case the system has one row fullof zeroes beside the last column then the system has no solution and the program terminates
    */
    printf("--------------------------------------------------------------\n\n");
    printf("The matrix provided is: \n\n");
    displayMatrix(m, n, matrix);
    printf("The vector b provided is: \n\n ");
    printf("b=[");
    for (int i = 0; i < m - 1; i++)
    {
        printf("%.2f,", b[i]);
    }
    printf("%.2f]^T \n\n", b[m - 1]);
    printf("--------------------------------------------------------------\n\n");
    rref(m, n, matrix, b, resultMatrix);
    putZeroRowsAtEnd(m, n + 1, resultMatrix);
    printf("The row reduced echelon form of the augmented matrix [A|b] is: \n\n");
    displayMatrix(m, n + 1, resultMatrix);
    printf("\n\n\n");

    // Check if the system has no solution, if yes exit else continue.
    for (int i = 0; i < m; i++)
    {
        if (checkZeroRow(m, n, i, resultMatrix) && resultMatrix[i][n] != 0)
        {
            printf("The system has no solution.\n\n\n");
            exit(0);
            return;
        }
    }
    //We put values to matrixU and vectorC according to the matrix resultMatrix.
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n + 1; j++)
        {
            if (j < n)
            {
                matrixU[i][j] = resultMatrix[i][j];
            }
        }
    }
    for (i = 0; i < m; i++)
    {
        vectorC[i] = resultMatrix[i][n];
    }
    /*
        Next we find the free variables and the pivot variables. Free variables are the one who have no no 1  in a column.
        So we define a matrix y with n elements with values 1 and 0 and if a variable is free the value is 0 else the value is 1
    */

    int* y = (int*)malloc(n * sizeof(int));
    if (y == NULL)
    {
        printf("Memory allocation failed. Exiting...\n");
        exit(1);
    }
    for (j = 0; j < n; j++)// We put values to the matrix y equal to 1
    {
        y[j] = 1;
    }
    j = 0;
    for (i = 0; i < m; i++)
    {
        while (matrixU[i][j] == 0 && j < n)
        {
            j++;
        }
        if (matrixU[i][j] == 1)
        {
            y[j] = 0;
        }
        j = 0;

    }

    // The matrix nullmatrix is m x (n-r), where r is the rank of the matrix matrixU
    int r;
    r = matrixRank(m, n, matrixU);
    printf("--------------------------------------------------------------\n\n");


    float** nullmatrix = (int**)malloc(n * sizeof(float*));

    // Check if memory allocation was successful
    if (nullmatrix == NULL) 
    {
        printf("Memory allocation failed.\n");
        return 1;
    }

    // Allocate memory for each row of the matrix
    for (int i = 0; i < n; i++) {
        nullmatrix[i] = (float*)malloc((n - r) * sizeof(float));

        // Check if memory allocation was successful for the current row
        if (nullmatrix[i] == NULL) {
            printf("Memory allocation failed for row %d.\n", i);
            return 1;
        }
    }
    float* nullvector;
    nullvector = (float*)malloc(n * sizeof(float));
    // Check if memory allocation was successful
    if (nullvector == NULL)
    {
        printf("Memory allocation failed.\n");
        return 1;
    }
    float* vectorR;
    vectorR = (float*)malloc(m * sizeof(float));
    // Check if memory allocation was successful
    if (vectorR == NULL)
    {
        printf("Memory allocation failed.\n");
        return 1;
    }

    // We check if the system is square and if it is full rank. If yes then it has a unique solution.
    if (m == n && m == r) 
    {
        printf("The system has a unique solution:x=[");
        backSubstitution(m, n, matrixU, nullvector, vectorC);
        for (i = 0; i < n-1; i++)
        {
            printf("%.2f,", nullvector[i]);
        }
        printf("%.2f]^T \n\n\n ", nullvector[n-1]);
        exit(0);
    }

    int count;
    count = 0;
    int a;
    for (a = 0; a < n; a++)
    {
        if (y[a] == 1)
        {
            //We put the value -1 to the nullvector matrix
            for (i = 0; i < n; i++)
            {
                nullvector[i] = -1;
            }
            nullvector[a] = 1;
            //We add values 0,1 to the nullvector
            for (i = 0; i < n; i++)
            {
                if (y[i] == 1)
                {
                    for (j = 0; j < n; j++)
                    {
                        if (y[j] == 1 && j != a)
                            nullvector[j] = 0;
                    }
                }
            }
            // We add values to the vectorR with the column which is the free var for the a of the present step
            if (y[a] == 1)
            {
                for (i = 0; i < m; i++)
                {
                    vectorR[i] = matrixU[i][a];
                }
            }
            int findPivot = 0;
            // We send the solution vector vectorR to the nullvector in the elements with values -1
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if (matrixU[i][j] == 1)
                    {
                        nullvector[j] = vectorR[i] * (-1);
                        break;
                    }
                }
            }
            // We store each nullvector to the nullmatrix as a column.
            for (j = 0; j < n; j++)
            {
                nullmatrix[j][count] = nullvector[j];
            }
            count += 1;
        }

    }

    printf("\nSolving the system Ax=0, the solutions are the column vectors of the below matrix: \n\n");
    displayMatrix(n, n - r, nullmatrix);
    // Now we find a special solution x for the system Ax=b, when each free variable take the value 0.

    for (a = 0; a < n; a++)
    {
        if (y[a] == 1)
        {
            //We put the value -1 to the nullvector matrix
            for (i = 0; i < n; i++)
            {
                nullvector[i] = -1;
            }
            nullvector[a] = 0;
            //We add values 0 to the nullvector
            for (i = 0; i < n; i++)
            {
                if (y[i] == 1)
                {
                    for (j = 0; j < n; j++)
                    {
                        if (y[j] == 1)
                            nullvector[j] = 0;
                    }
                }
            }
            // We add values to the vectorR with the column which is the free var for the a of the present step


            for (i = 0; i < m; i++)
            {
                vectorR[i] = vectorC[i];
            }

            int findPivot = 0;
            // We send the solution vector vectorR to the nullvector in the elements with values -1
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if (matrixU[i][j] == 1)
                    {
                        nullvector[j] = vectorR[i];
                        break;
                    }
                }
            }
        }

    }
    printf("--------------------------------------------------------------\n\n");
    printf("Solving the system Ax=b, for the ginen A,b the particular solution placing zeroes to the free variables are:\n");
    printf("[");
    for (i = 0; i < n - 1; i++)
    {
        printf("%.2f,", nullvector[i]);
    }
    printf("%.2f]\n\n", nullvector[n - 1]);

    // Now we print the solution of the system Ax=b using the nullMatrix we found solving Ax=0 and the nullvector solving Ax=b.
    printf("--------------------------------------------------------------\n\n");

    printf("\n\n");
    printf(" The general solution of the system is:\n x=[");
    for (i = 0; i < n - 1; i++)
    {
        printf("%.2f,", nullvector[i]);
    }
    printf("%.2f]^T +", nullvector[n-1]);
    for (i = 0; i < n - r; i++)
    {
        printf("k%d", i);
        displayColumn(n, n - r, i, nullmatrix);
        printf("^T ");

        if (i < n - r - 1)
        {
            printf("+");
        }
    }
    printf("\n\n\n");
    printf("------");

    return 0;

}




