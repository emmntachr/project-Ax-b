#ifndef ALGOS_H
#define ALGOS_H

void displayColumn(const int m, const int n, const int j, float** matrix);
void giveValuesToMatrices(int m, int n, float** matrix);
void freeMatrix(const int m, float*** matrix);
void displayMatrix(int m, int n, float** matrix);
void matrix_add(int m, int n, const float** matrix1, const float** matrix2, float** matrix3);
void matrix_mult(int m, int n, const float** matrix1, const float** matrix2, float** matrix3);
void identical_matrix(int n, float** identical);
void augmentedMatrix(const int m, const int n, float** matrix, float* b, float** result);
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
float findPivot(const int m, const int n, const int j, float** matrix);
int checkZeroRow(const int m, const int n, const int k, float** matrix);
void changeRow(const int m, const int n, const int i, const int j, float** matrix);
void putZeroRowsAtEnd(const int m, const int n, float** matrix);
void gaussianElimination(const int m, const int n, float** matrix);
void rref(const int m, const int n, float** matrix, float* b, float** resultMatrix);
int matrixRank(const int m, const int n, float** matrix);
float* multipleAx(const int m, const int n, float** A, float* x);
void backSubstitution(const int m, const int r, float** A, float* x, float* vectorb);





#endif


