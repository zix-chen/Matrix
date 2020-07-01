#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "Matrix.h"

#define N 4
#define DEBUG 1

void matrix_inverse_LU(double a[][N],double a_inverse[N][N]);
int main() {
    double a_real[N][N],a_i[N][N];
    double a_imag[N][N];
    a_real[0][0] = 4;a_real[0][1] = 3;a_real[0][2] = 4;a_real[0][3] = 5;
    a_real[1][0] = 1;a_real[1][1] = 8;a_real[1][2] = 2;a_real[1][3] = 9;
    a_real[2][0] = 4;a_real[2][1] = 5;a_real[2][2] = 1;a_real[2][3] = 7;
    a_real[3][0] = 6;a_real[3][1] = 7;a_real[3][2] = 1;a_real[3][3] = 2;

    a_imag[0][0] = 2;a_imag[0][1] = 1;a_imag[0][2] = 3;a_imag[0][3] = 5;
    a_imag[1][0] = 7;a_imag[1][1] = 2;a_imag[1][2] = 2;a_imag[1][3] = 3;
    a_imag[2][0] = 4;a_imag[2][1] = 6;a_imag[2][2] = 7;a_imag[2][3] = 2;
    a_imag[3][0] = 1;a_imag[3][1] = 8;a_imag[3][2] = 4;a_imag[3][3] = 5;

    double arr[N * N];
    double arr2[N * N];
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            arr[i * N + j] = a_real[i][j];
            arr2[i * N + j] = a_imag[i][j];
        }
    }
    Matrix matrix(N, N, arr, arr2);
//    matrix.inverseOfMatrixImag();
    matrix.inverse();
    double arr3[] = {1, 2, 3, 1, 2, 3, 1, 2, 3};
    double arr4[] = {2, 2, 2, 1, 1, 1, 3, 3, 3};
    auto one = Matrix(3, 3, arr3);
    auto two = Matrix(3, 3, arr4);
    Matrix m = one.conv2(two);
    m.print();
    return 0;

}


