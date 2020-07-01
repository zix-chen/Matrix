#ifndef PROJECT_MATRIX_H
#define PROJECT_MATRIX_H

#include <iostream>
#include <complex>
#include "vector"

using namespace std;

/*class Complex{
private:
    double imag;
    double real;
public:
    Complex();
    Complex(double, double);
    Complex operator+(const Complex&);
    Complex operator-(const Complex&);
    Complex operator*(const Complex&);
    friend Complex operator*(double, const Complex&);
    Complex operator~();
    int operator==(const Complex&);
    int operator!=(const Complex&);
    friend std::ostream& operator<<(std::ostream& out, const Complex& u);
    friend std::istream& operator >> (std::istream& in, Complex& u);
    ~Complex();
};*/
//warning: functions without note "//tested" still needs testing (might have bugs), take care
class Matrix {
private:
    int row;
    int col;
    complex<double> *value;
    int sparse;


public:
    bool cross_multiple_Judge(const Matrix &);

    bool dot_multiple_Judge(const Matrix &);


    bool inverse_Judge();//是否方阵

    Matrix inverse();//复数矩阵求逆

    Matrix inverseOfMatrixReal();//求矩阵实数部分的逆矩阵
    Matrix inverseOfMatrixImag();//求矩阵虚数部分的的逆矩阵
    Matrix getFromImagPartOfMatrix();//将虚数部分放到实数部分构成一个新矩阵
    Matrix getFromRealPartOfMatrix();//将虚数部分删除构成一个新矩阵
    double *getReal();

    double *getImag();

    Matrix conv2(Matrix &y); //卷积


    Matrix();

    Matrix(int, int, complex<double> *);

    Matrix(int, int, double *);

    Matrix(Matrix &);

    Matrix(int, int, int, complex<double> *);

    Matrix(int row, int col, double *initial_real, double *initial_imag);

    Matrix operator+(const Matrix);

    Matrix operator-(const Matrix);

    virtual Matrix operator*(Matrix &matrix);

    Matrix operator*(complex<double>);

    friend Matrix operator*(complex<double>, Matrix);

    Matrix operator/(double a);

    Matrix operator/(complex<double>);


    int getRow();

    int getCol();

    void SetMatrixEle(complex<double> ele);

    complex<double> getValue(int row, int col);

    complex<double> *getValue();


    Matrix elementMulti(Matrix);

    Matrix transposition();

    Matrix conjugation();

    Matrix reshape(int Newrow, int Newcol);

    complex<double> det();

    complex<double> eleAverage();

    complex<double> eleMax();

    complex<double> eleSum();

    complex<double> *slice(int startRow, int startCol, int endRow, int endCol);

    complex<double> *slice(int startRow, int startCol, int endRow, int endCol, int step);

    void print();

    Matrix dec(Matrix);

    complex<double> trace();

    void expand();

    void deflation();

    void QRdecomposition(Matrix A, double *Q, double *R);

    complex<double> *REAL_E_value(Matrix A);

    complex<double> *CPLX_E_value(Matrix A);

    Matrix REAL_E_vector(Matrix A);

    Matrix CPLX_E_vector(Matrix A);


};

template<typename T>
class others : public Matrix {
private:
    int row;
    int col;
    T *value;
public:
    others();

    others(int, int, T *);

    others operator+(others<T> &);
};


#endif //PROJECT_MATRIX_H