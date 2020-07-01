#include "Matrix.h"
#include <complex>
#include <memory.h>

#define DEBUG 1
using namespace std;

bool Matrix::cross_multiple_Judge(const Matrix &matrix) {
    if (this->col == matrix.row)
        return true;
    else
        return false;
}

bool Matrix::dot_multiple_Judge(const Matrix &matrix) {
    if (this->col == matrix.col && this->row == matrix.row)
        return true;
    else
        return false;
}

//-------------------------------------Constructor-----------------------------------//
Matrix::Matrix() {//tested
    row = 1;
    col = 1;
    complex<double> initial[1] = {0};
    value = initial;
    sparse = 0;
}

Matrix::Matrix(int row, int col, complex<double> *initial) {//tested
    this->row = row;
    this->col = col;
    this->value = initial;
    sparse = 0;
}

Matrix::Matrix(int row, int col, double *initial) {//tested
    complex<double> *value1 = new complex<double>[row * col];
    for (int i = 0; i < col * row; ++i) {
        //cout << initial[i];
        value1[i] = complex<double>(initial[i], 0);
        //cout << value1[i];
    }
    //cout<<endl;
    this->value = value1;
    this->row = row;
    this->col = col;
    sparse = 0;
}

Matrix::Matrix(int row, int col, double *initial_real, double *initial_imag) {//tested
    complex<double> *value1 = new complex<double>[row * col];
    for (int i = 0; i < col * row; ++i) {
        //cout << initial[i];
        value1[i] = complex<double>(initial_real[i], initial_imag[i]);
        //cout << value1[i];
    }
    //cout<<endl;
    this->value = value1;
    this->row = row;
    this->col = col;
    sparse = 0;
}

Matrix::Matrix(Matrix &matrix) {//tested
    this->row = matrix.row;
    this->col = matrix.col;
    this->value = matrix.value;
    this->sparse = matrix.sparse;
}

Matrix::Matrix(int sparse, int row, int col, complex<double> *value) {
    this->value = value;
    this->row = row;
    this->col = col;
    this->sparse = 1;
}


//---------------------------------------------matrix operator-------------------------------------------//
Matrix Matrix::operator+(const Matrix matrix) {
    this->expand();
    if (!this->dot_multiple_Judge(matrix)) {
        cerr << "Invalid matrix shape!" << endl;
        Matrix defaultMat(row, col, value);
        return defaultMat;
    }
    complex<double> *Newvalue = new complex<double>[row * col];
    for (int i = 0; i < this->row * this->col; ++i) {
        Newvalue[i] = this->value[i] + matrix.value[i];
    }
    this->deflation();
    Matrix result(row, col, Newvalue);
    return result;
}

Matrix Matrix::operator-(const Matrix matrix) {
    this->expand();
    if (!this->dot_multiple_Judge(matrix)) {
        cerr << "Invalid matrix shape!" << endl;
        Matrix defaultMat(row, col, value);
        return defaultMat;
    }
    complex<double> *Newvalue = new complex<double>[row * col];
    for (int i = 0; i < this->row * this->col; ++i) {
        Newvalue[i] = this->value[i] - matrix.value[i];
    }
    this->deflation();
    Matrix result(row, col, Newvalue);
    return result;
}

Matrix Matrix::operator*(Matrix &matrix) {//tested
    this->expand();
    matrix.expand();
    if (!this->cross_multiple_Judge(matrix)) {
        cerr << "Invalid matrix shape!" << endl;
        Matrix defaultMat(row, col, value);
        return defaultMat;
    }

    complex<double> *result = new complex<double>[row * matrix.col];
    int num = 0;
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < matrix.col; ++j) {
            complex<double> a(0, 0);
            for (int p = 0; p < matrix.row; ++p) {
                //cout << (col*i)+p << ", " << (matrix.col*p)+j << endl;
                //cout << this->value[(col*i)+p] <<", "<< matrix.value[(matrix.col*p)+j] << endl;
                a = a + this->value[(col * i) + p] * matrix.value[(matrix.col * p) + j];
                //a = a + this->reverse()[i][p] * matrix.reverse()[p][i];
                //cout << a << endl;
            }
            result[num] = a;
            num += 1;
        }
    }
    this->deflation();
    matrix.deflation();
    Matrix matrix_result(row, matrix.col, result);
    return matrix_result;
}

Matrix Matrix::operator*(complex<double> a) {
    this->expand();
    complex<double> *Newvalue = new complex<double>[row * col];
    for (int i = 0; i < row * col; ++i) {
        Newvalue[i] = this->value[i] * a;
    }
    this->deflation();
    Matrix matrix(row, col, Newvalue);
    return matrix;
}

Matrix operator*(complex<double> a, Matrix matrix) {
    matrix.expand();
    complex<double> *Newvalue = new complex<double>[matrix.row * matrix.col];
    for (int i = 0; i < matrix.row * matrix.col; ++i) {
        Newvalue[i] = matrix.value[i] * a;
    }
    matrix.deflation();
    Matrix result(matrix.row, matrix.col, Newvalue);
    return result;
}

Matrix Matrix::operator/(double a) {
    this->expand();
    complex<double> *Newvalue = new complex<double>[row * col];
    for (int i = 0; i < row * col; ++i) {
        Newvalue[i] = this->value[i] / a;
    }
    this->deflation();
    Matrix matrix(row, col, Newvalue);
    return matrix;
}

Matrix Matrix::operator/(complex<double> a) {
    this->expand();
    complex<double> *Newvalue = new complex<double>[row * col];
    for (int i = 0; i < row * col; ++i) {
        Newvalue[i] = this->value[i] / a;
    }
    this->deflation();
    Matrix matrix(row, col, Newvalue);
    return matrix;
}

//------------------------------------Matrix Value getting function--------------------------//
int Matrix::getRow() {
    return row;
}

int Matrix::getCol() {
    return col;
}

void Matrix::SetMatrixEle(complex<double> ele) {
    for (int i = 0; i < row * col; ++i) {
        value[i] = ele;
    }
}

complex<double> Matrix::getValue(int Vrow, int Vcol) {
    return value[(Vrow - 1) * col + Vcol - 1];
}

complex<double> *Matrix::getValue() {
    return value;
}

//------------------------------------------------matrix functions---------------------------//
Matrix Matrix::elementMulti(Matrix matrix) {
    this->expand();
    matrix.expand();
    if (!this->dot_multiple_Judge(matrix)) {
        cerr << "Invalid matrix shape!" << endl;
        Matrix defaultMat(row, col, value);
        return defaultMat;
    }

    complex<double> *Newvalue = new complex<double>[row * col];
    for (int i = 0; i < row * col; ++i) {
        Newvalue[i] = this->value[i] * matrix.value[i];
    }
    this->deflation();
    matrix.deflation();
    Matrix mat(row, col, Newvalue);
    return mat;
}

Matrix Matrix::transposition() {//tested
    this->expand();
    int newrow = col;
    int newcol = row;
    complex<double> *newmatrix = new complex<double>[col * row];
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            newmatrix[row * j + i] = value[col * i + j];
        }
    }
    this->deflation();
    Matrix result(newrow, newcol, newmatrix);
    return result;
}

Matrix Matrix::conjugation() {//队友写的
    this->expand();
    complex<double> *newmatrix = new complex<double>[col * row];
    for (int i = 0; i < row * col; ++i) {
        newmatrix[i] = conj(value[i]);
    }
    this->deflation();
    Matrix result(row, col, newmatrix);
    return result;
}

complex<double> Matrix::eleAverage() {//tested
    this->expand();
    double real = 0;
    double imag = 0;
    for (int i = 0; i < row * col; ++i) {
        real = real + value[i].real();
        imag = imag + value[i].imag();
    }
    this->deflation();
    double avgReal = real / (row * col);
    double avgImag = imag / (row * col);
    return complex<double>(avgReal, avgImag);
}

complex<double> Matrix::eleMax() {//tested
    this->expand();
    double absolute = 0;
    complex<double> a = value[0];
    for (int i = 0; i < row * col; ++i) {
        if (absolute < abs(value[i])) {
            absolute = abs(value[i]);
            a = value[i];
        }
    }
    this->deflation();
    return a;
}

complex<double> Matrix::eleSum() {//tested
    this->expand();
    complex<double> a = 0;
    for (int i = 0; i < row * col; ++i) {
        a = a + value[i];
    }
    this->deflation();
    return a;
}

Matrix Matrix::reshape(int Newrow, int Newcol) {//tested
    this->expand();
    Matrix result(row, col, value);
    if (row * col == Newcol * Newrow) {
        result.row = Newrow;
        result.col = Newcol;
    }
    this->deflation();
    return result;
}

complex<double> *Matrix::slice(int startRow, int startCol, int endRow, int endCol) {//tested
    //参数中的行数和列数对应了真实的行数和列数，即从1开始
    this->expand();
    int startindex = startCol + (col * (startRow - 1)) - 1;
    int endindex = endCol + (col * (endRow - 1)) - 1;
    int arraylength = endindex - startindex + 1;
    if (startindex >= endindex) {
        cerr << "Invalid Index number" << endl;
        return value;
    }
    if (endindex >= row * col) {
        cerr << "Index out of range" << endl;
        return value;
    }
    //cout<<arraylength<<" "<<startindex<<" "<<endindex;
    complex<double> *slicevalue = new complex<double>[arraylength];
    for (int i = 0; i < arraylength; i++) {
        slicevalue[i] = value[startindex + i];
        //cout<<slicevalue[i];
    }
    this->deflation();
    return slicevalue;
}

complex<double> *Matrix::slice(int startRow, int startCol, int endRow, int endCol, int step) {//tested
    //参数中的行数和列数对应了真实的行数和列数，即从1开始
    this->expand();
    int startindex = startCol + (col * (startRow - 1)) - 1;
    int endindex = endCol + (col * (endRow - 1)) - 1;
    int arraylength = endindex - startindex + 1;
    if (startindex >= endindex) {
        cerr << "Invalid Index number" << endl;
        return value;
    }
    if (endindex >= row * col) {
        cerr << "Index out of range" << endl;
        return value;
    }
    //cout<<arraylength<<" "<<startindex<<" "<<endindex;
    complex<double> *slicevalue = new complex<double>[arraylength];

    int i = 0;
    while (i <= arraylength / step) {
        slicevalue[i] = value[startindex + i];
        //cout<<slicevalue[i];
        i = i + step;
    }
    this->deflation();
    return slicevalue;
}

complex<double> Matrix::det() {//队友写的
    if (row != col) {
        cerr << "Matrix without full rank, det = 0!" << endl;
        return 0;
    }
    this->expand();
    complex<double> result = 0;
    if (this->row == 2) {
        result = value[0] * value[3] - value[1] * value[2];
    } else {
        for (int i = 0; i < this->row; ++i) {
            complex<double> *submatrix_value[(row - 1) * (col - 1)];
            int num = 0;
            for (int j = 0; j < row * col; ++j) {
                if ((j - i) % row == 0)
                    continue;
                submatrix_value[num] = &this->value[j];
                num += 1;
            }
            Matrix submatrix(row - 1, col - 1, *submatrix_value);
            result = result + this->value[i] * submatrix.det() * pow(-1, i);
        }
    }
    this->deflation();
    return result;
}

complex<double> Matrix::trace() {//tested
    if (row != col) {
        cerr << "Invalid shape of matrix!" << endl;
        return complex<double>(0, 0);
    }
    this->expand();
    complex<double> result = (0, 0);
    for (int i = 0; i < row; i++) {
        result = result + getValue(i + 1, i + 1);
    }
    return result;
}

void Matrix::print() {
    if (sparse == 1) {
        int count = 0;
        for (int i = 0;; ++i) {
            cout << value[2 * i].real() << "*" << value[2 * i + 1] << ", ";
            count += value[2 * i].real();
            if (count == row * col)
                break;
        }
        cout << "end";
        //this->expand();
        //this->print();
    } else {
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                cout << (this->value[i * col + j]).real() << "+" << (this->value[i * col + j]).imag() << "I ";
            }
            cout << endl;
        }
    }
}

void Matrix::expand() {//队友写的
    if (sparse == 1) {
        complex<double> change[row * col];
        //value = new complex<double>[row*col];
        int num = 0;
        int pointer_move = 0;
        while (num <= row * col - 1) {
            //cout << pointer_move << endl;
            for (int i = 0; i < this->value[pointer_move].real(); ++i) {
                change[num] = value[pointer_move + 1];
                num += 1;
            }
            pointer_move += 2;
        }
        for (int j = 0; j < row * col; ++j) {
            this->value[j] = change[j];
        }
        sparse = 2;//sparse=2代表矩阵可还原为稀疏矩阵存储
    }
}

void Matrix::deflation() {//队友写的
    if (sparse == 2) {
        complex<double> *value = new complex<double>[2 * row * col];
        int num = 0;
        int pointer_move = 0;
        while (true) {
            int count = 1;
            while (this->value[num + 1].real() == this->value[num].real() &&
                   this->value[num + 1].imag() == this->value[num].imag()) {
                count += 1;
                num += 1;
            }
            const double a = (double) count;
            value[pointer_move] = complex<double>(a, 0);
            //value[pointer_move].imag(0);
            value[pointer_move + 1] = this->value[num];
            pointer_move += 2;
            if (num == row * col - 1)
                break;
            num += 1;
        }
        for (int i = 0; i <= pointer_move; ++i) {
            this->value[i] = value[i];
        }
        //this->value = value;
        sparse = 1;
    }
}

void Matrix::QRdecomposition(Matrix A, double *Q, double *R)//队友写的
{
    //进行单次QR分解的模块
    unsigned int i, j, k, m;
    const unsigned N = A.row;
    double temp;
    double *a = new double[N];
    double *b = new double[N];
    for (j = 0; j < N; ++j) {
        for (i = 0; i < N; ++i) {
            a[i] = b[i] = A.value[i * A.col + j].real();
        }
        for (k = 0; k < j; ++k) {
            R[k * N + j] = 0;
            for (m = 0; m < N; ++m) {
                R[k * N + j] += a[m] * Q[m * N + k];//
            }
            for (m = 0; m < N; ++m) {
                b[m] -= R[k * N + j] * Q[m * N + k];
            }
        }
        double norm = 0;
        for (i = 0; i < N; ++i) {
            norm += b[i] * b[i];
        }
        temp = (double) sqrt(norm);
        R[j * N + j] = temp;
        for (i = 0; i < N; ++i) {
            Q[i * N + j] = b[i] / temp;
        }
    }
    delete a, b;
}

complex<double> *Matrix::REAL_E_value(Matrix A) {//队友写的
    int C_num = A.col;
    complex<double> *eigenValue = new complex<double>[C_num];

    if (C_num != A.row) {
        cerr << "Invalid Matrix shape, cannot get E-Value!" << endl;
        return eigenValue;
    }

    const int iteratorUPPER_BOUND = 50;
    double *testQ = new double[C_num * C_num];
    double *testR = new double[C_num * C_num];
    Matrix Q(C_num, C_num, testQ);
    Matrix R(C_num, C_num, testR);
    Matrix iteratingMAT = A;
    for (int k = 0; k < iteratorUPPER_BOUND; ++k) {
        A.QRdecomposition(iteratingMAT, testQ, testR);
        Q = Matrix(3, 3, testQ);
        R = Matrix(3, 3, testR);
        iteratingMAT = R * Q;
    }
    //iteratingMAT.print();
    //R.print();
    //Q.print();

    for (int k = 1; k <= C_num; ++k) {
        eigenValue[k - 1] = iteratingMAT.getValue(k, k);
        //cout<<eigenValue[k-1]<<endl;
    }
    ///*
    int index = 0;
    complex<double> max;
    for (int i = 0; i < C_num; i++) {
        max = eigenValue[i];
        for (int j = i; j < C_num; j++) {
            if (eigenValue[j].real() > max.real()) {
                max = eigenValue[j];
                index = j;
            }//find the biggest after it
        }
        eigenValue[index] = eigenValue[i];
        eigenValue[i] = max;//then exchange
    }
    //*/
    //Qresult = testQ;
    return eigenValue;
}

Matrix Matrix::REAL_E_vector(Matrix A) {//队友写的
    int C_num = A.col;
    complex<double> *eigenValue = new complex<double>[C_num];
    double *eigenVector = new double[C_num * C_num];
    if (C_num != A.row) {
        cerr << "Invalid Matrix shape, cannot get E-Value!" << endl;
        return A;
    }
    eigenValue = A.REAL_E_value(A);

    double Nth_evalue = 0;
    double *iterator = new double[C_num * C_num];
    double tempSum, tempELE;
    //for(int k = 0;k < C_num*C_num;++k)
    //{
    //    cout<<eigenVector[k]<<endl;
    //}//测试模块
    for (int count = 0; count < C_num; ++count) {
        //拉取一个eValue，得到其对应evector的系数矩阵iterator(这里其实命名错误了，跟迭代没啥关系的)
        Nth_evalue = eigenValue[count].real();
        for (int i = 0; i < C_num * C_num; i++) {
            iterator[i] = (A.value[i]).real();
        }
        for (int i = 0; i < C_num; ++i) {
            iterator[i * C_num + i] -= Nth_evalue;
        }

        //把iterator化成阶梯阵后导出每一个evalue对应的evector到result矩阵中
        for (int i = 0; i < C_num - 1; ++i) {
            tempELE = iterator[i * C_num + i];
            for (int j = i; j < C_num; ++j) {
                iterator[i * C_num + j] /= tempELE;
            }

            for (int j = i + 1; j < C_num; ++j) {
                tempELE = iterator[j * C_num + i];
                for (int q = i; q < C_num; ++q) {
                    iterator[j * C_num + q] -= tempELE * iterator[i * C_num + q];
                }
            }
        }

        //
        tempSum = eigenVector[(C_num - 1) * C_num + count] = 1;
        //eigenVector[(C_num - 1) * C_num + count] = 1;
        //
        for (int m = C_num - 2; m >= 0; --m) {
            double sum = 0;
            for (int j = m + 1; j < C_num; ++j) {
                sum += iterator[m * C_num + j] * eigenVector[j * C_num + count];
            }
            sum = -sum / iterator[m * C_num + m];
            tempSum += sum * sum;
            eigenVector[m * C_num + count] = sum;
        }

        tempSum = sqrt(tempSum);
        for (int i = 0; i < C_num; ++i) {
            eigenVector[i * C_num + count] /= tempSum;
        }
        Matrix result(C_num, C_num, eigenVector);
        //result.print();
    }
    Matrix result(C_num, C_num, eigenVector);
    return result;
}


complex<double> *Matrix::CPLX_E_value(Matrix A) {//队友写的
    int C_num = A.col;
    complex<double> *eigenValue = new complex<double>[C_num];
    complex<double> *imag = new complex<double>[C_num];
    double *A_IMAG = new double[C_num * C_num];
    if (C_num != A.row) {
        cerr << "Invalid Matrix shape, cannot get E-Value!" << endl;
        return eigenValue;
    }
    for (int i = 0; i < C_num * C_num; i++) {
        A_IMAG[i] = A.value[i].imag();
    }
    Matrix A_CPLX(C_num, C_num, A_IMAG);
    eigenValue = A.REAL_E_value(A);
    imag = A_CPLX.REAL_E_value(A_CPLX);
    for (int i = 0; i < C_num; i++) {
        eigenValue[i].imag(imag->real());
    }
    return eigenValue;
    /*
    double *ABBA = new double[4 * C_num * C_num];
    for (int i = 0; i < C_num; ++i){
        for (int j = 0; j < C_num; ++j) {
            ABBA[i*2*C_num + j] = A.value[i*C_num + j].real();
            ABBA[i*2*C_num + j+C_num] = -A.value[i*C_num + j].imag();
            ABBA[(i+C_num)*2*C_num + j] = -ABBA[i*2*C_num + j+C_num];
            ABBA[(i+C_num)*2*C_num + j+C_num] = ABBA[i*2*C_num + j];
        }
    }
     */
}

Matrix Matrix::CPLX_E_vector(Matrix A) {//队友写的
    int C_num = A.col;
    double *imag_eigenVector = new double[C_num * C_num];
    double *imag_A = new double[C_num * C_num];
    Matrix eigenVector, eigenVector_CPLX;
    if (C_num != A.row) {
        cerr << "Invalid Matrix shape, cannot get E-Value!" << endl;
        return A;
    }
    for (int i = 0; i < C_num * C_num; i++) {
        imag_A[i] = A.value[i].imag();
    }
    Matrix A_IMAG(C_num, C_num, imag_A);
    eigenVector = A.REAL_E_vector(A);
    eigenVector_CPLX = A_IMAG.REAL_E_vector(A);
    for (int i = 0; i < C_num * C_num; i++) {
        eigenVector.value[i].imag(eigenVector_CPLX.value[i].real());
    }
    return eigenVector;
}

template<typename T>
others<T>::others() {//队友写的
    row = 0;
    col = 0;
    value = NULL;
}

template<typename T>
others<T>::others(int row, int col, T *value) {//队友写的
    this->row = row;
    this->col = col;
    this->value = value;
}

template<typename T>
others<T> others<T>::operator+(others<T> &a) {//队友写的
    T *result = new T[row * col];
    for (int i = 0; i < row * col; ++i) {
        result[i] = this->value[i] + a.value[i];
    }
    others<T> resultt(row, col, result);
    return resultt;
}

Matrix Matrix::inverseOfMatrixReal() {
    int n = row;
    double l[n][n], u[n][n];
    double l_inverse[n][n], u_inverse[n][n];
    double_t a_inverse[n * n];

    memset(l, 0, sizeof(l));
    memset(u, 0, sizeof(u));
    memset(l_inverse, 0, sizeof(l_inverse));
    memset(u_inverse, 0, sizeof(u_inverse));
    memset(a_inverse, 0, sizeof(a_inverse));

    double s, t;
    for (int i = 0; i < n; ++i) {
        l[i][i] = 1;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            s = 0;
            for (int k = 0; k < i; ++k) {
                s += l[i][k] * u[k][j];
            }
            u[i][j] = getValue(i + 1, j + 1).real() - s;
        }
        for (int j = 0; j < n; ++j) {
            s = 0;
            for (int k = 0; k < i; ++k) {
                s += l[j][k] * u[k][i];
            }
            l[j][i] = (getValue(j + 1, i + 1).real() - s) / u[i][i];
        }
    }


    for (int i = 0; i < n; ++i) {
        l_inverse[i][i] = 1;
    }
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            s = 0;
            for (int k = 0; k < i; ++k) {
                s += l[i][k] * l_inverse[k][j];
            }
            l_inverse[i][j] = -s;
        }
    }
#if DEBUG
    printf("test l_inverse:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            s = 0;
            for (int k = 0; k < n; k++) {
                s += l[i][k] * l_inverse[k][j];
            }

            printf("%f ", s);
        }
        putchar('\n');
    }
#endif


    for (int i = 0; i < n; ++i) {
        u_inverse[i][i] = 1 / u[i][i];
    }
    for (int i = 0; i < n; ++i) {
        for (int j = i - 1; j >= 0; --j) {
            s = 0;
            for (int k = j + 1; k <= i; ++k) {
                s += u[j][k] * u_inverse[k][i];
            }
            u_inverse[j][i] = -s / u[j][j];
        }
    }
#if DEBUG
    printf("test u_inverse:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            s = 0;
            for (int k = 0; k < n; k++) {
                s += u[i][k] * u_inverse[k][j];
            }

            printf("%f ", s);
        }
        putchar('\n');
    }
#endif

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                a_inverse[i * n + j] += u_inverse[i][k] * l_inverse[k][j];
            }
        }
    }
#if DEBUG
    printf("test a:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            s = 0;
            for (int k = 0; k < n; k++) {
                s += getValue(i + 1, k + 1).real() * a_inverse[k * n + j];
            }

            printf("%f ", s);
        }
        putchar('\n');
    }
#endif
    Matrix ans(n, n, a_inverse);
    return ans;
}

Matrix Matrix::inverseOfMatrixImag() {
    int n = row;
    double l[n][n], u[n][n];
    double l_inverse[n][n], u_inverse[n][n];
    double_t a_inverse[n * n];

    memset(l, 0, sizeof(l));
    memset(u, 0, sizeof(u));
    memset(l_inverse, 0, sizeof(l_inverse));
    memset(u_inverse, 0, sizeof(u_inverse));
    memset(a_inverse, 0, sizeof(a_inverse));

    double s, t;
    for (int i = 0; i < n; ++i) {
        l[i][i] = 1;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            s = 0;
            for (int k = 0; k < i; ++k) {
                s += l[i][k] * u[k][j];
            }
            u[i][j] = getValue(i + 1, j + 1).imag() - s;
        }
        for (int j = 0; j < n; ++j) {
            s = 0;
            for (int k = 0; k < i; ++k) {
                s += l[j][k] * u[k][i];
            }
            l[j][i] = (getValue(j + 1, i + 1).imag() - s) / u[i][i];
        }
    }


    for (int i = 0; i < n; ++i) {
        l_inverse[i][i] = 1;
    }
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            s = 0;
            for (int k = 0; k < i; ++k) {
                s += l[i][k] * l_inverse[k][j];
            }
            l_inverse[i][j] = -s;
        }
    }
#if DEBUG
    printf("test l_inverse:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            s = 0;
            for (int k = 0; k < n; k++) {
                s += l[i][k] * l_inverse[k][j];
            }

            printf("%f ", s);
        }
        putchar('\n');
    }
#endif


    for (int i = 0; i < n; ++i) {
        u_inverse[i][i] = 1 / u[i][i];
    }
    for (int i = 0; i < n; ++i) {
        for (int j = i - 1; j >= 0; --j) {
            s = 0;
            for (int k = j + 1; k <= i; ++k) {
                s += u[j][k] * u_inverse[k][i];
            }
            u_inverse[j][i] = -s / u[j][j];
        }
    }
#if DEBUG
    printf("test u_inverse:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            s = 0;
            for (int k = 0; k < n; k++) {
                s += u[i][k] * u_inverse[k][j];
            }

            printf("%f ", s);
        }
        putchar('\n');
    }
#endif

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                a_inverse[i * n + j] += u_inverse[i][k] * l_inverse[k][j];
            }
        }
    }
#if DEBUG
    printf("test a:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            s = 0;
            for (int k = 0; k < n; k++) {
                s += getValue(i + 1, k + 1).imag() * a_inverse[k * n + j];
            }

            printf("%f ", s);
        }
        putchar('\n');
    }
#endif
    Matrix ans(n, n, a_inverse);
    return ans;
}

Matrix Matrix::inverse() {
    if (!inverse_Judge()) {
        cerr << "Invalid matrix inverse!" << endl;
        Matrix defaultMat(row, col, value);
        return defaultMat;
    }
    auto A_inverse = inverseOfMatrixReal();
    auto A = getFromRealPartOfMatrix();
    auto B = getFromImagPartOfMatrix();
    auto BmulA_inverse = B * A_inverse;
    auto BmulAmulB = BmulA_inverse * B;
    auto ApluBmulAmulB = A + BmulAmulB;
    auto ApluBmulAmulB_inverse = ApluBmulAmulB.inverseOfMatrixReal();
    auto A_inversemulB = A_inverse * B;
    auto A_inversemulBmulApluBmulAmulB_inverse = A_inversemulB * ApluBmulAmulB_inverse;
    double *a = new double[row * row];
    double *b = new double[row * row];
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < row; ++j) {
            a[i * row + j] = ApluBmulAmulB_inverse.getValue(i + 1, j + 1).real();
            b[i * row + j] = A_inversemulBmulApluBmulAmulB_inverse.getValue(i + 1, j + 1).real() * -1;
        }
    }
    Matrix matrix(row, row, a, b);
#if DEBUG
    //输出a的逆矩阵的实部矩阵
    printf("test a_inv_real:\n");
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < row; j++) {
            printf("%f ", matrix.getValue(i + 1, j + 1).real());
        }
        printf("\n");
    }
    //输出a的逆矩阵的虚部矩阵
    printf("test a_inv_imag:\n");
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < row; j++) {
            printf("%f ", matrix.getValue(i + 1, j + 1).imag());
        }
        printf("\n");
    }
#endif
    return matrix;
}

Matrix Matrix::getFromImagPartOfMatrix() {
    Matrix B(getRow(), getCol(), getImag());
    return B;
}

Matrix Matrix::getFromRealPartOfMatrix() {
    Matrix A(getRow(), getCol(), getReal());
    return A;
}


double *Matrix::getReal() {
    double *Reals = new double[row * col];
    for (int i = 0; i < row * col; ++i) {
        Reals[i] = value[i].real();
    }
    return Reals;
}

double *Matrix::getImag() {
    double *Imags = new double[row * col];
    for (int i = 0; i < row * col; ++i) {
        Imags[i] = value[i].imag();
    }
    return Imags;
}

bool Matrix::inverse_Judge() {
    if (row != col)
        return false;
    return true;
}

Matrix Matrix::conv2(Matrix &y) {
    int N1 = getRow();
    int N2 = y.getRow();
    int M1 = getCol();
    int M2 = y.getCol();
    int fullRow = N1 + N2 - 1, fullCol = M1 + M2 - 1;
    int i, j;
    int n, m;
    complex<double> *z = new complex<double>[fullRow * fullCol];
    complex<double> *z2 = new complex<double>[N1 * M1];
    for (i = 0; i < N1 + N2 - 1; i++)
        for (j = 0; j < M1 + M2 - 1; j++) {
            complex<double> temp = complex<double>(0, 0);
            for (m = 0; m < N1; m++)
                for (n = 0; n < M1; n++)
                    if ((i - m) >= 0 && (i - m) < N2 && (j - n) >= 0 && (j - n) < M2)
                        temp = temp + getValue(m + 1, n + 1) * y.getValue(i - m + 1, j - n + 1);
            z[i * fullCol + j] = temp;
        }
    for (i = 0; i < N1; i++) {
        for (j = 0; j < M1; j++) {
            z2[i * M1 + j] = z[(i + (N2 - 1) / 2) * fullCol + (j + (M2 - 1) / 2)];
        }
    }
#if DEBUG
    Matrix(fullRow, fullCol, z).print();
#endif
    return Matrix(N1, M1, z2);
}


