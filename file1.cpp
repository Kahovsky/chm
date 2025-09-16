#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

struct Matrix {
    int n, k;
    vector<double> AD;
    vector<vector<double>> AN;
    vector<vector<double>> AL;

};

Matrix readMatrix(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Не удалось открыть файл: " + filename);
    }
    int k, n;

    Matrix matrix;
    file >> n >> k;
    matrix.n = n;
    matrix.k = k;

    
    matrix.AD.resize(n);

    for (int i = 0; i < n; i++) {
        file >> matrix.AD[i];
    }

    matrix.AN.resize(n, vector<double>(k+1));

    //верхняя часть
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            file >> matrix.AN[i][j];
            
        }

    }

    matrix.AL.resize(n, vector<double>(k));
    //нижняя часть
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            file >> matrix.AL[i][j];
        }

    }
    return matrix;
}

double getelem(const Matrix& A, int i, int j) {
    if (i == j) {
        return A.AD[i];
    }
    else if (j > i) { // Верхний треугольник
        int band_index = A.k - j - i - 2;
        if (band_index < A.k && band_index >= 0) {
            return A.AN[i][band_index];
        }
    }
    else { // Нижний треугольник
        int band_index = A.k - j - i - 2;
        if (band_index < A.k && band_index >= 0) {
            return A.AL[j][band_index];
        }
    }
    return 0.0;
}

void printMatrix(const Matrix& matrix) {
    cout << matrix.n << " " << matrix.k << "\n";
    for (int i = 0; i < matrix.n; i++) {
        cout << matrix.AD[i] << " ";
    }
    cout << "\n";
    //верхняя часть
    for (int i = 0; i < matrix.n; i++) {
        for (int j = 0; j < matrix.k; j++) {
            cout <<  matrix.AN[i][j] << " ";
        }
        cout << "\n";

    }
    cout << "\n";
    //нижняя часть
    for (int i = 0; i < matrix.n; i++) {
        for (int j = 0; j < matrix.k; j++) {
            cout << matrix.AL[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

void LUdec(const Matrix& A, Matrix& L, Matrix& U) {
    int n = A.n;
    int k = A.k;

    L = A;
    U = A;

    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            L.AN[i][j] = 0.0; // Верхняя часть L = 0
            U.AL[i][j] = 0.0; // Нижняя часть U = 0
        }
    }

    // Диагональ U = 1
    for (int i = 0; i < n; i++) {
        U.AD[i] = 1.0;
    }
   


    for (int i = 0; i < n; i++) {
        // Вычисляем элементы L для текущей строки
        for (int j = max(0, i - k); j < i; j++) {
            double sum = 0.0;
            for (int p = max(0, max(i - k, j - k)); p < j; p++) {
                sum += getelem(L, i, p) * getelem(U, p, j);
            }

            // L[i][j] = (A[i][j] - sum) / U[j][j]
            double a_ij = getelem(A, i, j);
            double l_ij = (a_ij - sum);

            if (i - j - 1 < k && i - j - 1 >= 0) {
                L.AL[j][i - j - 1] = l_ij;
            }
        }

        // Вычисляем элементы U для текущей строки
        for (int j = i + 1; j < min(n, i + k + 1); j++) {
            double sum = 0.0;
            for (int p = max(0, i - k); p < i; p++) {
                sum += getelem(L, i, p) * getelem(U, p, j);
            }

            // U[i][j] = (A[i][j] - sum) / L[i][i]
            double a_ij = getelem(A, i, j);
            double u_ij = (a_ij - sum) / getelem(L, i, i);

            if (j - i - 1 < k && j - i - 1 >= 0) {
                U.AN[i][j - i - 1] = u_ij;
            }
        }

        // Вычисляем диагональный элемент L
        double sum = 0.0;
        for (int p = max(0, i - k); p < i; p++) {
            sum += getelem(L, i, p) * getelem(U, p, i);
        }
        L.AD[i] = getelem(A, i, i) - sum;
    }
}


// ТОЛКО ДЛЯ ПРОВЕРОК!!! спизженно и не обязательно, мы слишком тупые что бы самим до такого догадаться
void checkLU(const Matrix & A, const Matrix & L, const Matrix & U) {
    int n = A.n;
    cout << "Check A = L * U:\n";

    double max_error = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double a_ij = getelem(A, i, j);

            // Вычисляем (L * U)[i][j]
            double lu_ij = 0.0;
            for (int p = 0; p < n; p++) {
                lu_ij += getelem(L, i, p) * getelem(U, p, j);
            }

            double error = abs(a_ij - lu_ij);
            max_error = max(max_error, error);

            if (error > 1e-6) {
                cout << "A[" << i << "][" << j << "] = " << a_ij
                    << ", (L*U)[" << i << "][" << j << "] = " << lu_ij
                    << ", err = " << error << "\n";
            }
        }
    }

    cout << "max err: " << max_error << "\n";
}


int main() {
    Matrix matrix = readMatrix("matrix.txt");
    printMatrix(matrix);
    cout << "*****" << "\n";
    Matrix L, U;
    LUdec(matrix, L, U);
    printMatrix(L);
    cout << "*****" << "\n";
    printMatrix(U);
    checkLU(matrix, L, U);
}
