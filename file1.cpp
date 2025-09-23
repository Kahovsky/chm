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

/*Matrix readMatrix(const string& filename) {
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

matrix.AN.resize(n, vector<double>(k + 1));

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
}*/

vector<vector<double>> readMatrix(const string& filename) {
    ifstream file(filename);
    int n;
    file >> n;

    vector<vector<double>> matrix;
    matrix.resize(n, vector<double>(n));


    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            file >> matrix[i][j];
        }
    }
    return matrix;
}

Matrix denseToBand(const vector<vector<double>>& dense) {
    int n = dense.size();
    int k = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (dense[i][j] != 0) {
                int dist = abs(i - j);
                if (dist > k) k = dist;
            }
        }
    }

    Matrix band;
    band.n = n;
    band.k = k;

    band.AD.resize(n);
    band.AN.resize(n, vector<double>(k, 0.0));
    band.AL.resize(n, vector<double>(k, 0.0));

    for (int i = 0; i < n; i++) {
        band.AD[i] = dense[i][i];
        for (int l = 1; l <= k; l++) {
            if (i - l >= 0) {
                band.AL[i][l - 1] = dense[i][i - l];
                band.AN[i][l - 1] = dense[i - l][i];
            }
        }
    }

    return band;
}

vector<double> readVector(const string& filename) {
    ifstream file(filename);
    int n;
    file >> n;

    vector<double> vec(n);
    for (int i = 0; i < n; i++) {
        file >> vec[i];
    }
    file.close();
    return vec;
}

double getelem(const Matrix& A, int i, int j) {
    if (i == j) {
        return A.AD[i];
    }
    else if (j > i) { // Верхний треугольник - элемент над диагональю
        int dist = j - i - 1;
        if (dist >= 0 && dist < A.k) {
            return A.AN[i][dist];  // Берем из строки i
        }
    }
    else { // Нижний треугольник - элемент под диагональю
        int dist = i - j - 1;
        if (dist >= 0 && dist < A.k) {
            return A.AL[i][dist];  // Берем из строки i
        }
    }
    return 0.0;
}

void printLUMatrix(const Matrix& matrix) {
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
                cout << "@" << getelem(L, i, p) << " " << getelem(U, p, j) << "\n";
            }

            // L[i][j] = (A[i][j] - sum) / U[j][j]
            double a_ij = getelem(A, i, j);
            double l_ij = (a_ij - sum);
            cout << "#" << a_ij << " " << l_ij << "\n";

            if (i - j - 1 < k && i - j - 1 >= 0) {
                L.AL[i][i - j - 1] = l_ij;
                cout << "%" << getelem(L, i, j) << "\n";
            }
        }

        // Вычисляем диагональный элемент L
        double sum = 0.0;
        for (int p = max(0, i - k); p < i; p++) {
            sum += getelem(L, i, p) * getelem(U, p, i);
        }
        L.AD[i] = getelem(A, i, i) - sum;

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

        
    }
}

void printDenseMatrix(const vector<vector<double>>& matrix) {
    int n = matrix.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << "\n";
    }
}

void checkLU(const Matrix& A, const Matrix& L, const Matrix& U) {
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
    vector<vector<double>> matrix = readMatrix("matrix.txt");
    printDenseMatrix(matrix);

    cout << "*****" << "\n";
    
    Matrix matrixBand = denseToBand(matrix);
    printLUMatrix(matrixBand);
    Matrix L, U;
    LUdec(matrixBand, L, U);
    printLUMatrix(L);
    cout << "*****" << "\n";
    printLUMatrix(U);

    checkLU(matrixBand, L, U);
}
