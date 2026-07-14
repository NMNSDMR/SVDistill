#include <iostream>
#include <vector>
#include <cmath>

using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;

struct SVDComponent {
    Vector u;
    double sigma;
    Vector v;
};

Matrix transpose(const Matrix& matrix){

    size_t rows = matrix.size();
    size_t cols = matrix[0].size();

    Matrix trans_matr(cols,Vector(rows));

    for(size_t i = 0; i < rows; i++){
        for(size_t j = 0; j < cols; j++){
            trans_matr[j][i] = matrix[i][j];
        }
    }
    return trans_matr;
}

Vector multiply(const Matrix& matrix, const Vector& vec){

    size_t rows = matrix.size();
    size_t cols = matrix[0].size();

    Vector res_vec(rows);

    for(size_t i = 0; i < rows; i++){
        for(size_t j = 0; j < cols; j++){
            res_vec[i] += vec[j] * matrix[i][j];
        }
    }

    return res_vec;

}

double getNorm(const Vector& v){

    double result = 0.0;

    for(size_t i = 0; i < v.size(); i++ ){
        result += v[i] * v[i];
    }

    return std::sqrt(result);
}



Vector normalize(const Vector& v){
    
    double norm = getNorm(v);
    size_t rows = v.size();

    Vector new_vec(rows);

    for(size_t i = 0; i < rows; i++){
        new_vec[i] = v[i] / norm;
    }
    return new_vec;
}

void print_matr(const Matrix& matrix){
    
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();
    
    for(size_t i = 0; i < rows; i++){
        for(size_t j = 0; j < cols; j++){
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void print_vec(const Vector& vec){

    size_t rows = vec.size();

    for(size_t i = 0; i < rows; i++){
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}

SVDComponent getFirstSVDComponent(const Matrix& A, int maxIterations = 100) {
    
    size_t rows = A.size();
    size_t cols = A[0].size();

    Matrix transpor_A = transpose(A);

    Vector v(cols,1.0); // случайный вектор 
    v = normalize(v);

    for(size_t i = 0; i < maxIterations; i++){

        Vector y = multiply(A,v); // поворачиваем матрицей А вектор v
        Vector x = multiply(transpor_A,y); // это как раз A^T * A * y. перемножение матриц усиливает растсяжение вдоль главной оси 
        v = normalize(x); // вектор х как то изменился в размерах, нам нужно ток направление 

    }
    
    Vector u_raw = multiply(A,v);
    double sigma = getNorm(u_raw);
    Vector u = normalize(u_raw);

    return {u,sigma,v};
}

struct SVDResult {
    Matrix U;
    Vector Sigma;
    Matrix VT;
};

SVDResult computeSVD(Matrix A, int k){
    size_t m = A.size();
    size_t n = A[0].size();

    Matrix U_res;
    Vector Sigma_res;
    Matrix V_res;

    for(size_t i = 0; i < k; i++){
        SVDComponent comp = getFirstSVDComponent(A,100);
        
        U_res.push_back(comp.u);
        Sigma_res.push_back(comp.sigma);
        V_res.push_back(comp.v);

        for(size_t row = 0; row < m; j++){
            for(size_t col = 0; col < n; k++){
                A[row][col] -= (comp.sigma * comp.u[row] * comp.v[col]);
            }
        }
    }
    
    Matrix U = transpose(U_res);
    Matrix VT = V_res;

    return {U,Sigma_res,VT};
}

Matrix compressChannel(const Matrix& A, int k) {
    SVDResult svd = computeSVD(A, k);
    size_t m = A.size();
    size_t n = A[0].size();
    
    Matrix compressed(m, Vector(n, 0.0));
    
    for(size_t i = 0; i < k; i++){
        for(size_t row = 0; row < m; row++){
            for(size_t col = 0; col < n; col++){
                compressed[row][col] += svd.Sigma[i] * svd.U[row][i] * svd.VT[i][col];
            }    
        }
    }
    
    return compressed;
}


void start_math_test(){
    Matrix test_matr = {{1,2},{3,4}};
    Vector test_vect = {3,4};
    std::cout << "///// Тестирование математики //////" << std::endl;

    std::cout << "-> Изначальный вектор и матрица:" << std::endl;
    std::cout << "-> Вектор: " << std::endl;
    print_vec(test_vect);

    std::cout << "-> Матрица: " << std::endl;
    print_matr(test_matr);

    std::cout << "-> Транспонированная матрица:" << std::endl;
    Matrix trans_matr = transpose(test_matr);
    print_matr(trans_matr);

    std::cout << "-> Перемножение матрицы на вектор:" << std::endl;
    Vector multiply_vec = multiply(test_matr,test_vect);
    print_vec(multiply_vec);

    double norm = getNorm(test_vect);
    std::cout << "-> Норма вектора: " << norm << std::endl;
    
    std::cout << "-> Нормализированный вектор:" << std::endl;
    Vector normalized_vec = normalize(test_vect);
    print_vec(normalized_vec); 

}

void start_power_iteration_test(){
    std::cout << "///// Тестирование степенного метода //////" << std::endl;
    Matrix test_svd_matrix = {{3.0,2.0},{2.0,3.0}};

    std::cout << " * Изначальная матрица:" << std::endl;
    print_matr(test_svd_matrix);

    SVDComponent comp = getFirstSVDComponent(test_svd_matrix,100);

    std::cout << " * Результаты SVD для первого приближения :" << std::endl;
    std::cout << " * Сигма (ожидаемые результаты около 5.0): " << comp.sigma << std::endl << std::endl;
    
    std::cout << " * Вектор U (ожидамые результаты ~0.707, ~0.707):" << std::endl;
    print_vec(comp.u);
    
    std::cout << " * Вектор V (ожидаемые результаты ~0.707, ~0.707):" << std::endl;
    print_vec(comp.v);

    // начинаем пытаться восстановить матрицу
    size_t rows = comp.u.size();
    size_t cols = comp.v.size();
    Matrix A_approx(rows, Vector(cols, 0.0));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            A_approx[i][j] = comp.sigma * comp.u[i] * comp.v[j];
        }
    }

    std::cout << "--- Восстановленная матрица (должна быть близка к оригиналу) ---" << std::endl;
    print_matr(A_approx);
}
int main(){

    start_math_test();
    start_power_iteration_test();
    
    return 0;
}