#include <iostream>
#include <vector>
#include <cmath>

using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;

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

int main(){
    Matrix test_matr = {{1,2},{3,4}};
    Vector test_vect = {3,4};
    std::cout << "///// Тестирование //////" << std::endl;

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

    return 0;
} 