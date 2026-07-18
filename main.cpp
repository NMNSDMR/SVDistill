#include <iostream>
#include <vector>
#include <cmath>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h" // ЧИТАТЬ картинки

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h" // ЗАПИСЫВАТЬ картинки

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

        for(size_t row = 0; row < m; row++){
            for(size_t col = 0; col < n; col++){
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

unsigned char clamp(double val) { // пишем это чтобы если значение вылетело за разрешенное значение там не получился какой то цвет другой
        if (val < 0.0) return 0;
        if (val > 255.0) return 255;
        return static_cast<unsigned char>(std::round(val));
    }

struct ImageRGB {
    int width = 0;
    int height = 0;
    
    Matrix R;
    Matrix G;
    Matrix B;

    bool load(const std::string& filename) {
        int channels;
        unsigned char* data = stbi_load(filename.c_str(), &width, &height, &channels, 3);
        
        if (!data) {
            std::cerr << "Не удалось загрузить картинку: " << filename << std::endl;
            return false;
        }
        R = Matrix(height, Vector(width, 0.0));
        G = Matrix(height, Vector(width, 0.0));
        B = Matrix(height, Vector(width, 0.0));
        
        size_t index = 0;
        for(size_t i = 0; i < height;i++){
            for(size_t j = 0; j < width;j++){
                R[i][j] = data[index++];
                G[i][j] = data[index++];
                B[i][j] = data[index++];
            }
        }
        stbi_image_free(data); 
        return true;
    }

    void save(const std::string& filename) {

        std::vector<unsigned char> data(width * height * 3);
        size_t index = 0;

        for(size_t i = 0; i < height; i++){
            for(size_t j = 0; j < width; j++){
                data[index++] = clamp(R[i][j]);
                data[index++] = clamp(G[i][j]);
                data[index++] = clamp(B[i][j]);
            }
        }

        stbi_write_png(filename.c_str(), width, height, 3, data.data(), width * 3);
    }
};

int main(){

    ImageRGB img;
    std::string input_file = "imdir/input.jpg"; 
    
    if (!img.load(input_file)) {
        return -1;
    }
    
    std::cout << "Размер изображения: " << img.width << "x" << img.height << std::endl;

    int k = 20; // можно менять для повышения/понижения ранга

    img.R = compressChannel(img.R, k);
    img.G = compressChannel(img.G, k);
    img.B = compressChannel(img.B, k);

    std::string output_file = "imdir/output_k" + std::to_string(k) + ".png";
    img.save(output_file);
    return 0;
}