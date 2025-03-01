#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

class Complex {
    private:
        double real;
        double imag;
    public:
        Complex(){
            real = 0;
            imag = 0;
        };
        Complex(double real, double imag){
            this->real = real;
            this->imag = imag;
        }
        double getReal(){return real;}
        double getImag(){return imag;}
        void display(){
            cout << real << " + " << imag << "i";
        }

        // Programming Exercise 2.1.6 Write four functions that perform the addition, subtraction, multiplication, and division operations on complex numbers. The functions should accept and return the general form of the complex numbers. When diving, make sure the divisor is not zero.
        
        Complex operator+(Complex const& c){
            Complex res;
            res.real = real + c.real;
            res.imag = imag + c.imag;
            return res;
        }
        Complex operator-(Complex const& c){
            Complex res;
            res.real = real - c.real;
            res.imag = imag - c.imag;
            return res;
        }
        Complex operator*(Complex const& c){
            Complex res;
            res.real = (real * c.real) - (imag * c.imag);
            res.imag = (real * c.imag) + (imag * c.real);
            return res;
        }

        Complex operator/(Complex const& c){
            Complex res;
            res.real = ((real * c.real) + (imag * c.imag))/(pow(c.real,2) + pow(c.imag,2));
            res.imag = ((c.real * imag) - (c.imag * real))/(pow(c.real,2) + pow(c.imag,2));
            return res;
        }
};

// Programming Exercise 2.1.12(i) Write a function that accepts a complex number and returns its complex conjugate

Complex complexConjugate(Complex c){
    Complex res(c.getReal(),-c.getImag());
    return res;
}

// Programming Exercise 2.1.12(ii) Write a function that accepts a complex number and returns its modulus

double complexModulus(Complex c){
    return sqrt(pow(c.getReal(),2) + pow(c.getImag(),2));
}

// Programming Exercise 2.1.12(iii) Write a function that accepts a complex number and returns the normal form of the complex number

Complex complexNormalized(Complex c){
    double mod = complexModulus(c);
    double real_norm = c.getReal() / mod;
    double imag_norm = c.getImag() / mod;
    Complex norm(real_norm,imag_norm);
    return norm;
}

// Programming Exercise 2.1.12(iv) Write a function that accepts two moduli and the modulus of the sum. That is, the functions accepts c and c’ and returns |c| + |c’| - |c - c’|

double diffOfModSums(Complex c1, Complex c2){
    return complexModulus(c1) + complexModulus(c2) - complexModulus(c1 + c2);
}

// Programming Exercise 2.1.15(i) Write a function that accepts two complex matrices that are of the same size. The function should return a matrix which is the sum of the matrices.

vector<vector<Complex>> matrixSum(vector<vector<Complex>> matrix1, vector<vector<Complex>> matrix2){
    int rows1 = matrix1.size();
    int cols1 = matrix1[0].size();
    int rows2 = matrix2.size();
    int cols2 = matrix2[0].size();    

    if (rows1 != rows2 && cols1 != cols2){
        cout << "Error: Incompatible matrix dimensions for multiplication." << endl;
        return {};    
    }

    vector<vector<Complex>> sum(rows1, vector<Complex>(cols1,Complex(0,0)));

    for (int i = 0; i < rows1; i++){
        for (int j = 0; j < rows1; j++){
            sum[i][j] = sum[i][j] + matrix1[i][j] + matrix2[i][j];        
        }    
    }

    return sum;
}

// Programming Exercise 2.1.15(ii) Write a function that accepts a complex number and a complex matrix. The function should return a matrix which is the scalar multiplication of the complex number with the matrix.
 
vector<vector<Complex>> matrixScalar(vector<vector<Complex>> matrix, Complex r){
    int rows = matrix.size();
    int cols = matrix[0].size();

    vector<vector<Complex>> product(rows, vector<Complex>(cols,Complex(0,0)));
    
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            product[i][j] = product[i][j] + (r * matrix[i][j]);      
        }
    } 

    return product;

}

// Programming Exercise 2.1.15(iii) Write a function that accepts two matrices with complex entries such that the number of columns in the first matrix is the same as the number of rows in the second matrix. The function should return a multiplication of the two matrices.

vector<vector<Complex>> matrixMultiply(vector<vector<Complex>> matrix1, vector<vector<Complex>> matrix2){
    int rows1 = matrix1.size();
    int cols1 = matrix1[0].size();
    int rows2 = matrix2.size();
    int cols2 = matrix2[0].size();

    if (cols1 != rows2){
        cout << "Error: Incompatible matrix dimensions for multiplication." << endl;
        return {};    
    }

    vector<vector<Complex>> product(rows1, vector<Complex>(cols2,Complex(0,0)));
    
    for (int i = 0; i < rows1; i++){
        for (int j = 0; j < cols2; j++){
            for (int k = 0; k < cols1; k++){
                product[i][j] = product[i][j] + (matrix1[i][k] * matrix2[k][j]);            
            }        
        }
    }

    return product;
}

// Programming Exercise 2.1.15(iv) Write a function that accepts a complex matrix. The function should return the transpose of the matrix.

vector<vector<Complex>> matrixTranspose(vector<vector<Complex>> matrix){
    int cols = matrix.size();
    int rows = matrix[0].size();

    vector<vector<Complex>> transpose(rows, vector<Complex>(cols,Complex(0,0)));

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
              transpose[j][i] = matrix[i][j];
        }
    }

    return transpose;
}

int main(){

    return 0;

}
