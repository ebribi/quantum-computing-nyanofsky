#include <vector>
#include <iostream>
#include <cmath>

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

        Complex operator/(double d){
            Complex res;
            res.real = real / d;
            res.imag = imag / d;
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

// Programming Exercise 2.1.20(i) Write a function that accepts a complex matrix and returns its conjugate.

vector<vector<Complex>> matrixConj(vector<vector<Complex>> matrix){
    int rows = matrix.size();
    int cols = matrix[0].size();

    vector<vector<Complex>> conj(rows, vector<Complex>(cols,Complex(0,0)));

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            conj[i][j] = complexConjugate(matrix[i][j]);
        }
    }

    return conj;
}

// Programming Exercise 2.1.20(ii) Write a function that accepts a complex matrix and returns the dagger of that matrix.

vector<vector<Complex>> matrixDagger(vector<vector<Complex>> matrix){
    return matrixConj(matrixTranspose(matrix));
}

// Programming Exercise 2.1.20(iii) Write a function that accepts a complex matrix and returns the trace of the matrix.

Complex matrixTrace(vector<vector<Complex>> matrix){
    int rows = matrix.size();
    int cols = matrix[0].size();

    if (rows != cols){
        cout << "Error: Trace is undefined." << endl;
        return {};
    }

    Complex trace;

    for (int i = 0; i < matrix.size(); i++){
        trace = trace + matrix[i][i];
    }

    return trace;
    
}

// Programming Exercise 2.2.18(i) Write a function that accepts a vector and returns the norm of the vector

double vectorNorm(vector<Complex> vec){
    
    Complex norm(0,0);    
    
    for (int i = 0; i < vec.size(); i++){
        norm = norm + (vec[i] * complexConjugate(vec[i]));
    }

    return sqrt(norm.getReal());

}


// Programming Exercise 2.2.18(ii) Write a function that accepts a vector and returns the normalized form of the vector
vector<Complex> vectorNormalized(vector<Complex> vec){
    
    double norm = vectorNorm(vec);
    vector<Complex> vec_normalized(vec.size(),Complex(0,0));  
    
    for (int i = 0; i < vec.size(); i++){
        vec_normalized[i] = vec[i] / norm;
    }

    return vec_normalized;

}

// Programming Exercise 2.2.18(iii) Write three functions which accept a basis (a string or set of column vectors). The three functions should determine if the basis is normal, orthogonal, and orthonormal, respectively.

bool isNormal(vector<vector<Complex>> basis){

    for (int i = 0; i < basis.size(); i++){
        if (vectorNorm(basis[i]) != 1)
            return false;
    }

    return true;
}

double innerProduct(vector<Complex> vec1, vector<Complex> vec2){

    if (vec1.size() != vec2.size()){
        cout << "Effor: Sizes incompatible for inner product." << endl;
        return {};
    }

    Complex sum(0,0);
    for (int i = 0; i < vec1.size(); i++){
        sum = sum + (vec1[i] * vec2[i]);
    }
    return sqrt(sum.getReal());    
}

bool isOrthogonal(vector<vector<Complex>> basis){

    for (int i = 0; i < basis.size(); i++){
        for (int j = 0; j < basis.size(); j++){   
            if (i != j && innerProduct(basis[i],basis[j]) != 0)
                return false;
        }
    }

    return true;
}

bool isOrthonormal(vector<vector<Complex>> basis){
    return isNormal(basis) && isOrthogonal(basis);
}

