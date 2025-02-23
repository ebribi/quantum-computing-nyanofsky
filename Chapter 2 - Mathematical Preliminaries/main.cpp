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
            cout << real << " + " << imag << "i" << endl;
        }
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

Complex complexConjugate(Complex c){
    Complex res(c.getReal(),-c.getImag());
    return res;
}

double complexModulus(Complex c){
    return sqrt(pow(c.getReal(),2) + pow(c.getImag(),2));
}

Complex complexNormalized(Complex c){
    double mod = complexModulus(c);
    double real_norm = c.getReal() / mod;
    double imag_norm = c.getImag() / mod;
    Complex norm(real_norm,imag_norm);
    return norm;
}

double diffOfModSums(Complex c1, Complex c2){
    return complexModulus(c1) + complexModulus(c2) - complexModulus(c1 + c2);
}

int main(){


}