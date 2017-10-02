//
//  main.cpp
//  Fast Fourier Transform
//
//  Created by Niwatori on 2017/4/23.
//  Copyright © 2017 Niwatori. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <ctime>
using namespace std;

const double pi = acos(-1.0);

/*  Complex class  */
class Complex {
    double r, i;
public:
    Complex(double r_ = 0, double i_ = 0): r(r_), i(i_) {}
    double getReal() {return r;}
    double getImag() {return i;}
    Complex operator+(const Complex &c) {
        return Complex(r + c.r, i + c.i);
    }
    Complex operator-(const Complex &c) {
        return Complex(r - c.r, i - c.i);
    }
    Complex operator*(double c) {
        return Complex(r * c, i * c);
    }
    Complex operator*(const Complex &c) {
        return Complex(r * c.r - i * c.i, i * c.r + r * c.i);
    }
    Complex operator/(double c) {
        return Complex(r / c, i / c);
    }
    Complex operator/(const Complex &c) {
        return Complex(r * c.r + i * c.i, i * c.r - r * c.i) / (c.r * c.r + c.i * c.i);
    }
    friend ostream& operator<<(ostream &o, const Complex &c) {
        o << c.r;
        if (fabs(c.i) > 1e-10)
            o << (c.i >= 0 ? "+" : "") << c.i << "i";
        return o;
    }
};

/*  Binary inverse of int x of n digits  */
inline int rev(int x, int n) {
    int ans = 0;
    while (n--) {
        ans = (ans << 1) + (x & 1);
        x >>= 1;
    }
    return ans;
}

/*
 *  Fast Fourier Transform
 *  Input: a -- vector of length 2^n
 *         flag -- 1 for DFT and -1 for DFT inverse
 *  Output: c -- vector transformed from a
 */
void FFT(Complex *a, Complex *c, int n, int flag) {
    int N = 1 << n;
    for (int i = 0; i < N; ++i)
        c[i] = a[rev(i, n)];
    for (int p = 1; p <= n; ++p) {
        int stride = 1 << p;
        for (int start = 0; start < N; start += stride) {
            Complex w(1.0);
            Complex omega(cos(2 * pi / stride), -flag * sin(2 * pi / stride));
            for (int k = 0; k < stride / 2; ++k) {
                Complex x = c[start + k];
                Complex y = c[start + k + stride / 2];
                c[start + k] = x + w * y;
                c[start + k + stride / 2] = x - w * y;
                w = w * omega;
            }
        }
    }
    if (flag < 0)
        for (int i = 0; i < N; ++i)
            c[i] = c[i] / N;
}

void ReportTime(clock_t &start_time, clock_t &end_time) {
    cout << "Running time: "
         << static_cast<double> (end_time - start_time) / CLOCKS_PER_SEC * 1000
         << "ms" << endl;
}

inline double f(double t) {
    return exp(-t * t / 10) * (sin(2 * t) + 2 * cos(4 * t) + 0.4 * sin(t) * sin(50 * t));
}

int main(int argc, char **argv) {

{
    /* Problem A: Compute convolution */
    
    const int M = 500, Q = 200, N = 1024;
    Complex x[N] = {0}, x_[N], h[N] = {0}, h_[N];
    Complex y[N] = {0}, y_[N];
    for (int i = 1; i < M; ++i)
        x[i] = Complex(sin(i / 2.0));
    for (int i = 1; i < Q; ++i)
        h[i] = Complex(exp(1.0 / i));
    
    /* Computing with FFT */
    // clock_t start_time = clock();
    FFT(x, x_, 10, 1);
    FFT(h, h_, 10, 1);
    for (int i = 0; i < N; ++i)
        y_[i] = x_[i] * h_[i];
    FFT(y_, y, 10, -1);
    // clock_t end_time = clock();
    // ReportTime(start_time, end_time);
    
    /* Computing with brute force */
    // start_time = clock();
    memset(y, 0, sizeof(y));
    for (int i = 0; i < N; ++i)
        for (int j = 1; j < i; ++j)
            y[i] = y[i] + x[j] * h[i - j];
    // end_time = clock();
    // ReportTime(start_time, end_time);
}
    
{
    /* Problem B: Remove higher frequencies */
    
    const int N = 512, m = 18;
    Complex y[N] = {0}, y_[N] = {0}, y__[N] = {0};
    for (int k = 0; k <= 256; ++k)
        y[k] = Complex(f(2 * k * pi / 256));
    FFT(y, y_, 9, 1);
    for (int k = m; k <= 256 - m; ++k)
        y_[k] = 0;
    FFT(y_, y__, 9, -1);
    
    ofstream fout("data.csv");
    for (int k = 0; k <= 256; ++k)
        fout << y[k].getReal() << "," << y[k].getImag() << endl;
    for (int k = 0; k <= 256; ++k)
        fout << y__[k].getReal() << "," << y__[k].getImag() << endl;
}

{
    /* Problem C: Solve linear cyclic system with FFT */
    
    const int N = 1024;
    Complex c[N], b[N], x[N], c_[N], b_[N], x_[N];
    c[0] = Complex(4);
    c[1] = c[N - 1] = Complex(-1);
    for (int i = 0; i < N; ++i)
        b[i] = Complex(1);
    
    FFT(b, b_, 10, 1);
    FFT(c, c_, 10, 1);
    for (int i = 0; i < N; ++i)
        x_[i] = b_[i] / c_[i];
    FFT(x_, x, 10, -1);
    for (int i = 0; i < N; ++i)
        cout << x[i] << endl;
}
    
{
    /* Problem D: Solve differential equations with FFT */
    
    const int n = 8;  /* n = 4, 6, 8 */
    const int N = 1 << n;
    const double h = pi / (3 * N);
    Complex f[N], f_[N], u[N], u_[N];
    for (int j = 0; j < N; ++j)
        f[j] = 3 * cos(6 * j * h);
    FFT(f, f_, n, 1);
    
    for (int j = 0; j < N; ++j)
        u_[j] = (f_[j] * h * h) / (Complex(cos(2 * pi * j / N), sin(2 * pi * j / N)) * (1 + h)
                                    + Complex(cos(2 * pi * j / N), -sin(2 * pi * j / N)) * (1 - h)
                                    - 2 * (1 - h * h));
    FFT(u_, u, n, -1);
    for (int j = 0; j < N; ++j)
        cout << u[j] << "," << endl;
}
    
    return 0;
}
