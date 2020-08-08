
#include "pch.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include  <iomanip>

using namespace std;

int N = 20;

double a = 0.0, b = 1.0; // [a, b]
double alfa1 = -1.0, alfa2 = 1.0, A = 0.6;                   // alfa1 * y(a) + alfa2 * y'(a) = A
double beta1 = 1.0, beta2 = 1.0, B = 4.0*exp(3.0) + exp(4.0); // beta1 * y(b) + beta2 * y'(b) = B

//Analytic solution for given function
double ExactSolution(double x)
{
    return exp(-x) + exp(3 * x) + 0.2*exp(4 * x);
}

//Next functions returns appropriate parts of equation
// y'' + p*y' + q*y + f = 0
double p(double x)
{
    return -2;
}

double q(double x)
{
    return -3;
}

double f(double x)
{
    return exp(4 * x);
}

//Functions ai, bi, ci, di coefficients of a finite difference scheme
double ai(int n, double h, double pi, int i) {
    if (i < n) {
        return 1 - ((h*pi) / 2);
    }
    else {
        return -beta2;
    }
}

double bi(int n, double h, double qi, int i) {
    if (i < n && i > 0) {
        return h * h * qi - 2;
    }
    else if (i == 0) {
        return h * alfa1 - alfa2;
    }
    else if (i == n) {
        return h * beta1 + beta2;
    }

}

double ci(int n, double h, double pi, int i) {
    if (i < n && i > 0) {
        return 1 + ((h*pi) / 2);
    }
    else if (i == 0) {
        return alfa2;
    }

}

double di(int n, double h, double fi, int i) {
    if (i < n && i > 0) {
        return h * h * fi;
    }
    else if (i == 0) {
        return h * A;
    }
    else if (i == n) {
        return h * B;
    }
}


//Second possible implementation for Tridiagonal matrix algorithm
//You just need to uncomment code below and ### in code
/*
void forward(int n, double **a, double *b)
{
    double v;
    for (int k = 0, i, j, im; k < n - 1; k++)
    {
        im = k;
        for (i = k + 1; i < n; i++)
        {
            if (fabs(a[im][k]) < fabs(a[i][k]))
            {
                im = i;
            }
        }
        if (im != k)
        {
            for (j = 0; j < n; j++)
            {
                v = a[im][j];
                a[im][j] = a[k][j];
                a[k][j] = v;
            }
            v = b[im];
            b[im] = b[k];
            b[k] = v;
        }
        for (i = k + 1; i < n; i++)
        {
            v = 1.0*a[i][k] / a[k][k];
            a[i][k] = 0;
            b[i] = b[i] - v * b[k];
            if (v != 0)
                for (j = k + 1; j < n; j++)
                {
                    a[i][j] = a[i][j] - v * a[k][j];
                }
        }
    }
}

void reverse(int n, double **a, double *b, double *x)
{
    double s = 0;
    x[n - 1] = 1.0*b[n - 1] / a[n - 1][n - 1];
    for (int i = n - 2, j; 0 <= i; i--)
    {
        s = 0;
        for (j = i + 1; j < n; j++)
        {
            s = s + a[i][j] * x[j];
        }
        x[i] = 1.0*(b[i] - s) / a[i][i];
    }
}
*/

//Call Gauss method (Gaussian Elimination) to solve A*x = b system
vector<double> GaussianEliminationSSP(vector<vector<double>>& matrix, int n) {
    //Initialize vector s
    vector<double> s;
    for (int i = 0; i < n; i++) {
        //Consider all entries except the last one, which is storing value of vector b
        double max = *max_element(matrix[i].begin(), matrix[i].end() - 1);
        s.push_back(max);
    }

    for (int i = 0; i < n - 1; i++) {
        //Finding the row which the relative pivot element is the largest.
        double max = matrix[i][i] / s[i];
        int max_index = i;
        for (int j = i + 1; j < n; j++) {
            if (max < matrix[j][i] / s[j]) {
                max = matrix[j][i] / s[j];
                max_index = j;
            }
        }
        //Switching two rows, as well as two elements in vector s
        swap(matrix[i], matrix[max_index]);
        swap(s[i], s[max_index]);

        //Do the forward elimination
        for (int j = i + 1; j < n; j++) {
            double factor = matrix[j][i] / matrix[i][i];
            matrix[j][i] = 0;
            for (int k = i + 1; k <= n; k++) {
                matrix[j][k] = matrix[j][k] - factor * matrix[i][k];
            }
        }
    }
    //Computing the final results
    vector<double> result;
    for (int i = n - 1; i >= 0; --i) {
        double x = matrix[i][n] / matrix[i][i];
        for (int j = 0; j < i; j++) {
            matrix[j][n] -= x * matrix[j][i];
        }
        result.insert(result.begin(), x);
    }

    return result;
}

double *finite(int n) {
    //define matrix
    double **matrA = new double*[n + 1];
    for (int i = 0; i < n + 1; i++) {
        matrA[i] = new double[n + 1];
    }
    double *matrB = new double[n + 1];
    double *X = new double[n + 1]; //vector X[n+1]
    //eval h
    double h = (b - a) / n;
    
    //eval pi, qi, fi
    double *mp = new double[n];
    double *mq = new double[n];
    double *mf = new double[n];

    double x = a;
    for (int i = 0; i < n; i++) {
        mp[i] = p(x);
        mq[i] = q(x);
        mf[i] = f(x);
        x += h;
    }

    // zeroing
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n; j++) {
            matrA[i][j] = 0;
        }
    }

    // fill matrix free element
    for (int i = 0; i <= n; i++) {
        matrB[i] = di(n, h, mf[i], i);
    }
    // fill 0 row
    matrA[0][0] = bi(n, h, mq[0], 0);
    matrA[0][1] = ci(n, h, mp[0], 0);

    //Find coefficients of a finite difference scheme (matrix A)
    //and column vector of f(x) (right side - matrix B) and represent 
    //differential equation as a system of linear equations: A*x = b

    // fill i row
    for (int i = 1; i <= n - 1; i++) {
        matrA[i][i - 1] = ai(n, h, mp[i], i);
        matrA[i][i] = bi(n, h, mq[i], i);
        matrA[i][i + 1] = ci(n, h, mp[i], i);
    }
    // fill n row
    matrA[n][n - 1] = ai(n, h, mp[n], n);
    matrA[n][n] = bi(n, h, mq[n], n);
    /*
        cout << "System of linear equations: \n";

        for(int i = 0; i <= n; i++){
            for(int j = 0; j <= n; j++){
                cout << fixed << setprecision(3) <<  matrA[i][j] << "\t";
            }
            cout << " |= "<< matrB[i] << endl;
        }
    */

    ////////
    //### Second possible implementation for Tridiagonal matrix algorithm
    //You just need to uncomment code below and comment @@@
    /*
    //forward(n + 1, matrA, matrB);
    //reverse(n + 1, matrA, matrB, X);
    */
    ////////

    //Copy matrix A and b to std::vector
    vector<vector<double>> matrix;
    for (int i = 0; i <= n; i++) {
        matrix.push_back({});
        for (int j = 0; j <= n; j++) {
            matrix[i].push_back(matrA[i][j]);
        }
        matrix[i].push_back(matrB[i]);
    }

    //Call Gauss method (Gaussian Elimination) to solve A*x = b system
    vector<double> result = GaussianEliminationSSP(matrix, n+1);

    for (int i = 0; i <= n; ++i) {
        X[i] = result[i];
    }

    ////////

    delete[] mp;
    delete[] mq;
    delete[] mf;

    return X;
}

int main()
{
    cout << "Count of mesh intervals in space domain (N): " << endl;
    cin >> N;
    cout << "Set left boundary a in space domain [a, b]: " << endl;
    cin >> a;
    cout << "Set right boundary b in space domain [a, b]: " << endl;
    cin >> b;

    cout << "Tridiagonal matrix algorithm \n"; // Progonki method 

    vector<double> x(N + 1);
    vector<double> y(N + 1);
    vector<double> A(N + 1);
    vector<double> B(N + 1);
    vector<double> C(N + 1);
    vector<double> F(N + 1);
    vector<double> aa(N + 1);
    vector<double> bb(N + 1);

    double h, xx;
    int n = N;
    int i;

    //disrete step size in space domain [a, b]
    h = (b - a) / N;

    //disrete position in space domain [a, b]
    for (i = 0; i <= N; i++)
    {
        x[i] = a + h * i;
    }

    //coefficients of a finite difference scheme
    for (i = 0; i < N; i++)
    {
        A[i] = 1.0 / (h*h) - p(x[i]) / (2.0 * h);
        C[i] = 1.0 / (h*h) + p(x[i]) / (2.0 * h);
        B[i] = -2.0 / (h*h) + q(x[i]);
        F[i] = f(x[i]);
    }

    B[0] = -h - 1.0;
    C[0] = 1.0;
    F[0] = 0.6*h;
    B[N] = 1.0 + h;
    A[N] = -1.0;
    F[N] = h * (4.0 * exp(3.0) + exp(4.0));
    aa[0] = -C[0] / B[0];
    bb[0] = F[0] / B[0];

    //Forward phase of Tridiagonal matrix algorithm
    for (i = 1; i <= N; i++)
    {
        aa[i] = -C[i] / (A[i] * aa[i - 1] + B[i]);
        bb[i] = (F[i] - A[i] * bb[i - 1]) / (A[i] * aa[i - 1] + B[i]);
    }

    //Backward phase of Tridiagonal matrix algorithm
    //Using recurrent equation
    y[n] = (F[n] - bb[n - 1] * A[n]) / (B[n] + aa[n - 1] * A[n]);

    for (i = N - 1; i >= 0; i--)
    {
        y[i] = aa[i] * y[i + 1] + bb[i];
    }

    cout << "  x\t  y" << endl;
    for (i = 0; i <= N; i++)
    {
        cout << fixed << setprecision(5) << x[i] << '\t' << y[i] << endl;
    }


    cout << "_______________________________" << endl;
    cout << "_______________________________" << endl;
    cout << "_______________________________" << endl;

    cout << "Finite difference method \n"; // Progonki method 


    int flag = 0;
    double eps = 0.0001;
    int counter = 0;

    double *X = new double[n + 1];
    //Call implicit FDM solver 
    X = finite(n);
    double xa = a;
    //Show results
    cout << "  x\t  y" << endl;
    for (int i = 0; i <= n; i++) {
        cout << fixed << setprecision(5) << xa << '\t' << X[i] << endl;
        xa += (b - a) / n;
    }
    delete[] X;


    cout << "_______________________________" << endl;
    cout << "_______________________________" << endl;
    cout << "_______________________________" << endl;
    cout << "Exact Solution \n";
    cout << "  x\t  y" << endl;
    for (i = 0; i <= N; i++)
    {
        cout << fixed << setprecision(5) << x[i] << '\t' << ExactSolution(x[i]) << endl;
    }

}
