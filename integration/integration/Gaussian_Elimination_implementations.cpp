#include <iostream>
#include <vector>
#include <cmath>

using namespace std;




// Function to solve the Ax=b system
/*
    K: a reference to a vector of doubles representing the coefficients 
    of the system of linear equations to be solved.

    b: a reference to a vector of doubles representing the constant terms 
    of the system of linear equations.

    x: a reference to a vector of doubles representing the unknowns of the system 
    of linear equations, which will be computed by the function.
*/


void GaussianElimination(vector<double>& K, vector<double>& b, vector<double>& x)
{
    const int NEQ = b.size();
    double p, t, pivot;
    int ii, j, k, kpivot;

    // Forward elimination.
    for (j = 1; j < (NEQ - 2); j++)
    {
        // Find pivot.
        pivot = fabs(K[j * NEQ + j]);
        kpivot = j;

        for (ii = (j + 1); ii < NEQ; ii++)
        {
            if (fabs(K[ii * NEQ + j]) > pivot)
            {
                pivot = fabs(K[ii * NEQ + j]);
                kpivot = ii;
            }   // if
        }   // for

        // Exchange rows "j" and "kpivot" of [K:b].
        for (ii = 0; ii < NEQ; ii++)
        {
            t = K[j * NEQ + ii];
            K[j * NEQ + ii] = K[kpivot * NEQ + ii];
            K[kpivot * NEQ + ii] = t;
        }   // for
        t = b[j];
        b[j] = b[kpivot];
        b[kpivot] = t;

        // Proceed.
        for (ii = (j + 1); ii < NEQ; ii++)
        {
            p = K[ii * NEQ + j] / K[j * NEQ + j];
            K[ii * NEQ + j] = 0.0;
            b[ii] -= p * b[j];
            for (k = (j + 1); k < NEQ; k++)
            {
                K[ii * NEQ + k] -= p * K[j * NEQ + k];
            }   // for
        }   // for
    }   // for

    // Back substitution.
    for (ii = (NEQ - 1); ii >= 0; ii--)
    {
        x[ii] = b[ii];
        for (j = (ii + 1); j < NEQ; j++)
        {
            x[ii] -= K[ii * NEQ + j] * x[j];
        }   // for
        x[ii] /= K[ii * NEQ + ii];
    }   // for
}   

vector<double> gauss_elimination(vector<vector<double>>& A, vector<double>& b) {
    int n = A.size();

    // Gaussian elimination with partial pivoting and row scaling
    vector<int> index(n);
    vector<double> scale(n);
    for(int i = 0; i < n; i++) {
        index[i] = i;
        scale[i] = 0;
        double max_val = 0;
        for(int j = 0; j < n; j++) {
            if(abs(A[i][j]) > max_val) {
                max_val = abs(A[i][j]);
            }
        }
        if(max_val == 0) {
            throw runtime_error("Matrix is singular.");
        }
        scale[i] = 1 / max_val;
    }
    for(int k = 0; k < n-1; k++) {
        double max_val = 0;
        int max_index = k;
        for(int i = k; i < n; i++) {
            double val = scale[index[i]] * abs(A[index[i]][k]);
            if(val > max_val) {
                max_val = val;
                max_index = i;
            }
        }
        if(max_val == 0) {
            throw runtime_error("Matrix is singular.");
        }
        swap(index[k], index[max_index]);
        for(int i = k+1; i < n; i++) {
            double factor = A[index[i]][k] / A[index[k]][k];
            A[index[i]][k] = factor;
            for(int j = k+1; j < n; j++) {
                A[index[i]][j] -= factor * A[index[k]][j];
            }
            b[index[i]] -= factor * b[index[k]];
        }
    }

    // Back substitution
    vector<double> x(n);
    for(int i = n-1; i >= 0; i--) {
        double sum = 0;
        for(int j = i+1; j < n; j++) {
            sum += A[index[i]][j] * x[j];
        }
        x[i] = (b[index[i]] - sum) / A[index[i]][i];
    }

    return x;
}





int main()
{
    // Option 1
    //vector<double> x=gauss_elimination(A, B);
   
    
    // Option 2
    /*
    vector<double> x(B.size());
    
    vector<double> oneDVec(64); // initialize 1D vector with size 64

    // copy elements from 2D vector to 1D vector
    int index = 0;
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            oneDVec[index] = A[i][j];
            index++;
        }
    }

    GaussianElimination(oneDVec, B,x);
    */
}