#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Gaussian elimination method to solve Ax=b
vector<double> gauss_elim(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    // Augment the coefficient matrix with the right-hand side vector
    for (int i = 0; i < n; i++) {
        A[i].push_back(b[i]);
    }

    // Perform Gaussian elimination
    for (int i = 0; i < n; i++) {
        // Find pivot row and swap if necessary
        int max_row = i;
        for (int j = i+1; j < n; j++) {
            if (abs(A[j][i]) > abs(A[max_row][i])) {
                max_row = j;
            }
        }
        swap(A[i], A[max_row]);

        // Perform row operations to eliminate coefficients
        for (int j = i+1; j < n; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < n+1; k++) {
                A[j][k] -= factor * A[i][k];
            }
        }
    }

    // Back-substitution to solve for x
    vector<double> x(n, 0.0);
    for (int i = n-1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i+1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (A[i][n] - sum) / A[i][i];
    }

    return x;
}

int main() {
    // Example usage
    vector<vector<double>> A = {{3, 1, -1}, {2, 4, 2}, {-1, 2, 5}};
    vector<double> b = {4, 1, 1};

    vector<double> x = gauss_elim(A, b);

    cout << "Solution:" << endl;
    for (int i = 0; i < x.size(); i++) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

    return 0;
}
