#include <vector>
#include <cmath>


using namespace std;
/*
This implementation assumes that the knot vector is non-decreasing and has a length of n + degree + 1, where n is the number of control points. The knot vector is used to specify the parametric domain of the curve and determines how the control points influence the shape of the curve.
To use this function, you would need to provide the degree of the B-spline curve (e.g. 3 for cubic B-splines), the knot vector, and the control points as input. The control points are a 2D array of size n x d, where n is the number of control points and d is the dimension of the curve (e.g. 2 for a 2D curve). The function will return the coordinates of the point on the curve corresponding to the parameter value t.
*/

// Compute the B-spline basis functions for a given knot vector and degree
vector<double> bspline_basis(int degree, const vector<double>& knots, double t) {
  int n = knots.size() - degree - 1;
  vector<double> basis(n);
  vector<double> left(degree + 1);
  vector<double> right(degree + 1);
  basis[0] = 1.0;
  for (int j = 1; j <= degree; j++) {
    left[j] = t - knots[n - j];
    right[j] = knots[n + j] - t;
    double saved = 0.0;
    for (int r = 0; r < j; r++) {
      double temp = basis[r] / (right[r + 1] + left[j - r]);
      basis[r] = saved + right[r + 1] * temp;
      saved = left[j - r] * temp;
    }
    basis[j] = saved;
  }
  return basis;
}

// Evaluate a B-spline curve at a given parameter value
vector<double> bspline_eval(int degree, const vector<double>& knots, const vector<vector<double>>& control_points, double t) {
  vector<double> basis = bspline_basis(degree, knots, t);
  vector<double> point(control_points[0].size(), 0.0);
  for (int i = 0; i < control_points.size(); i++) {
    for (int j = 0; j < control_points[i].size(); j++) {
      point[j] += basis[i] * control_points[i][j];
    }
  }
  return point;
}
