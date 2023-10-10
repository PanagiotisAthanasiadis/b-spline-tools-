#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <set>
#include <algorithm>

using namespace std;


bool isEqual(double a, double b) {
    const double epsilon=1e-10;
    return fabs(a - b) < epsilon;
}


double gdiv(double a,double b)
{
    if(a==0.0&&b==0.0)
    {
        return(0.0);
    }    
    else
    { 
        return(a/b);
    }
}


// Function to calculate the B-spline value
/*
    i is the index of the B-spline basis function in the B-spline basis set. 
    This is a zero-based index, so the first B-spline basis function in the set has index 0.

    ord is the order of the B-spline basis function. 
    This is also known as the degree of the B-spline basis function. 
    The order is an integer greater than or equal to 1 that determines the number of 
    knots used to define the B-spline basis function.

    x is the value at which the B-spline basis function should be evaluated.

    nk is the number of knots in the knot vector kns.

    kns is a vector of knots that define the B-spline basis function. 
    The knots are sorted in ascending order, 
    and the B-spline basis function is defined over the range between the first and last knot. 
    The size of kns should be equal to nk.

*/
double bsp(int i, int ord, double x, int nk, vector<double> kns){
  // Check for illegal value of i
  if(i<0 || i>nk-ord-1){
    cout<<"illegal i value: i="<<i<<"; nk-ord="<<nk<<"-"<<ord<<"="<<nk-ord<<endl;
    return numeric_limits<double>::quiet_NaN();
  }

  // Return 0 if x is outside the interval defined by the knots
  if(x<kns[i] || x>kns[i+ord]) return 0.0;

  // Remove repeated knots
  int k=nk-1; 
  while(kns[k]==kns[k-1]) k--; 
  k--; 

  // If ord is 1, return 1 if x is within the interval defined by the knots
  if(ord == 1){
    if(i != k){
      return((kns[i]<=x && x<kns[i+1]) ? 1.0 : 0.0);
    }else{
      if(i == k){
        return((kns[i]<=x && x<=kns[i+1]) ? 1.0 : 0.0);
      }else return numeric_limits<double>::quiet_NaN();;
    }
  }
  // If ord is greater than 1, return the sum of two recursive calls
  else{
    return(gdiv((x-kns[i])*bsp(i, ord-1, x, nk, kns), kns[i+ord-1]-kns[i]) +
           gdiv((kns[i+ord]-x)*bsp(i+1, ord-1, x, nk, kns), kns[i+ord]-kns[i+1])
    );
  }
}

/*
This code is an implementation of the Gauss-Legendre quadrature method. 

The function has several parameters:



    n: This is the number of sample points that will be used to approximate the definite integral. 
    The function supports up to 4 sample points, but more could be added to the switch statement.

    a: This is the lower bound of the interval of integration.

    b: This is the upper bound of the interval of integration.

    order: The order of the B-spline curve.

    nk: The number of knots in the B-spline curve.

    kns: A vector of knots in the B-spline curve.

    f: This is a pointer to the bsp function, 
    which will be used to evaluate the B-spline at each of the sample points.

The function sets up arrays x and w to store the abscissas and weights, respectively.
The abscissas and weights are computed based on the number of sample points n specified by the user. 
Then, a loop is used to evaluate the B-spline at each of the sample points using the bsp function, 
and the results are summed and weighted by the weights. 
The final result is the approximation of the definite integral.

*/


// Gauss-Legendre quadrature
double gauss_legendre(int n, double a, double b,int order,int nk,vector<double> kns ,double (f)(int,int,double,int,vector<double>)) {
double *x = new double[n];
double* w = new double[n];

switch(n) {
case 1:
x[0] = 0;
w[0] = 2;
break;
case 2:
x[0] = -0.5773502691896257;
x[1] = 0.5773502691896257;
w[0] = 1;
w[1] = 1;
break;
case 3:
x[0] = -0.7745966692414834;
x[1] = 0;
x[2] = 0.7745966692414834;
w[0] = 0.5555555555555556;
w[1] = 0.8888888888888888;
w[2] = 0.5555555555555556;
break;
case 4:
x[0]=-0.8611363115940526;
x[1]=-0.3399810435848563;
x[2]=0.3399810435848563;
x[3]=0.8611363115940526;
w[0]=0.34785484513745385;
w[1]=0.6521451548625461;
w[2]=0.6521451548625461;
w[3]=0.34785484513745385;
break;
default:
cout << "Error";
break;
}

double local_sum = 0; // Sum of a B-spline basis index (ex B1 at [0,1]) in a specific interval
double bspline_sum=0; // Sum of all B-splines basis index(ex B1 at [0,1],[0,2]...) in a all intervals 
double total_sum=0; //  Sum of all B-splines basis of every index in all intervals

vector<double> v_bspline_sum; // Vector of sum of all B-splines basis index(ex B1 at [0,1],[0,2]...) in a all intervals 


set<pair<double, double>> intervals; // Set of pair of intervals(lower bound,upper bound)
vector<vector<double>> xq;
vector<double> coefficient_intervals;

// Single time claculation for intervals
for (int i = 0; i < kns.size() - 1; i++) 
{
  double start = kns[i];
  double end = kns[i+1];

  if (isEqual(start, end)) 
  {
    continue; // Skip over intervals that start and end with the same value
  }

  intervals.insert(make_pair(start, end));
}

// Single time claculation for quadrature points and interval coefficients and insert into 2d vector(intervals,xi)

for (auto interval : intervals) {
  //cout << interval.first << "-" << interval.second << endl;
  vector<double> xi;
  for (int i = 0; i < n; i++) 
  {
    xi.push_back((interval.second-interval.first)/2 * x[i] + (interval.second+interval.first)/2);
  }
  xq.push_back(xi);
  coefficient_intervals.push_back((interval.second-interval.first)/2);
}




for(int j=0;j<nk-order; j++)
{
  bspline_sum=0; // Reset the bpline sum 
  for(int h= 0; h < xq.size(); h++)
  {
    local_sum = 0; //  Reset the local sum
    for (int i = 0; i < xq[h].size() && i < n; ++i)
    {
      local_sum += w[i]*f(j, order, xq[h][i], kns.size(), kns);
      i++;
    }
    local_sum *= coefficient_intervals[h];
    bspline_sum += local_sum; 
  } 
  v_bspline_sum.push_back(bspline_sum);
  total_sum+=bspline_sum;
}

for (int h= 0; h < v_bspline_sum.size(); h++)
{
  cout << "B("<< h+1 <<","<< order-1 << ")" << ":"<<v_bspline_sum[h] << endl;
}



delete[] x; //free memory
delete[] w; //free memory
return total_sum;
}


int main()
{
    // Define the knots vector
    vector<double> kns = {0.0, 0.0, 0.0, 1.0, 2.0, 3.0,4.0,4.0, 5.0,5.0,5.0 };
    //vector<double> kns = {0.0,  1.0, 2.0, 3.0,4.0, 5.0};
  
    // Define the order of the B-spline
    int order = 3;
    double result = gauss_legendre(4, 0, 5,order,kns.size(),kns,&bsp);
    cout << "Result: " << result;

    
    return 0;
}