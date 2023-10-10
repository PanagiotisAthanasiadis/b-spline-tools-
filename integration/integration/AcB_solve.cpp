#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <set>
#include <algorithm>
#include <fstream>

using namespace std;


// Function which we will aproximate 
double f(double t)
{
    return sin(2*M_PI * t);
}

double f1(double t) //  This aproximates well 
{
    return sqrt(t+4);
}

double f2(double t) //  This aproximates well 
{
    return t*t;
}

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


vector<double> create_Pn(double a,double b,double dt)
{
    vector<double> pn;
    pn.reserve((b-a)/dt + 1);
    //cout << "DEBUG" << (b-a)/dt + 1 ;
    while(a < b || isEqual(a,b) )
    {
        pn.push_back(a);
        a += dt;

    }
    /*for (auto i = pn.begin(); i != pn.end(); ++i)
        cout << *i << " "; */
    return pn;
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

vector<vector<double>> create_a_gauss_legendre(int n, double a, double b,int order,int nk,vector<double> kns ,double (f)(int,int,double,int,vector<double>)) {
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
vector<vector<double>> A;


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



for(int u=0;u<nk-order; u++)
{
  for(int j=0;j<nk-order; j++)
  {
    bspline_sum=0; // Reset the bpline sum 
    for(int h= 0; h < xq.size(); h++)
    {
      local_sum = 0; //  Reset the local sum
      for (int i = 0; i < xq[h].size() && i < n; ++i)
      {
        local_sum += w[i]*f(u, order, xq[h][i], kns.size(), kns) *
        f(j, order, xq[h][i], kns.size(), kns);
        //i++;
      }
      local_sum *= coefficient_intervals[h];
      bspline_sum += local_sum; 
    } 
    v_bspline_sum.push_back(bspline_sum);
  }
  A.push_back(v_bspline_sum); // Some error here 
  v_bspline_sum.clear(); // Empty the vector for next iteration 
}

delete[] x; //free memory
delete[] w; //free memory
return A;
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
      //i++;
    }
    local_sum *= coefficient_intervals[h];
    bspline_sum += local_sum; 
  } 
  v_bspline_sum.push_back(bspline_sum);
  total_sum+=bspline_sum;
}






delete[] x; //free memory
delete[] w; //free memory
return total_sum;
}


vector<double> create_b_gauss_legendre(int n, double a, double b,int order,int nk,vector<double> kns ,double (f)(int,int,double,int,vector<double>),double (g)(double)) {
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
      local_sum += w[i]*f(j, order, xq[h][i], kns.size(), kns)*g(xq[h][i]);
      //i++;
    }
    local_sum *= coefficient_intervals[h];
    bspline_sum += local_sum; 
  } 
  v_bspline_sum.push_back(bspline_sum);
  
}





delete[] x; //free memory
delete[] w; //free memory
return v_bspline_sum;
}


vector<double> create_Bf(vector<double> coefficients,vector<double> xi,int order,vector<double> kns)
{
  vector<double> results;
  double bf=0; // Temp variable 
  for(auto x : xi)
  {
    for(int i=0; i<coefficients.size(); i++)
    {
      //Error
      bf+= coefficients[i] * bsp(i, order, x, kns.size(), kns); // Calculate the aproximations
      //bf+= coefficients[i] * gaus_legendre(4,);
    }
    results.push_back(bf); // Add to the results vector  
    bf=0; // Reset 
  }
  return results;
}


double custom_create_Bf(vector<double> coefficients,double xi,int order,vector<double> kns)
{
  double bf=0; // Temp variable 
  for(int i=0; i<coefficients.size(); i++)
  {
  
    bf+= coefficients[i] * bsp(i, order, xi, kns.size(), kns); // Calculate the aproximations
  }
  return bf;
}


double error_calc(vector<double> coefficients,int n, double a, double b,int order,int nk,vector<double> kns ,double (bf)(vector<double>,double,int,vector<double>),double (f)(double)) {
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


  for(int h= 0; h < xq.size(); h++)
  {
    local_sum = 0; //  Reset the local sum
    for (int i = 0; i < xq[h].size(); ++i)
    {
      //local_sum += w[i]*pow(pow((bf(coefficients,xq[h][i],order,kns)-f(xq[h][i])),2),0.5);
      local_sum += w[i]*pow(pow((bf(coefficients,xq[h][i],order,kns)-f(xq[h][i])),2),0.5);
      //local_sum += w[i]*(bf(coefficients,xq[h][i],order,kns)-f(xq[h][i]));
    }
    local_sum *= coefficient_intervals[h];
    bspline_sum += local_sum; 
  } 


delete[] x; //free memory
delete[] w; //free memory
return bspline_sum;
}

// Utility function 
void print2DVectorToFile(const vector<vector<double>>& vec, ofstream& outputFile) {
    if (!outputFile.is_open()) {
        cerr << "Error: output file not open" << endl;
        return;
    }

    const int numRows = vec.size();
    const int numCols = vec[0].size();

    // Print column headers
    outputFile << setw(12) << " ";
    for (int col = 0; col < numCols; ++col) {
        outputFile << setw(12) << col;
    }
    outputFile << endl;

    // Print row headers and vector contents
    for (int row = 0; row < numRows; ++row) {
        outputFile << setw(12) << row;
        for (int col = 0; col < numCols; ++col) {
            outputFile << setw(12) << fixed << setprecision(6) << vec[row][col] << "(" << row << "," << col << ")";
        }
        outputFile << endl;
    }
}



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

int main()
{
    // Define the lower and upper bound 
    int a=0;
    int b=3;

    // Define the order of the B-spline
    int order = 3;

    // Precision vector 
    vector<double> dt_v ={0.2,0.1,0.05,0.025};

    vector<double>error_v;

    for(auto dt:dt_v)
    {
      // Define the knots vector and make it repeat the values at the start and the end based on the order
      vector<double>kns=create_Pn(a, b, dt);
      for(int i=0; i<order-1; i++)
      {
        kns.insert(kns.begin(),a);
        kns.push_back(b);
      }
      
      //vector<double> kns = {0.0, 0.0, 0.0,0.25, 0.5,0.6,0.75,1.0,1.1,1.25,1.5,1.75,2,2,2 };
      //vector<double> kns = {0.0,  1.0, 2.0, 3.0,4.0, 5.0};
      
      // Create the right hand side
      vector<double> B = create_b_gauss_legendre(4, a, b,order,kns.size(),kns,&bsp,&f);
    
      // Create the left hand side 
      vector<vector<double>> A = create_a_gauss_legendre(4,a,b,order,kns.size(),kns,&bsp);
      
      // Solve the system and store the results
      vector<double> x = gauss_elim(A, B);

      // Create the partition vector
      vector<double> xi = create_Pn(a, b, 0.01);

      // Aproximations vector

      vector<double> Bf = create_Bf(x, xi, order,  kns);

      double temp22=error_calc(x, 4,  a,  b,  order,  kns.size(),  kns,&custom_create_Bf,&f);
      cout <<"Real Error:" <<  temp22 << " Dt:" << dt << endl;
      error_v.push_back(temp22);
    }
    for(int i=0; i<error_v.size(); i++)
    {
      cout << "r:" <<log(error_v[i+1]/error_v[i]) / log(dt_v[i+1]/dt_v[i]) << endl;
    }
    return 0;
}