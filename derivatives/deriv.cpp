#include <iostream>
#include <system_error>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <set>

using namespace std;


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

bool isEqual(double a, double b) {
    const double epsilon=1e-10;
    return fabs(a - b) < epsilon;
}





double bsp(int i, int ord, double x, int nk, vector<double> kns){
  if(i<0||i>nk-ord-1){
    cout<<"illegal i value: i="<<i<<"; nk-ord="<<nk<<"-"<<ord<<"="<<nk-ord<<endl;
    return numeric_limits<double>::quiet_NaN();
  }
  if(x<kns[i]||x>kns[i+ord])return(0.0);
  int k=nk-1; while(kns[k]==kns[k-1])k--; k--; 
  if(ord==1){
    if(i!=k){
      return((kns[i]<=x && x<kns[i+1])? 1.0 : 0.0);
    }else{
      if(i==k){
	return((kns[i]<=x && x<=kns[i+1])? 1.0 : 0.0);
      }else return numeric_limits<double>::quiet_NaN();;
    }
  }else
    return(gdiv((x-kns[i])*bsp(i,ord-1,x,nk,kns), kns[i+ord-1]-kns[i])+
	   gdiv((kns[i+ord]-x)*bsp(i+1,ord-1,x,nk,kns),kns[i+ord]-kns[i+1])
	   );
}

double deriv(int k,vector<double> kns,int i,int t)
{
  return (k/(kns[i+k]-kns[i])) * (bsp(i, k-1,t, kns.size(), kns)-bsp(i+1, k-1,t, kns.size(), kns));
}


vector<double> product_rule_gauss_legendre(int n, double a, double b,int order,int nk,vector<double> kns ,double (f)(int,int,double,int,vector<double>)) {
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


for(int u=0;u<nk-order; u++)
{
    bspline_sum=0; // Reset the bpline sum 
    for(int h= 0; h < xq.size(); h++)
    {
      local_sum = 0; //  Reset the local sum
      for (int i = 0; i < xq[h].size() && i < n; ++i)
      {
        //local_sum += w[i]*f(u, order-1, xq[h][i], kns.size(), kns); // The derivative of the b-spline
        local_sum += w[i]*deriv(order,kns,u,xq[h][i]);
        for(int o=0; o<nk-order; o++)
        {
          if(u == o)
          {
            continue; //  Avoid calculating the original of the derivative B-spline 
          }
          local_sum *= f(o, order, xq[h][i], kns.size(), kns);
        }
        local_sum *= coefficient_intervals[h];
      }
    bspline_sum += local_sum;  
  }
  v_bspline_sum.push_back(bspline_sum);
}

double temp=0;
for(auto x:v_bspline_sum)
{
  cout << x<<endl;
  temp+=x;
}


cout << " Result:" << temp; 

delete[] x; //free memory
delete[] w; //free memory
return v_bspline_sum;
}






int main()
{
    int a=0;
    int b=5;
    int order = 3;
    vector<double> kns = {0, 0, 0, 1, 2, 3,4,4,5,5,5};

    vector<double> x=product_rule_gauss_legendre(4, a,  b,  order, kns.size(),  kns, &bsp);
    
    return 0;
}
