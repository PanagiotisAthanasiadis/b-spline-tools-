#include <iostream>
#include <vector>
#include <cmath>
#include <limits>


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


// Function to calculate the B-spline value
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





// Function to calculate the values of B-spline bases at given points
vector<double> bsbasesCpp(vector<double> xs, vector<double> kns, int order){
  int nx = xs.size(); // number of x values
  int nk = kns.size(); // number of knots
  int ord = order; // order of the B-spline
  int nb = nk - ord; // number of B-spline bases
  int ansSize = nb * nx; // size of the output vector
  vector<double> ans(ansSize); // output vector
  
  // Loop through all x values and calculate the value of each B-spline basis at each x
  for(int i = 0; i < nx; i++){
    for(int j = 0; j < nb; j++){
      ans[j + i * nb] = bsp(j, ord, xs[i], nk, kns);
    }
  }
  
  // Return the output vector
  return(ans);
}


vector<double> kns = {0, 0, 0, 1, 2, 3, 3, 3};
int ord = 2;

int main()
{
   vector<double> coef(kns.size(),1);

  vector<double> r=bsbasesCpp(kns, kns,ord);

  return 0;
}