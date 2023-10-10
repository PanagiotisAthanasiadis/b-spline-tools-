#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <tuple>
#include <fstream>
#include <iomanip> 
#include <string>


using namespace std;

// Helper functions for bold text and colors

ostream& bold_on(ostream& os)
{
    return os << "\e[1m";
}


ostream& bold_off(ostream& os)
{
    return os << "\e[0m";
}

ostream& red(ostream& os)
{
    return os << "\e[31m";
}

ostream& black(ostream& os)
{
    return os << "\e[30m";
}

ostream& blue(ostream& os)
{
    return os << "\e[34m";
}

ostream& yellow(ostream& os)
{
    return os << "\e[33m";
}

ostream& green(ostream& os)
{
    return os << "\e[32m";
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
tuple< double,vector<double>,vector<vector<double>>,vector<vector<double>> > gauss_legendre(int n, double a, double b,int order,int nk,vector<double> kns ,double (f)(int,int,double,int,vector<double>)) {
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

double sum = 0;


// Single time claculation for quadrature points
vector<double> xi;
for (int i = 0; i < n; i++) 
{
  xi.push_back((b-a)/2 * x[i] + (b+a)/2);
}

//cout << bold_on << "a: " << a << " b: " << b << endl <<bold_off; 

vector<vector<double>> Bxq;
vector<vector<double>> WBxq;


for(int j=0; j<nk-order; j++)
{
    // Create a new row vector
    vector<double> row;
    vector<double> row2;

    for (int i = 0; i < n; i++) 
    {
        row.push_back(f(j, order, xi[i], kns.size(), kns));
        //cout << red << "xi:"<< xi[i] << " n:"<< j << " : "<<row.back() << "\n";
        
        row2.push_back(w[i]*row.back());
        sum += row2.back(); // Avoiding calculating 2 times 
        //cout << blue << "w*f: " <<  row2.back() << endl;

       // cout<<yellow << "sum:"<< sum << endl;
    }
    Bxq.push_back(row);
    WBxq.push_back(row2);
}

sum *= (b-a)/2;
delete[] x; //free memory
delete[] w; //free memory
return {sum,xi,Bxq,WBxq};
}


int main()
{
    // Define the knots vector
    vector<double> kns = {0.0, 0.0, 0.0, 1.0, 2.0, 3.0,4.0, 5.0,5.0,5.0 };
    //vector<double> kns = {0.0,  1.0, 2.0, 3.0,4.0, 5.0};
  
    // Define the order of the B-spline
    int order = 3;


    // Define the number of quadrature points to use 
    int n = 4;

    // Total values 
    double total_result=0;

    // Create an empty 2D vector to store the quadrature points
    vector<vector<double>> xq;


     // Create an output file stream object
    ofstream outfile;

    // Open a file for writing
    outfile.open("Results.txt", ios::out);

    // Check if the file opened successfully
    if (!outfile.is_open()) {
      // Handle the error
      cout << "Error opening the file";
      return 1;
    }


    // Create a vector of intervals
    vector<string> intervals;

    // Create the intervals
    for(int i=0; i<kns.size()-1; i++) // kns.size()-1 because iterator overflows (i+1)
    {
      if(kns[i]==kns[i+1]) // Checking if the current knot and the next one is the same 
      {
        continue;
      } 
      intervals.push_back("["+to_string((int(kns[i]))) + "," + to_string(int(kns[i+1])) +"]" );
    }

    // Create a vector of labels
    vector<string> labels;

    // Create the labels
    for(int j=0; j<kns.size()-order; j++)
    {
     
      labels.push_back("B(" + to_string(j+1) + "," + to_string(order-1)+")"); // order-1 and j+1 because the index is zero based
    }


    // Create a vector of labels
    vector<string> labels_quadrature;

    // Create the labels
    for(int j=0; j<n; j++)
    {
     
      labels_quadrature.push_back("Xq"+ to_string(j+1));
    }

    // Create a 3D vector to store the B(Xqi)
    vector<vector<vector<double>>> Bxq;

    // Create a 3D vector to store the Wi*B(Xqi)
    vector<vector<vector<double>>> WBxq;

    // Create a vector to keep the sum of every interval
    vector<double> sum;

    // Code for finding every unique interval 
    for(int i=0; i<kns.size()-1; i++) // kns.size()-1 because iterator overflows (i+1)
    {
      if(kns[i]==kns[i+1]) // Checking if the current knot and the next one is the same 
      {
        continue;
      } 
      
      auto [temp,temp2,temp3,temp4] = gauss_legendre(n, kns[i], kns[i+1],order,kns.size(),kns,&bsp);


      sum.push_back(temp);
      xq.push_back(temp2); // Push back the quadrature points for every knot interval  
      Bxq.push_back(temp3); 
      WBxq.push_back(temp4);
      
      total_result+=temp;
      
    }

    cout  << "\e[31m Total_Result: " << total_result;


    outfile <<"Quadrature Points"<< endl;
    // Print the contents of the data vector into the file
    int i=0; // intervals iterator 

    outfile << "\t\t ";
      for (const string& s : labels_quadrature) 
      {
        outfile << s << "   "; 
      }
       outfile << endl;

    for (const auto& row : xq) {
      outfile << intervals[i] << " " ;
      i++;
      for (const auto& value : row) 
      {
        outfile <<fixed << setprecision(4) << value << " ";
      }
      outfile << endl;
    }

    
    outfile << endl;
    outfile <<"Bi,j(Xq)"<< endl <<endl;

    // Create label iterator
    auto it = begin(labels);


    // Reset intervals iterator
    i=0;
   
    // Print the contents of the data vector into the file
    for (const auto& row : Bxq) {
      outfile << intervals[i] << " " << endl ;
      i++;
      it=begin(labels); // Reset label iterator

      // Iterate over the vector of labels_quadrature and prints them to the file
     
      outfile << "\t\t ";
      for (const string& s : labels_quadrature) 
      {
        outfile << s << "   "; 
      }
       outfile << endl;
      for (const auto& value : row) 
      {
        outfile << *it << " ";
        it++;
        for (const auto& value1 : value) 
        {
          outfile <<fixed << setprecision(4) << value1 << " " ;
        }
        outfile << endl;
      }
      outfile << endl;
    }


    outfile << endl;
    outfile <<"Wi * Bi,j(Xq)"<< endl <<endl;
    // Reset intervals iterator
    i=0;

    for (const auto& row : WBxq) {
      outfile << intervals[i] << " " << endl ;
      i++;
      it=begin(labels); // Reset label iterator

      // Iterate over the vector of labels_quadrature and prints them to the file
     
      outfile << "\t\t ";
      for (const string& s : labels_quadrature) 
      {
        outfile << s << "   "; 
      }
       outfile << endl;
      for (const auto& value : row) 
      {
        outfile << *it << " ";
        it++;
        for (const auto& value1 : value) 
        {
          outfile <<fixed << setprecision(4) << value1 << " " ;
        }
        outfile << endl;
      }
      outfile << endl;
    }

    outfile << endl;
    outfile <<"Sum for every interval:"<< endl <<endl;
    // Reset intervals iterator
    i=0;
    for (const auto& x : sum) {
      outfile << intervals[i] << ": " << x << endl ;
      i++;
    }

    // Close the file
    outfile.close();
    
    
    return 0;
}