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


vector<double> bsbasesCpp(vector<double> xs, vector<double> kns, int order){
  int nx=xs.size(), nk=kns.size(),i,j;
  int ord=order, nb=nk-ord;
  int ansSize=nb*nx;
  vector<double> ans(ansSize);
  for(i=0;i<nx;i++){
    for(j=0;j<nb;j++){
      ans[j+i*nb]=bsp(j,ord,xs[i],nk,kns);
    }
  }
  return(ans);
}

vector<double> bsplineCpp(vector<double> xs, int ord, vector<double> kns, vector<double> coef){
  int nx=xs.size(), nk=kns.size(), i, j;
  vector<double> ans(nx,0.0);
  for(i=0;i<nx;i++)
    for(j=0; j<nk-ord; j++)
      ans[i] += coef[j]*bsp(j,ord,xs[i],nk,kns);
  return(ans);
}


vector<double> ibsCpp(vector<double> xs, int ord, vector<double> kns, vector<double> coef){
  int nx=xs.size(), nk=kns.size(), nc=coef.size(),i, j, k;
  double tmp;
  vector<double> ans(nx,0.0),kns_ext(nk+1),bet(nc);
  for(i=0;i<nk;i++){
    kns_ext[i]=kns[i];
  }
  kns_ext[nk]=kns_ext[nk-1];
  bet[0]=coef[0]*(kns_ext[ord]-kns[0])/ord;
  for(i=1;i<nc;i++){
    bet[i]=bet[i-1]+coef[i]*(kns_ext[i+ord]-kns_ext[i])/ord;
  }
  ans=bsplineCpp(xs,ord+1,kns_ext,bet);
  return(ans);
}





int main()
{
    int ord = 2;
    vector<double> kns = {0, 0, 0, 1, 2, 3, 3, 3};
    vector<double> coef(kns.size(),1);

    vector<double> r=ibsCpp(kns, ord,kns, coef);
    cout << "size of kns: " << kns.size() << " " << "size of r " << r.size() << endl ;
    
    auto x=kns.begin();
    auto y=r.begin();

    while( x != kns.end() && y!= r.end())
    {
        cout << "point of integration:" << *x << " " << "result:" << *y << " " <<  endl;
        x++;
        y++;
    }
    return 0;
}