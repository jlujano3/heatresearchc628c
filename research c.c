/* solid_problem2.c follows the theory and method of laser_newsystem.c to solve the solidification problem.
Includes:  New problem with analytic solution, New system from subtracting u and v (u with nontrivial IC) to avoid the newton potential, 
simple singularity subtraction, forward euler for r, and error approximation. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ONETHRD 0.333333333333333333333333333333
#define M_SQPI 1.77245385090552
#define M_SQPINV 0.564189583547756287
#ifndef MIN
#define MIN(A,B)  ( (A) > (B) ? (B) : (A) )
#endif

double rKnown(double t){
  double g;
  g=t*t*t;
  return g;
}

double dRKnown(double t){
  double g;
  g=3*t*t;
  return g;
}

/*qKnownA is for r(t)=t*/

double qKnownA(double xo,double t,double to){
  double g;
  g=exp(-(rKnown(t)-xo)*(rKnown(t)-xo)/(4.0*(t+to)))*(2.0*t-rKnown(t)+xo+2.0*to)/(2.0*(t+to)*sqrt(t+to));
  return g;
}

/*qKnownB is for r(t)=1*/

double qKnownB(double xo,double t,double to){
  double g;
  g=exp(-(rKnown(t)-xo)*(rKnown(t)-xo)/(4*(t+to)))*(-rKnown(t)+xo)/(2*(t+to)*sqrt(t+to));
  return g;
}

/*qKnownB is for r(t)=t^2*/

double qKnownC(double xo,double t,double to){
  double g;
  g=exp(-(rKnown(t)-xo)*(rKnown(t)-xo)/(4*(t+to)))*(-rKnown(t)+xo+4*t*(t+to))/(2*(t+to)*sqrt(t+to));
  return g;
}

double qKnownD(double xo,double t,double to){
  double g;
  g=exp(-(rKnown(t)-xo)*(rKnown(t)-xo)/(4*(t+to)))*(-rKnown(t)+xo+6*t*t*(t+to))/(2*(t+to)*sqrt(t+to));
  return g;
}

double kv(double x,double xo,double t,double to){
  double g;
  g=1.0/sqrt(4.0*M_PI)*exp(-(x-xo)*(x-xo)/(4.0*(t-to)));
  return g;
}

/*Used to test the single layer potential with trivial kv*/

double kv2(double x,double xo,double t,double to){
  double g;
  g=1.0;
  return g;
}

double kvtntn(){
  double g;
  g=1.0/sqrt(4.0*M_PI);
  return g;
}

/*Used to test the single layer potential with trivial kv*/

double kvtntn2(){
  double g;
  g=1.0;
  return g;
}

double kd(double x,double xo,double t,double to){
  double g;
  g=exp(-(x-xo)*(x-xo)/(4.0*(t-to)))*(x-xo)/(4.0*sqrt(M_PI)*(t-to));
  return g;
}

double kdtntn(double t){
  double g;
  g=dRKnown(t)/(4.0*sqrt(M_PI));
  return g;
}

double u(double x,double xo,double t,double to){
  double g;
  g=exp(-(x-xo)*(x-xo)/(4.0*(t+to)))/sqrt(t+to);
  return g;
}

double uo(double x,double xo,double to){
  double g;
  g=u(x,xo,0.0,to);
  return g;
}

double sLay(int n,double h,double *r,double *q){
  double tn=n*h;
  int i;
  double sum=0.5*kv(r[n],r[0],tn,0)*q[0]/sqrt(tn-0);
  for (i=1; i<n; i++){
    sum += kv(r[n],r[i],tn,i*h)*q[i]/sqrt(tn-i*h);
  }
  return sum*h;
}

double sLay2(int n,double h,double *r,double *q){
  double tn=n*h;
  int i;
  double sum=0.5*kv2(r[n],r[0],tn,0)*q[0]/sqrt(tn-0);
  for (i=1; i<n; i++){
    sum += kv2(r[n],r[i],tn,i*h)*q[i]/sqrt(tn-i*h);
  }
  return sum*h;
}

double dLay(int n,double h,double *r,double *u){
  double tn=n*h;
  int i;
  double sum=0.5*kd(r[n],r[0],tn,0)*u[0]/sqrt(tn-0);

  for (i=1; i<n; i++){
    sum += kd(r[n],r[i],tn,i*h)*u[i]/sqrt(tn-i*h);
  }
  return sum*h;
}

double mun(int n,double h){
  double tn=n*h;
  int i;
  double sum=0.5*h/sqrt(tn-0);
  for (i=1; i<n; i++){
    sum += h/sqrt(tn-i*h);
  }
  return 2.0*sqrt(tn)-sum;
}

double initPot(int n, double h, double xo,double to,double *r){
  int numInts=10000;
  double b=10.0, z;
  double tn=(n)*h, s4t=sqrt(4.0*tn);
  double a=r[n]/s4t;

  if (a<b){
    a=a;
  }
  else{
    a=b;
  }

  double delta=(b+a)/numInts;
  double sum=0.5*(exp(-a*a)*uo(r[n]-a*s4t,xo,to)+exp(-b*b)*uo(r[n]+b*s4t,xo,to));
  int i;
  for(i=1; i<numInts; i++){
    z = -a+delta*i;
    sum += exp(-z*z)*uo(r[n]+z*s4t,xo,to);
  }
  return sum*delta*M_SQPINV;
}

double qTestExact(double t){
  double g;
  g=exp(t)*sqrt(M_PI)*erf(sqrt(t));
  return g;
}

int main(int nargs, char *argv[]){
  char title1[20] = "kv";
  char title2[20] = "kd";
  char title3[20] = "u";
  char title4[20] = "uo";
  char title5[20] = "rKnown";
  char title6[20] = "dRKnown";
  char title7[20] = "qKnownA";
  char title8[20] = "qKnownB";
  char title9[20] = "sLay";
  char title10[20] = "dLay";
  char title11[20] = "mun";
  char title12[20] = "InitPot";
  char title13[20] = "Difference";
  char title14[20] = "LHS";
  char title15[20] = "qKnown";
  char title16[20] = "exact";
  char title17[20] = "t";
  char title18[20] = "munValsD";
  char title19[20] = "IP Test";
  char title20[20] = "Analytical";
  char title21[20] = "Numerical";
  char title22[20] = "sLay Test";

  double tmax=1.0;
  int numInts=1000;
  double dt=tmax/numInts;
  double xo=2.0;
  double to=2.0;
  int i;
  double *r;
  double *t;
  double *qKnownList;
  double *q;
  double *qTest;
  double *uPasser;
  double *diff;
  double *sLays;
  double *dLays;
  double *munValsD;
  double *initPote;
  double *denom;
  
  r=(double*)malloc( (numInts+1)*sizeof(double) );
  t=(double*)malloc( (numInts+1)*sizeof(double) );
  qKnownList=(double*)malloc( (numInts+1)*sizeof(double) );
  uPasser=(double*)malloc( (numInts+1)*sizeof(double) );
  qTest=(double*)malloc( (numInts+1)*sizeof(double) );
  for(i=0;i<numInts+1;i++){
    t[i]=0+i*dt;
    r[i]=rKnown(t[i]);
    qKnownList[i]=qKnownD(xo,t[i],to);
    uPasser[i]=u(r[i],xo,0+i*dt,to);
    qTest[i]=exp(t[i]);
  }
  q=(double*)malloc( (numInts+1)*sizeof(double) );
  q[0]=qKnownD(xo,t[0],to);
  diff=(double*)malloc( (numInts+1)*sizeof(double) );
  diff[0]=0;
  sLays=(double*)malloc( (numInts+1)*sizeof(double) );
  sLays[0]=0;
  dLays=(double*)malloc( (numInts+1)*sizeof(double) );
  dLays[0]=0;
  munValsD=(double*)malloc( (numInts+1)*sizeof(double) );
  munValsD[0]=0;
  initPote=(double*)malloc( (numInts+1)*sizeof(double) );
  initPote[0]=0;
  denom=(double*)malloc( (numInts+1)*sizeof(double) );
  denom[0]=0;
  
  int j;
  for(j=1;j<numInts+1;j++){
    sLays[j]=sLay(j,dt,r,q);
    dLays[j]=dLay(j,dt,r,uPasser);
    munValsD[j]=kdtntn(j*dt)*mun(j,dt)*uPasser[j];
    initPote[j]=initPot(j,dt,xo,to,r);
    denom[j]=mun(j,dt)*kvtntn();
    q[j]=(-0.5*uPasser[j]-sLays[j]+dLays[j]+munValsD[j]+initPote[j])/denom[j];
    diff[j]=(q[j]-qKnownList[j])/qKnownList[j]*100;
    printf("%s\n",title13);
    printf("%f\n",(q[j]-qKnownList[j])/qKnownList[j]*100);
  }
  
  /*Build An Initial Potential Testing Module*/
  double ans2=initPot(10,dt,xo,to,r);
  double analyt=0.622853;
  printf("%s\n",title19);
  printf("%s\n",title21);
  printf("%f\n",ans2);
  printf("%s\n",title20);
  printf("%f\n",analyt);
  
  /*This is proof that if we feed the single layer scheme the correct q value at each arrival, then the scheme is accurate*/
  
  int testpoint=10;
  double ans3=sLay2(testpoint,dt,r,qTest)+mun(testpoint,dt)*qTest[testpoint];
  double analyt2=qTestExact(testpoint*dt);
  printf("%s\n",title22);
  printf("%s\n",title21);
  printf("%10.9f\n",ans3);
  printf("%s\n",title20);
  printf("%10.9f\n",analyt2);
  
}


