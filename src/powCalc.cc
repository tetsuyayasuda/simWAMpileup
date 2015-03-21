#include "powCalc.hh"


double powConst(double x1, double y1, double x2, double y2){ // y = A*x**B
  double B = log10(y1/y2) / log10(x1/x2);
  double A = y1/pow(x1, B);
  return A;
}

double powIndex(double x1, double y1, double x2, double y2){ // y = A*x**B
  double B = log10(y1/y2) / log10(x1/x2);
  return B;
}

double powInte(double x1, double x2, double A, double B){ // y = A*x**B
  return A * (pow(x2, B+1.0) - pow(x1, B+1.0))/(B+1.0);
}

double powValue(double x, double A, double B){
  return A * pow(x, B);
}

double powInvValue(double y, double A, double B){
  return pow(y/A, 1.0/B);
}




double linearConst(double x1, double y1, double x2, double y2){ // y = A*x +B
  return (y2-y1)/(x2-x1);
}

double linearIndex(double x1, double y1, double x2, double y2){ // y = A*x +B
  double A = (y2-y1)/(x2-x1);
  double B1 = y1 - A*x1;
  double B2 = y2 - A*x2;
  return (B1+B2)*0.5;
}

double linearInte(double x1, double x2, double A, double B){ // y = A*x +B
  return ((A*0.5*x2*x2) + B*x2) - ((A*0.5*x1*x1) + B*x1);
}

double linearValue(double x, double A, double B){ // y = A*x +B
  return (A*x) + B;
}

double linearInvValue(double y, double A, double B){ // y = A*x +B
  return (y - B)/A;
}

