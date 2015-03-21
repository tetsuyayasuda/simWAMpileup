//-------------------------------
//     powCalc.hh
//-------------------------------
// 2013-09-13 yasuda

#ifndef POWCALC_HH
#define POWCALC_HH

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

double powConst(double x1, double y1, double x2, double y2);
double powIndex(double x1, double y1, double x2, double y2);
double powInte(double x1, double x2, double A, double B);
double powValue(double x, double A, double B);
double powInvValue(double y, double A, double B);

double linearConst(double x1, double y1, double x2, double y2);
double linearIndex(double x1, double y1, double x2, double y2);
double linearInte(double x1, double x2, double A, double B);
double linearValue(double x, double A, double B);
double linearInvValue(double y, double A, double B);

#endif
