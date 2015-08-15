/* 
 * File:   SH.h
 * Author: antonermakov
 *
 * Created on April 21, 2012, 4:21 PM
 */

#ifndef SH_H
#define	SH_H
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <random>

#include "local_math.h"

class SH
{
private:
    double RefRad;
    double R;
    double mu;
    double* C;
    double* S;
    unsigned int MaxHarmDeg;
    double* M;

    double VWDiagonalRecursion(int m,double x, double y, double rsquared, double PreviousVW,double PreviousWV);
    double VWVerticalRecursion(int n, int m, double z, double rsquared, double PreviousVW,double PrePreviousVW);
    double VerticalSum(int m,double z,double rsquared,double Vmm,double Wmm, int MaxDeg);
    double NormCoef(int,int);
    double factorial(int);

public:
    SH(void);
    SH(double, double, double, double);
    SH(char*);
    double getC(int,int);
    double getS(int,int);
    void setC(int,int,double);
    void setS(int,int,double);
    double getRefRad();
    double getmu();
    double ExpandShape(double lat, double lon, int MaxDeg);
    double Expand(double lat,double lon, double r, int MaxDeg);
    double* ExpandPotential(double* x, double* y, double* z,int NumberOfPoints, int MaxDeg);

    void GenerateRandom(double intercept, double slope, double R);
    void ComputePowerSpectrum(void);
    void PrintPowerSpectrum(char*);

    void PrintSH(char*);
    void PrintGrid(char* filename, unsigned int nlat, unsigned int nlon);

    int CreateGrid(double FiStep, double LambdaStep, double R, double* x, double* y, double* z);
    
};



#endif	/* SH_H */

