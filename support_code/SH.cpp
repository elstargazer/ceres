#include "SH.h"

SH::SH(void)
{

	RefRad = 1.0;
	R = 1.0;
	MaxHarmDeg = 0;
	mu = 1.0;

	M = new double[MaxHarmDeg+1];

	unsigned int HarmCoefNumel=(MaxHarmDeg+1)*(MaxHarmDeg+2)/2;
	C = new double[HarmCoefNumel];
	S = new double[HarmCoefNumel];

	for(unsigned int i = 0; i<HarmCoefNumel; i++)
	{
	    C[i]=0;
	    S[i]=0;
	}
}

SH::SH(double RefRadi, double Ri, double nmax, double GM)
{

	RefRad = RefRadi;
	R = Ri;
	MaxHarmDeg = nmax;
	mu = GM;

	M = new double[MaxHarmDeg+1];

	unsigned int HarmCoefNumel=(MaxHarmDeg+1)*(MaxHarmDeg+2)/2;
    C = new double[HarmCoefNumel];
    S = new double[HarmCoefNumel];

    for(unsigned int i = 0; i<HarmCoefNumel; i++)
    {
    	C[i]=0;
    	S[i]=0;
    }

}

SH::SH(char* filename)
{
    FILE* pFile=fopen(filename,"rt");
    int StringLength=140;
    char* str = new char[StringLength];
    fgets(str,StringLength,pFile);
    //printf(str);
    char * next;
    
    RefRad=strtod(str,&next)*1000;   
    mu=strtod(next,NULL)*pow(10,9);    
    
    MaxHarmDeg = 270;
    M = new double[MaxHarmDeg+1];
    
    int NumberOfStrings=36856;
    int HarmCoefNumel=(MaxHarmDeg+1)*(MaxHarmDeg+2)/2;
    
    C = new double[HarmCoefNumel];
    S = new double[HarmCoefNumel];
    
    C[0]=1;
    S[0]=0;      
    
    int n;
    int m;   
    
    for(int i=1;i<HarmCoefNumel;i++)
    {
        fgets(str,StringLength,pFile);
        
        n=strtod(str,&next);
        m=strtod(next,&next);
        C[i]=strtod(next,&next);
        S[i]=strtod(next,&next);
        
        //printf("i=%i, n=%i, m=%i, C=%LE, S=%LE\n",i,n,m,C[i],S[i]);
        //printf(str);
    };

    fclose(pFile);
};

double SH::getC(int n, int m)
{
    int ZonalTermIndex=(n+1)*(n+2)/2-n-1;
    return(C[ZonalTermIndex+m]);
};

void SH::setC(int n, int m, double a)
{
    int ZonalTermIndex=(n+1)*(n+2)/2-n-1;
    C[ZonalTermIndex+m] = a;
};

double SH::getS(int n, int m)
{
    int ZonalTermIndex=(n+1)*(n+2)/2-n-1;
    return(S[ZonalTermIndex+m]);    
};

void SH::setS(int n, int m, double a)
{
    int ZonalTermIndex=(n+1)*(n+2)/2-n-1;
    S[ZonalTermIndex+m] = a;
};

double SH::getRefRad()
{
    return RefRad;
};

double SH::getmu()
{
    return mu;
}

double SH::Expand(double lat, double lon, double r, int MaxDeg)
{
	double x = r*cos(lat)*cos(lon);
	double y = r*cos(lat)*sin(lon);
	double z = r*sin(lat);

    double Q=0;
    
    double V00=1.0/r;
    double W00=0;
    
    double rsquared = r*r;

    double V1;
    double W1;
    double V0=V00;
    double W0=W00;
    
    Q+=VerticalSum(0,z,rsquared,V0,W0,MaxDeg);    
    
    for(int m=1;m<=MaxDeg-2;m++)        
    {
        V1=VWDiagonalRecursion(m,x,y,rsquared,V0,W0);
        W1=VWDiagonalRecursion(m,x,y,rsquared,W0,-V0);
        
        //printf("V(%i,%i)=%20.15LE\n\n",m,m,V1);
        //printf("W(%i,%i)=%20.15LE\n\n",m,m,W1);    
        
        Q+=VerticalSum(m,z,rsquared,V1,W1,MaxDeg);
        
        W0=W1;
        V0=V1;        
    };
    
    V1=VWDiagonalRecursion(MaxDeg-1,x,y,rsquared,V0,W0);
    W1=VWDiagonalRecursion(MaxDeg-1,x,y,rsquared,W0,-V0);    
    
    V0=V1;
    W0=W1;
    
    double Vlast=VWDiagonalRecursion(MaxDeg,x,y,rsquared,V0,W0);
    double Wlast=VWDiagonalRecursion(MaxDeg,x,y,rsquared,W0,-V0);
    
    V1=VWVerticalRecursion(MaxDeg,MaxDeg-1,z,rsquared,V0,0);
    W1=VWVerticalRecursion(MaxDeg,MaxDeg-1,z,rsquared,W0,0);  
    
    Q+=(getC(MaxDeg-1,MaxDeg-1)*V0+getS(MaxDeg-1,MaxDeg-1)*W0)/NormCoef(MaxDeg-1,MaxDeg-1);
    Q+=(getC(MaxDeg,MaxDeg-1)*V1+getS(MaxDeg,MaxDeg-1)*W1)/NormCoef(MaxDeg,MaxDeg-1);
    Q+=(getC(MaxDeg,MaxDeg)*Vlast+getS(MaxDeg,MaxDeg)*Wlast)/NormCoef(MaxDeg,MaxDeg); 
    
    double a=(getC(MaxDeg,MaxDeg)*Vlast+getS(MaxDeg,MaxDeg)*Wlast)/NormCoef(MaxDeg,MaxDeg);
    
    return Q;
    
};

double SH::VWDiagonalRecursion(int m,double x, double y, double rsquared, double PreviousVW, double PreviousWV)
{    
   return (2*m-1)*RefRad*(x*PreviousVW-y*PreviousWV)/rsquared;
};


double SH::VWVerticalRecursion(int n, int m, double z, double rsquared, double PreviousVW, double PrePreviousVW)
{    
    return (2*n-1)*z*RefRad*PreviousVW/(n-m)/rsquared-(n+m-1)*RefRad*RefRad*PrePreviousVW/rsquared/(n-m);
};

double SH::VerticalSum(int m, double z, double rsquared, double Vmm, double Wmm, int MaxDeg)
{
    double Qv=0;
    
    double V0=Vmm;
    double W0=Wmm;
    
    double V1=VWVerticalRecursion(m+1,m,z,rsquared,V0,0);
    double W1=VWVerticalRecursion(m+1,m,z,rsquared,W0,0);
    
    Qv+=(getC(m,m)*V0+getS(m,m)*W0)/NormCoef(m,m);
    
    //printf("n=%i;m=%i; term=%LE\n",m,m,(getC(m,m)*V0+getS(m,m)*W0)/NormCoef(m,m)*mu/RefRad);
    
    Qv+=(getC(m+1,m)*V1+getS(m+1,m)*W1)/NormCoef(m+1,m);  
    
    //printf("n=%i;m=%i; term=%LE\n\n",m+1,m,(getC(m+1,m)*V1+getS(m+1,m)*W1)/NormCoef(m+1,m)*mu/RefRad);
    
    double Vnew;
    double Wnew;
        
    for(int n=m+2;n<=MaxDeg;n++)
    {
        Vnew=VWVerticalRecursion(n,m,z,rsquared,V1,V0);
        Wnew=VWVerticalRecursion(n,m,z,rsquared,W1,W0);

        Qv+=(getC(n,m)*Vnew+getS(n,m)*Wnew)/NormCoef(n,m);

        //printf("n=%i;m=%i; term=%LE\n\n",n,m,(getC(n,m)*Vnew+getS(n,m)*Wnew)/NormCoef(n,m)*mu/RefRad);
        
        V0=V1;
        W0=W1; 
        V1=Vnew;
        W1=Vnew;
   
        //printf("n=%i; m=%i\n",n,m);
        
    }  
    return Qv;
    
}

double SH::NormCoef(int n, int m)
{
    
    int delta;    
    if (m==0)
    delta=1;
    else
    delta=0;
    return sqrtl( factorial(n+m)/(2-delta)/(2*n+1)/factorial(n-m)    );
}

double SH::factorial(int num)
{
    if (num==1)
    return 1;
    else if (num==0)
    return 1;
    return factorial(num-1)*num;
}

void SH::GenerateRandom(double intercept, double slope, double Ri)
{
	R = Ri;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-0.5, 0.5);

	for(unsigned int n=2; n<=MaxHarmDeg; n++)
		for(unsigned int m = 0; m <=n; m++)
		{
			setC(n,m,dis(gen));
			setS(n,m,0.0);
		}

	ComputePowerSpectrum();

	double oldC;
	double newC;
	double oldS;
	double newS;
	double phase;

	for(unsigned int n=2; n<=MaxHarmDeg; n++)
		for(unsigned int m = 0; m <=n; m++)
		{
			oldC = getC(n,m);
			newC = oldC*sqrt(intercept*pow(n,slope))/sqrt(M[n]);
			setC(n,m,newC);

			// randomize phase
//			oldC = getC(n,m);
//			oldS = getS(n,m);
//
//			phase = dis(gen)*TWOPI;
//
//			newC = cos(phase)*oldC - sin(phase)*oldS;
//			newS = sin(phase)*oldC + cos(phase)*oldS;
//
//			setC(n,m,newC);
//			setS(n,m,newS);

		}
	ComputePowerSpectrum();
}

void SH::ComputePowerSpectrum()
{
	for(unsigned int n=0; n<=MaxHarmDeg; n++)
	{
		M[n] = 0;
		for(unsigned int m = 0; m <=n; m++)
			M[n]+=pow(getC(n,m),2) + pow(getS(n,m),2);
		M[n]/=(2*n+1);
	}
}

void SH::PrintPowerSpectrum(char* filename)
{
	FILE* pFile = fopen(filename,"w");
	for(unsigned int n = 0; n<=MaxHarmDeg; n++)
	{
		fprintf(pFile,"%i, %E\n",n,M[n]);
//		printf("%u, %E\n",n,M[n]);
	}
	fclose(pFile);
}

double SH::ExpandShape(double lat, double lon, int MaxDeg)
{
	double r;
	r = Expand(lat, lon, 1.0, MaxDeg);
	r += R;
	return r;
}

void SH::PrintGrid(char* filename, unsigned int nlat, unsigned int nlon)
{
	FILE* pFile = fopen(filename,"w");

	double lat = -PI/2;
	double lon;

	double dlat = PI/nlat;
	double dlon = TWOPI/nlon;

	double r;

	for(unsigned int i=0; i<nlat; i++)
	{
		lon = -PI;
		for(unsigned int j=0; j<nlon; j++)
		{
			r = ExpandShape(lat, lon, MaxHarmDeg);
			fprintf(pFile,"%13.6E ", r);
			lon+=dlon;
		}
		fprintf(pFile,"\n");
//		printf("%f\n",lat);
		lat+=dlat;
	}
	fclose(pFile);
}

