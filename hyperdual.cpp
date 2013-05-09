#include "./hyperdual.hpp"

#include <iostream>
#include <stdio.h>
#include <math.h>
using namespace std;

HyperDual::HyperDual()
{
	f0 = 0.0;
	f1 = 0.0;
	f2 = 0.0;
	f12 = 0.0;
}
HyperDual::HyperDual(double x1,double x2,double x3,double x4)
{
	f0 = x1;
	f1 = x2;
	f2 = x3;
	f12 = x4;
}
HyperDual::HyperDual(double x1)
{
	f0 = x1;
	f1 = 0.0;
	f2 = 0.0;
	f12 = 0.0;
}

void HyperDual::setvalues(double x1,double x2,double x3,double x4)
{
	f0 = x1;
	f1 = x2;
	f2 = x3;
	f12 = x4;
}
void HyperDual::view(void)
{
	printf("%g  +  %g i  +  %g j  +  %g ij\n",f0,f1,f2,f12);
}

double & HyperDual::realpart(void)
{
	return f0;
}
double & HyperDual::ipart(void)
{
	return f1;
}
double & HyperDual::jpart(void)
{
	return f2;
}
double & HyperDual::ijpart(void)
{
	return f12;
}
const double & HyperDual::realpart(void) const 
{
	return f0;
}
const double & HyperDual::ipart(void) const
{
	return f1;
}
const double & HyperDual::jpart(void) const
{
	return f2;
}
const double & HyperDual::ijpart(void) const
{
	return f12;
}

HyperDual HyperDual::operator+ () const
{
	return *this;
}
HyperDual HyperDual::operator+ (const HyperDual rhs) const
{
	HyperDual temp;
	temp.f0 = f0 + rhs.f0;
	temp.f1 = f1 + rhs.f1;
	temp.f2 = f2 + rhs.f2;
	temp.f12 = f12 + rhs.f12;
	return temp;
}

HyperDual operator+ (const double lhs, const HyperDual rhs)
{
	HyperDual temp;
	temp.f0 = lhs + rhs.f0;
	temp.f1 = rhs.f1;
	temp.f2 = rhs.f2;
	temp.f12 = rhs.f12;
	return temp;
}


HyperDual HyperDual::operator- () const
{
	HyperDual temp;
	temp.f0 = -f0;
	temp.f1 = -f1;
	temp.f2 = -f2;
	temp.f12 = -f12;
	return temp;
}

HyperDual HyperDual::operator- (const HyperDual rhs) const
{
	HyperDual temp;
	temp.f0 = f0 - rhs.f0;
	temp.f1 = f1 - rhs.f1;
	temp.f2 = f2 - rhs.f2;
	temp.f12 = f12 - rhs.f12;
	return temp;
}
HyperDual operator- (const double lhs, const HyperDual rhs)
{
	HyperDual temp;
	temp.f0 = lhs - rhs.f0;
	temp.f1 = -rhs.f1;
	temp.f2 = -rhs.f2;
	temp.f12 = -rhs.f12;
	return temp;
}

HyperDual HyperDual::operator* (const HyperDual rhs) const
{
	HyperDual temp;
	temp.f0 = f0*rhs.f0;
	temp.f1 = f0*rhs.f1 + f1*rhs.f0;
	temp.f2 = f0*rhs.f2 + f2*rhs.f0;
	temp.f12 = f0*rhs.f12 + f1*rhs.f2 + f2*rhs.f1 + f12*rhs.f0;
	return temp;
}

HyperDual operator* (const double lhs, const HyperDual rhs)
{
	HyperDual temp;
	temp.f0 = lhs*rhs.f0;
	temp.f1 = lhs*rhs.f1;
	temp.f2 = lhs*rhs.f2;
	temp.f12 = lhs*rhs.f12;
	return temp;
}
HyperDual operator/ (const HyperDual lhs, const HyperDual rhs)
{
	HyperDual temp,inv;
	inv = pow(rhs,-1);
	temp = lhs*inv;
	return temp;
}

HyperDual operator/ (const double lhs, const HyperDual rhs)
{
	HyperDual temp,inv;
	inv = pow(rhs,-1);
	temp = lhs*inv;
	return temp;
}

HyperDual operator/ (const HyperDual lhs, const double rhs)
{
	HyperDual temp;
	double inv;
	inv = 1.0/rhs;
	temp.f0 = inv*lhs.f0;
	temp.f1 = inv*lhs.f1;
	temp.f2 = inv*lhs.f2;
	temp.f12 = inv*lhs.f12;
	return temp;
}

HyperDual& HyperDual::operator+= (HyperDual rhs)
{
	f0 += rhs.f0;
	f1 += rhs.f1;
	f2 += rhs.f2;
	f12 += rhs.f12;
	return *this;
}
HyperDual& HyperDual::operator-= (HyperDual rhs)
{
	f0 -= rhs.f0;
	f1 -= rhs.f1;
	f2 -= rhs.f2;
	f12 -= rhs.f12;
	return *this;
}
HyperDual& HyperDual::operator*= (HyperDual rhs)
{
	double tf0,tf1,tf2,tf12;
	tf0=f0;
	tf1=f1;
	tf2=f2;
	tf12=f12;
	f0 = tf0*rhs.f0;
	f1 = tf0*rhs.f1 + tf1*rhs.f0;
	f2 = tf0*rhs.f2 + tf2*rhs.f0;
	f12 = tf0*rhs.f12 + tf1*rhs.f2 + tf2*rhs.f1 + tf12*rhs.f0;
	return *this;
}

HyperDual& HyperDual::operator*= (double rhs)
{
	f0 *= rhs;
	f1 *= rhs;
	f2 *= rhs;
	f12 *= rhs;
	return *this;
}
HyperDual& HyperDual::operator/= (double rhs)
{
	f0 /= rhs;
	f1 /= rhs;
	f2 /= rhs;
	f12 /= rhs;
	return *this;
}

HyperDual pow (HyperDual x, double a)
{
	HyperDual temp;
	double deriv,xval,tol;
	xval = x.f0;
	tol = 1e-15;
	if (fabs(xval) < tol)
	{
		if (xval >= 0)
			xval = tol;
		if (xval < 0)
			xval = -tol;
	}
	deriv = a*pow(xval,(a-1));
	//temp.f0 = pow(xval,a);
	temp.f0 = pow(x.f0,a);  //Use actual x value, only use tol for derivs
	temp.f1 = x.f1*deriv;
	temp.f2 = x.f2*deriv;
	temp.f12 = x.f12*deriv + a*(a-1)*x.f1*x.f2*pow(xval,(a-2));
	
	return temp;
}
HyperDual pow (HyperDual x, HyperDual a)
{
	return exp(a*log(x));
}
HyperDual exp(HyperDual x)
{
	HyperDual temp;
	double deriv;
	deriv = exp(x.f0);
	temp.f0 = deriv;
	temp.f1 = deriv*x.f1;
	temp.f2 = deriv*x.f2;
	temp.f12 = deriv*(x.f12 + x.f1*x.f2);
	return temp;
}
HyperDual log(HyperDual x)
{
	HyperDual temp;
	double deriv1,deriv2;
	deriv1 = x.f1/x.f0;
	deriv2 = x.f2/x.f0;
	temp.f0 = log(x.f0);
	temp.f1 = deriv1;
	temp.f2 = deriv2;
	temp.f12 = x.f12/x.f0 - (deriv1*deriv2);
	return temp;
}

HyperDual sin(HyperDual x)
{
	HyperDual temp;
	double funval,deriv;
	funval = sin(x.f0);
	deriv = cos(x.f0);
	temp.f0 = funval;
	temp.f1 = deriv*x.f1;
	temp.f2 = deriv*x.f2;
	temp.f12 = deriv*x.f12 - funval*x.f1*x.f2;
	return temp;
}

HyperDual cos(HyperDual x)
{
	HyperDual temp;
	double funval,deriv;
	funval = cos(x.f0);
	deriv = -sin(x.f0);
	temp.f0 = funval;
	temp.f1 = deriv*x.f1;
	temp.f2 = deriv*x.f2;
	temp.f12 = deriv*x.f12 - funval*x.f1*x.f2;
	return temp;
}
HyperDual tan(HyperDual x)
{
	HyperDual temp;
	double funval,deriv;
	funval = tan(x.f0);
	deriv  = funval*funval + 1.0;
	temp.f0 = funval;
	temp.f1 = deriv*x.f1;
	temp.f2 = deriv*x.f2;
	temp.f12 = deriv*x.f12 + x.f1*x.f2*(2*funval*deriv);
	return temp;
}

HyperDual asin(HyperDual x)
{
	HyperDual temp;
	double funval,deriv1,deriv;
	funval = asin(x.f0);
	deriv1 = 1.0-x.f0*x.f0;
	deriv = 1.0/sqrt(deriv1);
	temp.f0 = funval;
	temp.f1 = deriv*x.f1;
	temp.f2 = deriv*x.f2;
	temp.f12 = deriv*x.f12 + x.f1*x.f2*(x.f0/pow(deriv1,-1.5));
	return temp;
}
HyperDual acos(HyperDual x)
{
	HyperDual temp;
	double funval,deriv1,deriv;
	funval = acos(x.f0);
	deriv1 = 1.0-x.f0*x.f0;
	deriv = -1.0/sqrt(deriv1);
	temp.f0 = funval;
	temp.f1 = deriv*x.f1;
	temp.f2 = deriv*x.f2;
	temp.f12 = deriv*x.f12 + x.f1*x.f2*(-x.f0/pow(deriv1,-1.5));
	return temp;
}
HyperDual atan(HyperDual x)
{
	HyperDual temp;
	double funval,deriv1,deriv;
	funval = atan(x.f0);
	deriv1 = 1.0+x.f0*x.f0;
	deriv = 1.0/deriv1;
	temp.f0 = funval;
	temp.f1 = deriv*x.f1;
	temp.f2 = deriv*x.f2;
	temp.f12 = deriv*x.f12 + x.f1*x.f2*(-2*x.f0/(deriv1*deriv1));
	return temp;
}

HyperDual sqrt(HyperDual x)
{
	return pow(x,0.5);
}

HyperDual fabs(HyperDual x)
{
	HyperDual temp;
	if (x < 0.0)
		temp = -x;
	else
		temp = x;
	return temp;
}

HyperDual max(HyperDual x1, HyperDual x2)
{
	HyperDual temp;
	if (x1>x2)
		temp = x1;
	else
		temp = x2;
	return temp;
}
HyperDual max(HyperDual x1, double x2)
{
	HyperDual temp;
	if (x1>x2)
		temp = x1;
	else
		temp = x2;
	return temp;
}
HyperDual max(double x1, HyperDual x2)
{
	HyperDual temp;
	if (x1>x2)
		temp = x1;
	else
		temp = x2;
	return temp;
}

HyperDual min(HyperDual x1, HyperDual x2)
{
	HyperDual temp;
	if (x1<x2)
		temp = x1;
	else
		temp = x2;
	return temp;
}
HyperDual min(HyperDual x1, double x2)
{
	HyperDual temp;
	if (x1<x2)
		temp = x1;
	else
		temp = x2;
	return temp;
}
HyperDual min(double x1, HyperDual x2)
{
	HyperDual temp;
	if (x1<x2)
		temp = x1;
	else
		temp = x2;
	return temp;
}

bool operator> (HyperDual lhs, HyperDual rhs)
{
	return (lhs.f0 > rhs.f0);
}
bool operator> (double lhs, HyperDual rhs)
{
	return (lhs > rhs.f0);
}
bool operator> (HyperDual lhs, double rhs)
{
	return (lhs.f0 > rhs);
}
bool operator>= (HyperDual lhs, HyperDual rhs)
{
	return (lhs.f0 >= rhs.f0);
}
bool operator>= (double lhs, HyperDual rhs)
{
	return (lhs >= rhs.f0);
}
bool operator>= (HyperDual lhs, double rhs)
{
	return (lhs.f0 >= rhs);
}
bool operator< (HyperDual lhs, HyperDual rhs)
{
	return (lhs.f0 < rhs.f0);
}
bool operator< (double lhs, HyperDual rhs)
{
	return (lhs < rhs.f0);
}
bool operator< (HyperDual lhs, double rhs)
{
	return (lhs.f0 < rhs);
}
bool operator<= (HyperDual lhs, HyperDual rhs)
{
	return (lhs.f0 <= rhs.f0);
}
bool operator<= (double lhs, HyperDual rhs)
{
	return (lhs <= rhs.f0);
}
bool operator<= (HyperDual lhs, double rhs)
{
	return (lhs.f0 <= rhs);
}
bool operator== (HyperDual lhs, HyperDual rhs)
{
	return (lhs.f0 == rhs.f0);
}
bool operator== (double lhs, HyperDual rhs)
{
	return (lhs == rhs.f0);
}
bool operator== (HyperDual lhs, double rhs)
{
	return (lhs.f0 == rhs);
}
bool operator!= (HyperDual lhs, HyperDual rhs)
{
	return (lhs.f0 != rhs.f0);
}
bool operator!= (double lhs, HyperDual rhs)
{
	return (lhs != rhs.f0);
}
bool operator!= (HyperDual lhs, double rhs)
{
	return (lhs.f0 != rhs);
}

ostream& operator<<(ostream& output, const HyperDual& rhs)
{
        output << "(" << rhs.f0 << ","<< rhs.f1 << ","<< rhs.f2 << ","<< rhs.f12 << ")";
        return output;
}



