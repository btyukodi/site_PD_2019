#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//change functions to macros. it'll be faster


//returns an integer random number from {0,1, ... ,rmax-1}
long randint(long rmax){
	return (long)((double)rand()/(RAND_MAX+1.0)*rmax);
}

//returns an integer random number from {rmin, rmin+1, ... ,rmax-1}
long randint(long rmin, long rmax){
	return rmin + (long)((double)rand()/(RAND_MAX+1.0)*(rmax-rmin));
}

//overload rand()
double rand(double rmax){
	return (double)rand()/(RAND_MAX+1.0) * rmax;
}

double rand(double rmin, double rmax){
	return rmin + (double)rand()/(RAND_MAX+1.0) * (rmax-rmin);

}


