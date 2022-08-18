#ifndef RANDOM_H
#define RANDOM_H

//returns an integer random number from {0,1, ... ,rmax-1}
long randint(long rmax);

//returns an integer random number from {rmin, rmin+1, ... ,rmax-1}
long randint(long rmin, long rmax);

//returns a double random number from [0, rmax); overloads rand()
double rand(double rmax);

//returns a double random number from [rmin, rmax); overloads rand()
double rand(double rmin, double rmax);


#endif
