#ifndef HISTOGRAM_H
#define HISTOGRAM_H

int update_hist(array<long>& h, double new_xvalue);

int update_hist(array<long>& h, double new_xvalue, double new_yvalue);

//overload the same function to update values from an array (h.sizeY should be 1)
int update_hist(array<long>& h, array<double>& new_xvalue);

//overload the same function to update values from 2 arrays; the dimensions of the 2 arrays should be the same
int update_hist(array<long>& h, array<double>& new_xvalue, array<double>& new_yvalue);


#endif
