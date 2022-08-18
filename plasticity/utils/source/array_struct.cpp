#ifndef ARRAY_STRUCT_H
#define ARRAY_STRUCT_H

#include <stdio.h>
#include <stdlib.h>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>
#include<iomanip>
#include<math.h>
#include<omp.h>
#include"array.h"
#include"fft.h"
#include"random.h"
#include"fft2D.h"

//array type
//type T can be either int/long or float/double
/*usage examples:
definition:
 array<double> my_1D_array(128)
 array<double> sigma(256,256);
 array<double> G("green256_1_r2",256,256);
 array<long> histogram(10000,1,0.01)

accessing elements:
	- set pointer to values: T * line = sigma.values[line_nr];
	- use overloaded operators: sigma(i,j) = values[i*sizeX+j];
*/

template <typename T>
struct array{
	//dimensions along X and Y directions
	long sizeX, sizeY;

	long len; //total size, sizeX*sizeY

	//resolutions along X and Y directions; optional, used to convert between indices and real coordinates values correspond to
	double resX, resY;

	//lowest coordinates; optional, used to convert between indices and real coordinates values correspond to
	double minX, minY;

	//set it to 'p' to treat array as periodic, 'n' to set it non-periodic
	//array is periodic by default
	char periodic;

	//actual values
	T * values;

	//i0, j0; a shift-nek lehet ertelme, pl. skalar szorzasnal G.shift(2,5); Sigma+=eta*G;

	//constructor1; 
	array(long sx, long sy=1, double rx=1.0, double ry=1.0, double mx=0.0, double my=0.0, char pp ='p'){
		sizeX=sx;
		sizeY=sy;
		resX = rx;
		resY = ry;
		minX = mx;
		minY = my;	
		len = sizeX*sizeY;
		values = new T[len];//alloc<T>(sizeX, sizeY);
		//zeros<T>(values, sizeX, sizeY);
		periodic = pp;
	
	}


	//constructor2; load from file
	array(std::string file_name, long sx=1, long sy=1, double rx=1.0, double ry=1.0, double mx=0.0, double my=0.0, char pp='p'){
		sizeX=sx;
		sizeY=sy;
		resX = rx;
		resY = ry;
		minX = mx;
		minY = my;
		len = sizeX*sizeY;
		values = new T[len];//values = alloc<T>(sizeX, sizeY);
		loadtxt<T>(file_name, values, sizeX, sizeY);
		periodic = pp;	
		
	}

	//copy constructor
	array (const array<T> &m){
		//printf("--copyctr\n");
		sizeX=m.sizeX;
		sizeY=m.sizeY;
		resX=m.resX;
		resY=m.resY;
		minX=m.minX;
		minY=m.minY;
		len=m.len;
		periodic='p';
		values= new T[len]; //alloc<T>(sizeX, sizeY);
		long i; 
		#pragma omp parallel for if(len>parallel_threshold) private(i)
		for (i=0; i<len; i++) {
				values[i] = m.values[i];

		}		
	}

	//empty constructor?
	array(){

	}

	//init for the empty constructor
	void init(long sx, long sy=1, double rx=1.0, double ry=1.0, double mx=0.0, double my=0.0, char pp ='p'){
		sizeX=sx;
		sizeY=sy;
		resX = rx;
		resY = ry;
		minX = mx;
		minY = my;	
		len = sizeX*sizeY;
		values = new T[len];//alloc<T>(sizeX, sizeY);
		//zeros<T>(values, sizeX, sizeY);
		periodic = pp;
	
	}


	//------assignment operator---------------------------------------
	array<T>& operator = (const array<T>& m) {
		if (this == &m) { 
		// avoid self-assignment
		return *this;
		} 
		else{	
			long i; 
			#pragma omp parallel for private(i) if(len>parallel_threshold)
			for (i=0; i<len; i++) {				
				values[i] = m.values[i];				
			}

		return *this;
		}
	}
	

	//-----assignment operator----------------------------------------
	//fill in with number
	array<T>& operator = (const T m) {
			long i; 
			#pragma omp parallel for private(i) if(len>parallel_threshold)
			for (i=0; i<len; i++) {				
				values[i] = m;
				
			}

		return *this;
		
	}

	//-------operators to access elements as array(i,j)---------------
	 T& operator()(long i, long j) {
		return values[i*sizeX+j];
	}

	T operator()(long i, long j) const{
		return values[i*sizeX+j];
	}

	//for 1D arrays
	T& operator()(long j){
		return values[j];
	}

	T operator()(long j) const{
		return values[j];
	}

	//--------increment/decrement/divide/multiply by matrix ------------------------------------
	array<T>& operator += (const array<T>& m) {
			long i;
			#pragma omp parallel for private(i) if(len>parallel_threshold)
			for (i=0; i<len; i++) {				
				values[i] += m.values[i];				
		}
		return *this;
	}

	array<T>& operator += (const T m) {
			long i; 
			#pragma omp parallel for private(i) if(len>parallel_threshold)
			for (i=0; i<len; i++) {				
					values[i] += m;
				
		}
		return *this;
	}


	array<T>& operator -= (const array<T>& m) {
			long i; 
			#pragma omp parallel for private(i) if(len>parallel_threshold)
			for (i=0; i<len; i++) {
					values[i] -= m.values[i];
				
		}
		return *this;
	}

	array<T>& operator -= (const T m) {
			long i; 
			#pragma omp parallel for private(i) if(len>parallel_threshold)
			for (i=0; i<len; i++) {
					values[i] -= m;
			}
		return *this;
	}


	array<T>& operator /= (const array<T>& m) {
			long i;
			#pragma omp parallel for private(i) if(len>parallel_threshold)
			for (i=0; i<len; i++) {
					values[i] /= m.values[i];			
			}
		return *this;
	}

	array<T>& operator /= (const T m) {
			long i;
			#pragma omp parallel for private(i) if(len>parallel_threshold)
			for (i=0; i<len; i++) {
					values[i] /= m;				
			}
		return *this;
	}

		
	array<T>& operator *= (const array<T>& m) {
			long i; 
			#pragma omp parallel for private(i) if(len>parallel_threshold)
			for (i=0; i<len; i++) {
					values[i] *= m.values[i];
			}
		return *this;
	}


		
	array<T>& operator *= (const T m) {
			long i; long j;
			#pragma omp parallel for private(i) if(len>parallel_threshold)
			for (i=0; i<len; i++) {
					values[i] *= m;			
			}
		return *this;
	}


	//destructor
	~array(){

		delete[] values;
		//dealloc(values, sizeX, sizeY);
		//delete variables too??
	}

	//returns sum of elements of the array
	T sum(){
		return sum_array<T>(values, len);
	}

	


	//returns average of elements; this is always double
	double avg(){
		double s = sum();
		s = s/(sizeX*sizeY);
		return s;
	}

	//shifts all the elements such that their average is 0
	//only use for float/double arrays
	void remove_drift(){
		T s = avg();	
		long i;
		#pragma omp parallel for private(i) if(len>parallel_threshold)
		for (i=0; i<len; i++){
				 values[i] -= s;
		}
	}

	//fill up the array with uniformly distributed random values from [vmin, vmax)
	void randomize(T vmin, T vmax){
		long i;
		for (i=0; i<len; i++){
				 values[i] = vmin + (T)((double)rand()/(RAND_MAX+1.0)*(vmax-vmin)) ;
		}
	}

	void randomize(T vmax){
		long i;
		for (i=0; i<len; i++){
				 values[i] = (T)((double)rand()/(RAND_MAX+1.0)*vmax) ;
		}
	}

	void randomize(){
		long i;
		for (i=0; i<len; i++){
				 values[i] = (T)((double)rand()/(RAND_MAX+1.0)) ;
		}
	}


	//saves array values as a text file
	void save(std::string file_name){
		savetxt<T>(file_name, values, sizeX, sizeY, 4);
	}

	//!!not tested!!
	void load(std::string file_name){
		loadtxt<T>(file_name, values, sizeX, sizeY);
	}

	//!!not tested!!
	//appends array values to a text file
	void append(std::string file_name){
		appendtxt<T>(file_name, values, sizeX, sizeY, 4);
	}

	//to be implemented; save dimensions, values, resolutions etc.
	void save_array(char file_name[81]){

	}

	//periodic boundary conditions for the array; puts all indices back to the first copy
	//use some sort of inline for this to be faster; eventually reduce maximum indices to remove the modulo operators;
	void pbc(long &i, long &j){
		if (i>0 && i>sizeY-1){
			i = i % sizeY;
		}
		if (j>0 && j>sizeX-1){
			j = j % sizeX;
		}
		if (i<0){
			i = sizeY-1 + ( (i+1) % sizeY);	
		}
		if (j<0){
			j = sizeX-1 + ( (j+1) % sizeX);
		}
	}

	//maximum, minimum functions
	T max(long &i, long &j){
		long k;
		T m;
		m = findmax(values, len, k);
	        i=k/sizeX; 
		j=k%sizeX; 
		return m;
	}
	
	T max(long &j){
		long k;
		T m;
		m = findmax(values, len, k);
		j=k;
		return m;
	}

	T max(){
		long k;
		T m;
		m = findmax(values, len, k);
		return m;
	}

	T min(long &i, long &j){
		long k;
		T m;
		m = findmin(values, len, k);
	        i=k/sizeX; 
		j=k%sizeX; 
		return m;
	}
	
	T min(long &j){
		long k;
		T m;
		m = findmin(values, len, k);
		j=k;
		return m;
	}

	T min(){
		long k;
		T m;
		m = findmin(values, len, k);
		return m;
	}


	//return profiles under angles 0,45,90,135
	//returns row vectors
	//use it for square arrays
	//costy!!! includes a constructor
	array<T> profile(long angle){
		long i;
		array<T> p(sizeY); //hopefully destroy it at return;
		p=0;
		if (sizeX != sizeY){
			return p;
		}		

		switch(angle){
			case 0:	{
					for (i=0; i<sizeX; i++){
						p(i) = values[sizeY/2*sizeX + i];
					}
					break;
				}
			case 45:{
					for (i=0; i<sizeX; i++){
						p(i) = values[i*sizeX+i];
					}
					break;
				}
			case 90:{
					for (i=0; i<sizeY; i++){
						p(i) = values[i*sizeX + sizeX/2];
					}
					break;
				}
			case 135:
				{
					for (i=0; i<sizeX; i++){
						p(sizeX-1-i) = values[(sizeY-1-i)*sizeX+i];
					}
					break;
				}
		}
		return p;
	}


	//return profiles under angles 0,45,90,135
	//returns row vectors
	//use it for square arrays
	//costy!!! includes a constructor
	array<T> profile(long angle, long offset){
		long i;
		array<T> p(sizeY); //hopefully destroy it at return;
		p=0;
		if (sizeX != sizeY){
			return p;
		}		

		switch(angle){
			case 0:	{
					for (i=0; i<sizeX; i++){
						p(i) = values[offset*sizeX + i];
					}
					break;
				}
			case 45:{
					for (i=0; i<offset; i++){
						p(i) = values[i*sizeX + (sizeY-offset+i)];
					}
					for (i=offset; i<sizeY; i++){
						p(i) = values[i*sizeX + (i-offset)];
					}
					break;
				}
			case 90:{
					for (i=0; i<sizeY; i++){
						p(i) = values[i*sizeX + offset];
					}
					break;
				}
			case 135:
				{
					for (i=sizeY-1; i>offset; i--){
						p(i) = values[i*sizeX + (sizeY+offset-i)];
					}
					for (i=offset; i>=0; i--){
						p(i) = values[i*sizeX + (offset-i)];
					}
					break;
				}
		}
		return p;
	}


	
	//make it parallel if necessary
	//sums up elements along one axis and returns a line or column vector
	array<T> sum(long axis){
		long i,j;
		switch(axis){
			//sum along the y axis, returns row vector
			case 1:{
				array<T> p(sizeX, 1);
				p = 0;
				for (i=0; i<sizeY; i++){
					for (j=0; j<sizeX; j++){
						p(j)+=values[i*sizeX+j];
					}
				}
				return p;
			}
			//sum along the x axis, returns column vector
			case 0:{
				array<T> p(1, sizeY);
				p = 0;
				for (i=0; i<sizeY; i++){
					for (j=0; j<sizeX; j++){
						p(i,0)+=values[i*sizeX+j];
					}
				}
				return p;
			}
		}
	}

	
	//returns a two times longer array with the fourier transform
	//apply for row vectors
	array<T> rfft1(int isign=1){
		int nn = pow(2, 1 + (int)(log(sizeX-1)/log(2))); //=2**(1+int(log2(length)))
		array<T> p(2*sizeX+1);
		p(0)=0;
		long i;
		//feed real part
		for (i=1; i<2*sizeX+1; i+=2){
			p(i) = values[i/2];
		}

		//feed imaginary part
		for (i=2; i<2*sizeX+1; i+=2){
			p(i) = 0;
		}

		four1(p.values, nn, isign);

		return p;
	}

	//overwrite row vector with its power spectrum
	void rps1(){
		array<T> p(2*sizeX+1);
		p = rfft1();
		long i;
		for (i=1; i<2*sizeX+1; i+=2){
			values[i/2] = p(i)*p(i) + p(i+1)*p(i+1);
		}
	}

	void rps2(){

		complex * y;
		y = (complex *) malloc( len*sizeof(complex) );
		long i;
		for (i=0; i<len; i++){
			y[i].Re = values[i];
			y[i].Im = 0.0;
		}

			
		fft2D(y, sizeX, sizeY, 1);

		for (i=0; i<len; i++){
			values[i] = y[i].Re*y[i].Re + y[i].Im*y[i].Im;
		}

		free(y);
	}

//	void four1(double data[], int nn, int isign)


	//shuffle the elements of the array
	void shuffle(){
		 long n = sizeX*sizeY;
		 long sw;
		 T temp;
		  for (long i=n-1; i>0; --i) {
			sw = randint(n);
			temp = values[i];
			values[i] = values[sw];
			values[sw]=temp;
		  }
	}

};


template <typename T> 
void testf(array<T> a1){
	a1.resX=125.001;
}

#endif
