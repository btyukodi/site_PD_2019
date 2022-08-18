/*
General functions related to 1D and 2D array manipulation. 
*/



#include <stdio.h>
#include <stdlib.h>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>
#include<iomanip>
#include<omp.h>
//because of templates
#include"array.h"


#ifndef ARRAY_C
#define ARRAY_C


//put this shit in some global file
#define parallel_threshold 2000


//allocates memory for 2D array of type T
//1st index for row, 2nd for column (Row-major order)
//sizeY: number of rows
//sizeX: number of columns
//if sizeY=1 => array is a row vector
template <typename T> 
T **alloc(long sizeX, long sizeY/*=1*/){
      long i;	
      T **array;  
      array = new T*[sizeY];
      for(i=0; i<sizeY; i++){
	      array[i] = new T[sizeX];
      }
      return array;
}

//deallocates memory for 2D array
template <typename T> 
void dealloc(T **array, long sizeX, long sizeY/*=1*/){
      long i;	
 
      for(i=0; i<sizeY; i++){
	      delete[] array[i];
      }
     delete[] array;
}


//makes all the elements of a 2D array of type T zero
template <typename T>
void zeros(T **array, long sizeX, long sizeY/*=1*/){
	long i, j;
	for (i=0; i<sizeY; i++){
		for (j=0; j<sizeX; j++){
			array[i][j]=0;
		}
	}
}

//finds maximum value and position 
template <typename T>
T findmax(T **array, long sizeX, long sizeY, long &i, long &j){
	long km, lm, k, l;
	T m=array[0][0];
	km=0; lm=0;
		for (k=0; k<sizeY; k++){
			for (l=0; l<sizeX; l++){
				 if(array[k][l]>m){
					m=array[k][l];
					km = k; lm = l;
				 }
			}
		}
		i=km;
		j=lm;
	return m;	
}


//finds maximum value and position 
//OPTIMIZE; CHECK OUT FINDMIN
//1D array
template <typename T>
T findmax(T *array, long len, long &i){
	long km, k;
	T m=array[0];
	km=0;
		#pragma omp parallel for private(k) if(len>parallel_threshold)
		for (k=0; k<len; k++){
				 #pragma omp flush(m)
				 if(array[k]>m){
					#pragma omp critical
				 	if(array[k]>m){
						m=array[k];
						km = k; 
					}
				 }			
		}
		i=km;
	return m;	
}


//finds minimum value and position 
template <typename T>
T findmin(T **array, long sizeX, long sizeY, long &i, long &j){
	long km, lm, k, l;
	T m=array[0][0];
	km=0; lm=0;
		for (k=0; k<sizeY; k++){
			for (l=0; l<sizeX; l++){
				 if(array[k][l]<m){
					m=array[k][l];
					km = k; lm = l;
				 }
			}
		}
		i=km;
		j=lm;
	return m;	
}


//finds maximum value and position 
//1D array
/*
template <typename T>
T findmin(T *array, long len, long &i){
	long km, k;
	T m=array[0];
	km=0;
		#pragma omp parallel for private(k) if(len>parallel_threshold)
		for (k=0; k<len; k++){
				 #pragma omp flush(m)
				 if(array[k]<m){
					#pragma omp critical
				 	if(array[k]<m){
						m=array[k];
						km = k; 
					}
				 }			
		}
		i=km;
	return m;	
}
*/


template <typename T>
T findmin(T *array, long len, long &i){
	long km, k;
	T m=array[0];
	T mp=array[0];
	km=0;
		#pragma omp parallel private(mp, km) if(len>parallel_threshold)
{
		mp = array[0];
		km=0;
		#pragma omp for
		for (k=0; k<len; k++){
				 if(array[k]<mp){				
						mp=array[k];
						km = k; 					
				 }			
		}

	#pragma omp flush(m)
	if(mp<m){
		#pragma omp critical
		{
		if (mp<m){
			m=mp;
			i=km;
		}
		}
	}
}
	return m;	
}


/*
int largest = 0;
int lp
#pragma omp parallel private(lp)
{
  lp = 0;
  #pragma omp for
  for ( int i = 0; i < 1000; i++) {
    if ( data[i] > lp ) lp = data[i];
  }
  if ( lp > largest ) {
    #pragma critical
    {
      if ( lp > largest ) largest = lp;
    }
  }
}
*/


//saves 2D array of type T to text file
//transposed by default, so that row vectors are columns in the file
template <typename T>
void savetxt(std::string file_name, T **array, long sizeX, long sizeY/*=1*/, long precision/*=10*/, char transposed/*='y'*/){
	long i, j;
	std::ofstream file2; 
	file2.open(file_name.c_str()); 
	if (transposed != 'y'){
		for(i = 0; i<sizeY; i++){
			for(j = 0; j<sizeX; j++){
			     file2 << std::setprecision(precision)<< array[i][j] <<" ";
			}
			file2 << std::endl;
		} 
	}
	else{
		for(j = 0; j<sizeX; j++){
			for(i = 0; i<sizeY; i++){	
			     file2 << std::setprecision(precision)<< array[i][j] <<" ";
			}
			file2 << std::endl;
		} 
	}
	file2.close();
}


//saves 1D array of type T to text file
//!!treats it as a 2D array
//transposed by default, so that row vectors are columns in the file
template <typename T>
void savetxt(std::string file_name, T *array, long sizeX, long sizeY/*=1*/, long precision/*=10*/, char transposed/*='y'*/){
	long i, j;
	std::ofstream file2; 
	file2.open(file_name.c_str()); 
	if (transposed != 'y'){
		for(i = 0; i<sizeY; i++){
			for(j = 0; j<sizeX; j++){
			     file2 << std::setprecision(precision)<< array[i*sizeX+j] <<" ";
			}
			file2 << std::endl;
		} 
	}
	else{
		for(j = 0; j<sizeX; j++){
			for(i = 0; i<sizeY; i++){	
			     file2 << std::setprecision(precision)<< array[i*sizeX+j] <<" ";
			}
			file2 << std::endl;
		} 
	}
	file2.close();
}


//loads 2D array of type T from text file
//transposed by default, so that it's consistent with savetxt()
template <typename T>
void loadtxt(std::string file_name, T **array, long sizeX, long sizeY/*=1*/, char transposed/*='y'*/){
	long i, j;
	std::ifstream file2;
	file2.open(file_name.c_str()); 
	if (transposed != 'y'){
		for(i = 0; i<sizeY; i++){
			for(j = 0; j<sizeX; j++){
			     file2 >> array[i][j];
			}
		} 
	}
	else{
		for(j = 0; j<sizeX; j++){
			for(i = 0; i<sizeY; i++){
			     file2 >> array[i][j];
			}
		} 
	}
	file2.close();

}


//loads 1D array of type T from text file
//!!treats it as a 2D array
//transposed by default, so that it's consistent with savetxt()
template <typename T>
void loadtxt(std::string file_name, T *array, long sizeX, long sizeY/*=1*/, char transposed/*='y'*/){
	long i, j;
	std::ifstream file2;
	file2.open(file_name.c_str()); 
	if (transposed != 'y'){
		for(i = 0; i<sizeY; i++){
			for(j = 0; j<sizeX; j++){
			     file2 >> array[i*sizeX+j];
			}
		} 
	}
	else{
		for(j = 0; j<sizeX; j++){
			for(i = 0; i<sizeY; i++){
			     file2 >> array[i*sizeX+j];
			}
		} 
	}
	file2.close();

}



//append 2D array of type T to text file
//transposed by default, so that row vectors are columns in the file
template <typename T>
void appendtxt(std::string file_name, T **array, long sizeX, long sizeY/*=1*/, long precision/*=10*/, char transposed/*='y'*/){
	long i, j;
	std::ofstream file2; 
	file2.open(file_name.c_str(), std::ios::app); 
	if (transposed != 'y'){
		for(i = 0; i<sizeY; i++){
			for(j = 0; j<sizeX; j++){
			     file2 << std::setprecision(precision)<< array[i][j] <<" ";
			}
			file2 << std::endl;
		} 
	}
	else{
		for(j = 0; j<sizeX; j++){
			for(i = 0; i<sizeY; i++){	
			     file2 << std::setprecision(precision)<< array[i][j] <<" ";
			}
			file2 << std::endl;
		} 
	}
	file2.close();
}


//saves 1D array of type T to text file
//!!treats it as a 2D array
//transposed by default, so that row vectors are columns in the file
template <typename T>
void appendtxt(std::string file_name, T *array, long sizeX, long sizeY/*=1*/, long precision/*=10*/, char transposed/*='y'*/){
	long i, j;
	std::ofstream file2; 
	file2.open(file_name.c_str(), std::ios::app); 
	if (transposed != 'y'){
		for(i = 0; i<sizeY; i++){
			for(j = 0; j<sizeX; j++){
			     file2 << std::setprecision(precision)<< array[i*sizeX+j] <<" ";
			}
			file2 << std::endl;
		} 
	}
	else{
		for(j = 0; j<sizeX; j++){
			for(i = 0; i<sizeY; i++){	
			     file2 << std::setprecision(precision)<< array[i*sizeX+j] <<" ";
			}
			file2 << std::endl;
		} 
	}
	file2.close();
}




//returns sum of elements of a 2D array of type T
template <typename T>
T sum_array(T **G, long sizeX, long sizeY/*=1*/){
	long i, j;
	T s=0;
		for (i=0; i<sizeY; i++){
			for (j=0; j<sizeX; j++){
				s+=G[i][j];
			}
		}

	return s;
}


//returns sum of elements of a 1D array of type T
template <typename T>
T sum_array(T *G, long len){
	long i;
	T s=0;
		#pragma omp parallel for reduction(+:s) if(len>parallel_threshold)
		for (i=0; i<len; i++){			
				s+=G[i];

		}

	return s;
}




#endif
