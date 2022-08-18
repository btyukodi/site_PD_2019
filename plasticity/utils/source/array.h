/*
General functions related to 2D array manipulation. 
*/


#ifndef ARRAY_H
#define ARRAY_H


//allocate memory for 2D array
template <typename T> 
T **alloc(long sizeX, long sizeY=1);

//deallocate memory
template <typename T> 
void dealloc(T **array, long sizeX, long sizeY=1);

//fill up array with zeros
template <typename T>
void zeros(T **array, long sizeX, long sizeY=1);


//finds maximum value and position 
template <typename T>
T findmax(T **array, long sizeX, long sizeY, long &i, long &j);

//finds minimum value and position 
template <typename T>
T findmin(T **array, long sizeX, long sizeY, long &i, long &j);

//saves 2D array of type T to text file
//transposed by default, so that row vectors are columns in the file
template <typename T>
void savetxt(std::string file_name, T **array, long sizeX, long sizeY=1, long precision=10, char transposed='y');


//saves 2D array of type T to text file
//transposed by default, so that row vectors are columns in the file
template <typename T>
void savetxt(std::string file_name, T *array, long sizeX, long sizeY=1, long precision=10, char transposed='y');



//append 2D array of type T to text file
//transposed by default, so that row vectors are columns in the file
template <typename T>
void appendtxt(std::string file_name, T **array, long sizeX, long sizeY=1, long precision=10, char transposed='y');


//append 2D array of type T to text file
//transposed by default, so that row vectors are columns in the file
template <typename T>
void appendtxt(std::string file_name, T *array, long sizeX, long sizeY=1, long precision=10, char transposed='y');


//loads 2D array of type T from text file
//transposed by default, so that it's consistent with savetxt()
template <typename T>
void loadtxt(std::string file_name, T **array, long sizeX, long sizeY=1, char transposed='y');

//loads 2D array of type T from text file
//transposed by default, so that it's consistent with savetxt()
template <typename T>
void loadtxt(std::string file_name, T *array, long sizeX, long sizeY=1, char transposed='y');

//returns sum of elements of a 2D array of type T
template <typename T>
T sum_array(T **G, long sizeX, long sizeY=1);


//returns sum of elements of a 2D array of type T
template <typename T>
T sum_array(T *G, long len);

//recquired because of templates
#include"../source/array.cpp"

#endif
