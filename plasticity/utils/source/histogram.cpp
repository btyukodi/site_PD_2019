#include<stdio.h>
#include<stdlib.h>
#include"array_struct.cpp"
#include<omp.h>

//add element to 1D histogram
int update_hist(array<long>& h, double new_xvalue){
	double maxX = h.minX + h.sizeX * h.resX;

	//if out of range, drop value; slightly extend the first bin however so that the interval is closed
	if (new_xvalue > maxX || new_xvalue < h.minX-1e-20*h.resX){
		return 0;
	}

	//values along y axis are counted by row index
	h((long)((new_xvalue-h.minX)/h.resX))++;

return 1;
}

//add element to 2D histogram
int update_hist(array<long>& h, double new_xvalue, double new_yvalue){
	double maxX = h.minX + h.sizeX * h.resX;
	double maxY = h.minY + h.sizeY * h.resY;
	//if out of range, drop value
	if (new_xvalue > maxX || new_xvalue < h.minX-1e-20*h.resX || new_yvalue>maxY || new_yvalue<h.minY-1e-20*h.resY){
		return 0;
	}

	//values along y axis are counted by row index
	h((long)((new_yvalue-h.minY)/h.resY),(long)((new_xvalue-h.minX)/h.resX))++;
	
return 1;
}




//use this one for small histograms, i. e. h.len<new_xvalue.len
int update_hist(array<long>& h, array<double>& new_xvalue){
	double maxX = h.minX + h.sizeX * h.resX;
	long i, j;

//version #1: use this if new_xvalue.len > h.len; create individual copies for each thread; memory allocation becomes slower as histogram's size increases
//individual copies for threads may not worth it if the histogram's size is comparable to that of the array's
if(h.len<=new_xvalue.len){
	long * hista;
	#pragma omp parallel if(new_xvalue.len>parallel_threshold)
	{
		const int nthreads = omp_get_num_threads();
	    	const int ithread = omp_get_thread_num();

		#pragma omp single
		hista = new long[h.len*nthreads];

		#pragma omp for private(i)
		for (i=0; i<h.len*nthreads; i++){
		    hista[i]=0;
		}
		//fill up local histograms
		#pragma omp for private(i)
		for (i=0; i<new_xvalue.len; i++){
				//if out of range, drop value
				if (new_xvalue.values[i] > maxX || new_xvalue.values[i] < h.minX-1e-20*h.resX){
					//do nothing
				}
				else{
					//values along y axis are counted by row index (0 in this case)
					hista[h.len*ithread + ((long)((new_xvalue.values[i]-h.minX)/h.resX))]++;

				}

		
		}

		#pragma omp for private(i)
		for(i=0; i<h.len; i++) {
		  for(int t=0; t<nthreads; t++) {
			h(i) += hista[h.len*t+i];
		
		  }
		}

	}
	delete[] hista;
	return 1;
}

// version #2: use this for larger histograms; not as efficient as version #1, but still faster for h.len>new_xvalue.len
// try to keep h.len small
else{
	#pragma omp parallel for private(i) if(new_xvalue.len>parallel_threshold)
	for (i=0; i<new_xvalue.len; i++){
			//if out of range, drop value
			if (new_xvalue.values[i] > maxX || new_xvalue.values[i] < h.minX-1e-20*h.resX){
				//do nothing
			}
			else{
				//values along y axis are counted by row index (0 in this case)
				#pragma omp atomic
				h((long)((new_xvalue.values[i]-h.minX)/h.resX))++;

			}		
	}
return 1;

}

}



//overload the same function to update values from 2 arrays; the dimensions of the 2 arrays should be the same
int update_hist(array<long>& h, array<double>& new_xvalue, array<double>& new_yvalue){
	double maxX = h.minX + h.sizeX * h.resX;	
	double maxY = h.minY + h.sizeY * h.resY;
	long i;
//version #2:
//note that version #1 is unnecessary, for a 2d histogram it would take too long
	#pragma omp parallel for private(i) if(new_xvalue.len>parallel_threshold)
	for (i=0; i<new_xvalue.len; i++){
			//if out of range, drop value
			if (new_xvalue.values[i] > maxX || new_xvalue.values[i] < h.minX-1e-20*h.resX || new_yvalue.values[i]>maxY || new_yvalue.values[i]<h.minY-1e-20*h.resY){
				//do nothing
			}

			else{
			//values along y axis are counted by row index (0 in this case)
				#pragma omp atomic
				h((long)((new_yvalue.values[i]-h.minY)/h.resY),(long)((new_xvalue.values[i]-h.minX)/h.resX))++;
			}
		
	}
return 1;

}

