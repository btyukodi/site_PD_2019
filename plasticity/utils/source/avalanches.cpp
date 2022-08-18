#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"array_struct.cpp"
#include"histogram.h"



double spring_stress(double k, double deps, double last_sigma_eps){
	return last_sigma_eps- k*deps;
}


double end_of_avalanche(double k, double deps, double last_sigma_eps ,double sigma){

	double k_stress = spring_stress(k, deps, last_sigma_eps);
	if (k_stress<sigma){
		return deps;
	}	
	else {
		return -1;
	}
}




void avalanche(array<double>& ksprings, array<long>& hist, array<double>& deps, array<double>& last_sigma, double eps, double sigma){
	int nsprings = hist.sizeY;
	long len_h = hist.sizeX;
	double eps_res = hist.resX;
	

	int i;
	double dt;
/*	for (i=0; i<nsprings; i++){
			dt = end_of_avalanche(ksprings(i), deps(i), last_sigma(i), sigma);
			//if any of avalanches reached its end
			if (dt > -1){
				//update_1D_hist(hist.values[0][i], len_h, 0.0, eps_res, deps[i]);
				update_hist(hist,deps(i),i);
				deps(i)=0.0;
				last_sigma(i)=sigma;
				
			}	
			else{
				deps(i)+=eps;
			}

	}
*/
	for (i=0; i<nsprings; i++){
			deps(i)+=eps;
			dt = end_of_avalanche(ksprings(i), deps(i), last_sigma(i), sigma);
			//if any of avalanches reached its end
			if (dt > -1){
				//update_1D_hist(hist.values[0][i], len_h, 0.0, eps_res, deps[i]);
				update_hist(hist,deps(i),i);
				deps(i)=0.0;
				last_sigma(i)=sigma;
				
			}	
			
	}
//	printf("%e\n",eps);

}


