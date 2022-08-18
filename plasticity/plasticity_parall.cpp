/*
V2 : inital thresholds from fast/slow, renewals from inherent
*/

#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<string>
#include <time.h>
#include<omp.h>


#include"utils/source/array_struct.cpp"
#include"utils/source/random.h"
#include"utils/source/histogram.h"
#include"utils/source/avalanches.h"
#include"utils/source/performance.cpp"


//convert int to string
#define str(x) dynamic_cast< std::ostringstream & >((std::ostringstream() << std::dec << x )).str()


//define a few global parameters
//#define size_x 128
//#define size_y 256
#define max_local_yield_stress 1.0	
#define total_number_of_steps 1e9//1e10
#define trans_steps_coeff 50//50//50//1//0//100 //transient = trans_steps_coeff*(long)(float(size*size)/sqrt(d))



//simple linear interpolator assuming x = [0, dx, 2 dx, ...]
// y here will be the resampled inverse CDF
// new value is computed at xx
double interp_simple(double dx, array<double>& y,  double xx){
    long ii = (long)(xx/dx);
    //std::cout<< ii/y.len <<std::endl;
    if (ii< (y.len -1) ){
	
        return y(ii) + (y(ii+1) - y(ii)) * (xx- (ii*dx)) / (dx);
    }
    else{
        return y(y.len-1);
	}

}

//should be moved to utils/random
//rng with distribution p(x)~x^-alpha exp(-lamb*x)
double power_law_cutoff(double xmin, double lamb, double alpha){
	double x, r, p_acc, r_acc;
	p_acc = 0.0;
	r_acc = 1.0;
	while (r_acc > p_acc){
		r = rand(1.0);
		//std::cout<<r<<std::endl;

		x = xmin - 1.0/lamb*log(1-r);

		r_acc = rand(1.0);
		p_acc = pow( (x/xmin), -alpha);
	}
	return x;
}


//generate exponential distribution
double ran_expo(double lambda){
    double u;
    u = rand() / (RAND_MAX + 1.0);
    return -log(1- u) / lambda;
}


//for a given configuration, puts the indices of the weakest spot into i and j, passed by reference
//returns the weakest external stress value SIGMA_C
//this can be optimized by always keeping record of the spots' strength in increasing order
double find_weakest(long &i, long &j, array<double>& SigmaY, array<double>& Sigma){
	i=0;
	j=0;

	double SIGMA_C = SigmaY(0,0) - Sigma(0,0);
	double mp = SigmaY(0,0) - Sigma(0,0);
	long k;
	long kp=0, km=0;
	long len = SigmaY.len;
	
	//each thread finds its private minimum mp and the corresponding private index kp
	#pragma omp parallel private(mp, kp) if(len>parallel_threshold)
	{
		mp = SigmaY(0,0) - Sigma(0,0);
		kp = 0;
		#pragma omp for
		for (k=0; k<len; k++){
				 if(SigmaY(k)-Sigma(k)<mp){				
						mp=SigmaY(k)-Sigma(k);
						kp = k; 					
				 }			
		}

		#pragma omp flush(SIGMA_C, km)
		//compare minimums of each thread
		if(mp<SIGMA_C){
			#pragma omp critical
			{
				if (mp<SIGMA_C){
					SIGMA_C=mp;
					km=kp;
				}
			}
		}
	}
		i=km/SigmaY.sizeX; 
		j=km%SigmaY.sizeX; 
	
	return SIGMA_C;
}



long pbc_shift(long size, long dx){
	long dxret=dx;
	if ( fabs(dx)>0.5*size ){
		if (dx>0){
			dxret = dx - size;//size - i + k; 
		}
		else{
			dxret = dx + size;
		}
	}
	return dxret;
}


/*inline long abs(long x){
	if(x<0){
		x=-x;
	}
return x;
}*/

//updates residual stress Sigma[k][l] on site (k,l) when site (i,j) has slipped with a strain increment of eta given the response matrix G
void update_stress(long i, long j, long k, long l, double eta, array<double>& Sigma, array<double>& G){
	//Sigma[k][l] +=eta*G[size-i+k][size-j+l] ?
	long dx, dy;
	long sizeX=G.sizeX, sizeY=G.sizeY;
	//periodic boundary conditions; to be checked for +-1
	dy = k-i;
	dx = l-j;
	

/*
	if ( fabs(dx)>0.5*sizeX ){
		if (dx>0){
			dx = dx - sizeX;//size - i + k; 
		}
		else{
			dx = dx + sizeX;
		}
	}
	if ( fabs(dy)>0.5*sizeY ){
		if (dy>0){
			dy = dy - sizeY;
		}
		else{
			dy = dy + sizeY;
		}
	}

//wtf??
	if (dx>=sizeX/2){
		dx = -sizeX/2;//= size/2;
	}
	if (dy>=sizeY/2){
		dy = -sizeY/2;//size/2;
	}
*/

	//Sigma.values[k][l]+=eta*G.values[dy + sizeY/2][dx + sizeX/2];
	Sigma(k,l)+=eta*G(dy + sizeY/2, dx + sizeX/2);

}



//update residual stresses on all sites when site (i,j) has slipped
void update_all_stresses(long i, long j, double eta, array<double>& Sigma, array<double>& Ux, array<double>& Uy, array<double>& G, array<double>& G_Ux, array<double>& G_Uy, double alpha){
	long k, l;
	long len = Sigma.len;
	long GsizeX=G.sizeX, GsizeY=G.sizeY;
	long dx, dy;

	#pragma omp parallel for private(k,l, dx, dy) if(len>parallel_threshold)
	for (k=0; k<Sigma.sizeY; k++){
		for (l=0; l<Sigma.sizeX; l++){
			//update_stress(i,j,k,l,eta,Sigma,G);
			dy = k-i;
			dx = l-j;
			Sigma(k,l)+=eta*alpha*G(dy + GsizeY/2, dx + GsizeX/2);
			Ux(k,l)+=eta*G_Ux(dy + GsizeY/2, dx + GsizeX/2);
			Uy(k,l)+=eta*G_Uy(dy + GsizeY/2, dx + GsizeX/2);
		}
	}

}

//accomplishes one plastic event, returns the external stress (before the event)
//returns indices of plastified sites by reference in xpl, ypl
double propag(long &xpl, long &ypl, array<double>& Sigma, array<double>& Ux, array<double>& Uy, array<double>& SigmaY, array<double>& G, array<double>& G_Ux, array<double>& G_Uy, array<double>& Epsilon, double d, double &eta ,double sigmaY_max, double alpha){
	long i, j;
	double SIGMA_C;//, eta;

	//find weakest site
	SIGMA_C = find_weakest(i, j, SigmaY, Sigma);

	//palstify weakest site
	eta = rand(d);
	Epsilon(i,j)+=eta;
	update_all_stresses(i,j,eta,Sigma, Ux, Uy, G, G_Ux, G_Uy, alpha);

	//draw a new yield stress 
	SigmaY(i,j)=rand(sigmaY_max);

	//how exactly do we measure the distance with PBC??? ; ez nem ide kell, hanem az abs(xpl - xplnew) % (size/2)
	ypl = i; // (i % (size/2))
	xpl = j; 


	return SIGMA_C;
}

//find unstable sites; for the moment it's in serial
//array with unstable site indices passed by reference
//returns 0 if no unstable sites were found, 1 otherwise
int find_unstable(double Sigma_ext, array<long>& unstable_indices, array<double>& Sigma, array<double>& SigmaY){
	long len = SigmaY.len;
	long i;
	int flag=0;
	unstable_indices.len=0;
	//unstable_sign.len=0;
	//should be made parallel wisely; simplest way is to set unstable_indices(unstable_indices.len)=i; atomic
	for (i=0; i<len; i++){
		//forward slip
		if (Sigma_ext>=SigmaY(i)-Sigma(i)){
			unstable_indices(unstable_indices.len)=i;
			unstable_indices.len++;
			/*unstable_sign(unstable_sign.len)=0;
			unstable_sign.len++;*/
			flag=1;
//printf("%ld\n", i);
		}
		//backward slip
		/*if (Sigma_ext<=-SigmaY(i)-Sigma(i)){
			unstable_indices(unstable_indices.len)=i;
			unstable_indices.len++;
			unstable_sign(unstable_sign.len)=1;
			unstable_sign.len++;
			flag=1;

		}  */
	}
//printf("%ld\n", unstable_indices.len);
return flag;
}

//another kind of propagator which slips all unstable sites synchronously, until the system becomes stable
//then the applied strain is incremented; strain controlled experiment
double propag_sync(double &ep, double &Sigma_ext_start, double &Sigma_ext_final, double &deps_load, array<long>& unstable_indices, array<double>& Sigma, array<double>& Ux, array<double>& Uy, array<double>& SigmaY, array<double>& G, array<double>& G_Ux, array<double>& G_Uy, array<double>& Epsilon, double d, double &eta ,array<double>& late_dist_y, double dx_late, /*array<int>& Has_slipped,*/ double &eps_load, array<int>& SigmaY_type, double alpha, double mu_shear, double c0, double c1/*, double xmin*/){

int flag=1, nupdates;
int sizeX = Sigma.sizeX;
int sizeY = Sigma.sizeY;
long i, j, k;
double r;
//should mark the start of an avalanche
//int avalanche_start=0;

int flip_nr=0;
double Sigma_ext;
//load until the first event occurs
//possible issue from numerical accuracy: Sigma_ext may be slightly smaller than SigmaY(i)-Sigma(i); add a small number to it
Sigma_ext = find_weakest(i, j, SigmaY, Sigma);
//save the loading strain increment
deps_load =1.0/(2*mu_shear) *(Sigma_ext - Sigma_ext_final);
eps_load+=deps_load;
//start of an avalanche
nupdates=0;
//save stress at the beginning of an avalanche
Sigma_ext_start = Sigma_ext;




while (flag){
	//printf("%lf\n", Sigma_ext);
	flip_nr=0;
	//only flip one in the first round; the stress is tuned to flip this one
	if (nupdates==0){
		unstable_indices.len=0; //not necessary
		unstable_indices(unstable_indices.len)=i*sizeX+j;
		unstable_indices.len++;
		//start of the avalanche at the first update round
		//avalanche_start=1;
	}
	//for further rounds, flip synchronously
	else{
		flag = find_unstable(Sigma_ext, unstable_indices, Sigma, SigmaY);
		//avalanche_start=0;
	}
	//flip unstable sites one by one

	for(k=0; k<unstable_indices.len; k++){		
		//eta = ran_expo(1.0/d);//d;//power_law_cutoff(xmin, 1.0/0.026, 1.0);//d;//!!!!!!!!!!!!!!!!!!!!!!rand(d);



	//	printf("%ld %ld\n", unstable_indices(k),unstable_indices.len);
		i=unstable_indices(k)/SigmaY.sizeX; 
		j=unstable_indices(k)%SigmaY.sizeX; 
		//d is not relevant here
		eta = ran_expo(1.0/(c1*SigmaY(i,j)*SigmaY(i,j) +c0 ) );
		//printf("---- %lf\n", eta);
		//printf("xx---- %d %d %d\n", c0, c1, (c1*SigmaY(i,j)*SigmaY(i,j) +c0 ) );

		Epsilon(i,j)+=eta;
		ep+=eta/sizeX/sizeY;
		update_all_stresses(i,j,eta,Sigma, Ux, Uy, G, G_Ux, G_Uy, alpha);
		
		

		//draw a new yield stress 
		//here will need the condition based on whether it's soft site, has slipped yet, etc.

		if (SigmaY_type(i,j)==0){
			r = rand(1.0);
			SigmaY(i,j)= interp_simple(dx_late, late_dist_y, r);
			/*if (late_yield=="slow"){
				SigmaY(i,j)= interp_simple(dx_slow, slow_dist_y, r);
			}
			if (late_yield=="fast"){
				SigmaY(i,j)=interp_simple(dx_fast, fast_dist_y, r);
			}	
			if (late_yield=="inherent"){
				SigmaY(i,j)=interp_simple(dx_inherent, inherent_dist_y, r);
			}*/	
		}
		//else no need to redraw	


		flip_nr++;

		//in Botond's picture, eps_load = eps_tot; eps_el = eps_tot - ep =(1/mu)*Sigma_ext
		Sigma_ext-=2.0*mu_shear*eta/sizeX/sizeY;/* x shear modulus/sizeX^2?*/

	}

//printf("%lf\n", Sigma_ext);
nupdates++;
}


//save stress at the end of an avalanche
Sigma_ext_final = Sigma_ext;

//save stress at the end of an avalanche
//if there was at least one update round
/*if (nupdates>0){
avalanche_nr++;
}
//printf("***---***\n");*/
}


int main(int argc, char *argv[]){

const long sizeX=atol(argv[1]);
const long sizeY=atol(argv[2]);

double concentration = atof(argv[10]);
double sigmaY_inclusions = atof(argv[11]);

//stress prefactor tuning parameter
double alpha = atof(argv[12]);

//fast renewal distribution shift parameter
//double dm = atof(argv[13]);

//get ensemble
int ens_id=atoi(argv[5]);

//indicates initial yield stress distribution; either slow or fast corresponding to Sylvain's slow and fast quench
std::string init_yield = argv[8];
std::string late_yield = argv[9];

//seed RNG from the ensemble id
srand(100*ens_id);

//-----------------------------------------
//Setting up the slow and fast quench yield stress distributions
array<double> init_dist(("F_inv_"+init_yield+"_resampled.dat").c_str(), 100000, 2);
array<double> init_dist_y(100000);
array<double> late_dist(("F_inv_"+late_yield+"_resampled.dat").c_str(), 100000, 2);
array<double> late_dist_y(100000);
/*array<double> inherent_dist("F_inv_inherent_resampled.dat", 100000, 2);
array<double> inherent_dist_y(100000);*/

double dx_init = init_dist(0, 1) - init_dist(0,0);
double dx_late = late_dist(0, 1) - late_dist(0,0);
//double dx_inherent = inherent_dist(0, 1) - inherent_dist(0,0);

long i,j;
double r;
//double r;
for (i=0; i<init_dist.sizeX; i++){
	init_dist_y(i) = init_dist(1, i);
	late_dist_y(i) = late_dist(1, i);
	//inherent_dist_y(i) = inherent_dist(1,i);
}
//std::cout<<("F_inv_"+init_yield+"_resampled.dat").c_str()<<std::endl;
//-----------------------------------------
//track if site has already slipped
array<int> Has_slipped(sizeX, sizeY);
Has_slipped=0;

//initialize yield stresses (local stress thresholds)
array<double> SigmaY(sizeX, sizeY);
for (i=0; i<sizeY; i++){
	for(j=0; j<sizeX; j++){
		r = rand(1.0);
		SigmaY(i,j)=interp_simple(dx_init, init_dist_y, r);
		/*if (init_yield=="slow"){
			SigmaY(i,j)=interp_simple(dx_slow, slow_dist_y, r);
		}
		if (init_yield=="fast"){
			SigmaY(i,j)=interp_simple(dx_fast, fast_dist_y, r);
		}
		if (init_yield=="inherent"){
			SigmaY(i,j)=interp_simple(dx_inherent, inherent_dist_y, r);
		}*/		
	}
}
//-----------------------------------------
//spread soft sites here
//type is 1 for soft site, 0 for normal site
array<int> SigmaY_type(sizeX, sizeY);
SigmaY_type=0;
for (i=0; i<sizeY; i++){
	for(j=0; j<sizeX; j++){
		if (rand(1.0)<concentration){
			SigmaY(i,j) = sigmaY_inclusions;
			SigmaY_type(i,j) = 1;
			
		}
	}
}


//---------------------

//file name for the Green function from the first argument in the command line
std::string green_file_name(argv[3]);

std::string green_file_name_ux(argv[6]);
std::string green_file_name_uy(argv[7]);

//directory path to the green function
std::string green_path="green_functions/";
std::string green_path_ux="displ_green_functions/";
std::string green_path_uy="displ_green_functions/";
green_path+=green_file_name;
green_path_ux+=green_file_name_ux;
green_path_uy+=green_file_name_uy;

std::string directory_command="mkdir -p ";

//define a suffix; it will be added to each data file
std::string file_suff = "_"+green_file_name+"_d_"+argv[4]+"_inityield_"+argv[8]+"_lateyield_"+argv[9]+"_conc_"+argv[10]+"_sy_"+argv[11]+"_alpha_"+argv[12]+"_mu_"+argv[13]+"_c0_"+argv[14]+"_c1_"+argv[15]/*+"_xmin_"+argv[14]*/;

//define folder for data
std::string data_folder ="data/data";
data_folder+=argv[5];
data_folder+="/";
printf("%s\n", data_folder.c_str());

logfile = data_folder+"log/logfile"+file_suff+".log";

//==================================================================================================================//
//===================================== create directory path ======================================================//
//==================================================================================================================//
//remove previous logfile
log();	
log("Creating directory list...\n");
log("-", system( (directory_command + "data").c_str()));
log("-",system( (directory_command + data_folder).c_str()));
//log("-",system( (directory_command + data_folder+"avalanche").c_str()));
//log("-",system( (directory_command + data_folder+"avalanche_eps").c_str()));
log("-",system( (directory_command + data_folder+"log").c_str()));
//log("-",system( (directory_command + data_folder+"profiles").c_str()));
log("-",system( (directory_command + data_folder+"strain_map").c_str()));
log("-",system( (directory_command + data_folder+"stress_map").c_str()));
//log("-",system( (directory_command + data_folder+"stress_statistics").c_str()));
//log("-",system( (directory_command + data_folder+"stress_statistics/mean").c_str()));
//log("-",system( (directory_command + data_folder+"stress_statistics/sigma").c_str()));
//log("-",system( (directory_command + data_folder+"stress_statistics/Sigma").c_str()));
//log("-",system( (directory_command + data_folder+"stress_statistics/sigma_c").c_str()));
//log("-",system( (directory_command + data_folder+"strain_diff_PSD").c_str()));
//log("-",system( (directory_command + data_folder+"stress_statistics/SigmaY").c_str()));
//log("-",system( (directory_command + data_folder+"stress_statistics/w").c_str()));
log("-",system( (directory_command + data_folder+"stress_strain").c_str()));
//log("-",system( (directory_command + data_folder+"strain_variance").c_str()));
//log("-",system( (directory_command + data_folder+"stress_statistics/sigma_wyart").c_str()));
log("-",system( (directory_command + data_folder+"slipping_site").c_str()));
//log("-",system( (directory_command + data_folder+"diffusion").c_str()));
log("-",system( (directory_command + data_folder+"displacements").c_str()));
log("-",system( (directory_command + data_folder+"stress_variance").c_str()));



//==================================================================================================================//
//============================================ initialize playground ===============================================//
//==================================================================================================================//



//get maximum strain increment from second argument in the command line
double d=atof(argv[4]);

double mu_shear=atof(argv[13]);

//minimum of the plastic strain distribution (if a power law)
//double xmin = atof(argv[14]);

//add two more parameters for the threshold - slip increment correlation
//<ep> ~ c0 * c1*sy^2

double c0 = atof(argv[14]);
double c1 = atof(argv[15]);

//printf("init: %lf %lf", c0, c1);
//double delt = atof(argv[4]);

//current strain increment eta = rand(d)
double eta;




//set maximum locald stress treshold; 1.0 by default
double sigmaY_max=max_local_yield_stress;
long pi;



//SigmaY=1.0;
/*SigmaY.randomize();
SigmaY*=sigmaY_max;*/

//actual residual stresses
array<double> Sigma(sizeX, sizeY);
Sigma=0.0;


//stress barriers (dSigma = SigmaY - Sigma)
array<double> dSigma(sizeX, sizeY);
dSigma=SigmaY;
dSigma-=Sigma;

//plastic strain
array<double> Epsilon(sizeX, sizeY);
Epsilon=0.0;



//Sigma^2 for the variance
array<double> Sigma2(sizeX, sizeY);

//Epsilon^2 for the variance
array<double> Epsilon2(sizeX, sizeY);
//temp for the PR
array<double> Epsilon_tmp(sizeX, sizeY);
array<double> Sigma_tmp(sizeX, sizeY);


//interaction matrix; initialize from file 
//size 3 times larger than the system; rewrite it if out of memory
//G should be the real Green's function / mu_shear
array<double> G(green_path, 3*sizeX, 3*sizeY);

//x - y displacement Green functions
array<double> G_Ux(green_path_ux, 3*sizeX, 3*sizeY);
array<double> G_Uy(green_path_uy, 3*sizeX, 3*sizeY);

array<double> Ux(sizeX, sizeY);
Ux = 0;

array<double> Uy(sizeX, sizeY);
Uy = 0;

array<double> Ux2(sizeX, sizeY);
Ux2 = 0;

array<double> Uy2(sizeX, sizeY);
Uy2 = 0;


//auxiliary arrays for previous saves
array<double> Epsilon_prev(sizeX, sizeY);
Epsilon_prev=0.0;

array<double> Sigma_prev(sizeX, sizeY);
Sigma_prev=0.0;

array<double> SigmaY_prev(sizeX, sizeY);
SigmaY_prev=0.0;

double eps_load_prev;

//------- save parameters --------------//
//remove previous logfile		//
log();					//
log("----System parameters----\n");	//
log(green_file_name);			//
log("sizeX", sizeX);			//
log("sizeY", sizeY);			//
log("d",d);				//
log("SigmaY_inclusions", sigmaY_inclusions);
log("concentration", concentration);
log("initial yield distribution", init_yield);
log("late yield distribution", late_yield);
//--------------------------------------//



//-------------------------- synchronous dynamics with incremental strain---------------

double Sigma_ext_start, Sigma_ext_final;
double deps_load = 0.0;
double eps_load=0.0;

double deps_load_incr=1e-4;

//reserve memory but initialize its length to 0
array<long> unstable_indices(sizeX*sizeY);
unstable_indices.len = 0;

array<int> unstable_sign(sizeX*sizeY);
unstable_sign.len = 0;
//double mu_shear = 0.2; //get it from argv[]
//adimensionalize GF; not necessary when G/mu_shear is feeded:
//G/=mu_shear;

//total average plastic strain
double ep=0.0;

//auxiliary parameter for data recording
double ep_last=0.0;

long avalanche_nr=0;

double sample_window = 0.005;
double max_sampled_up_to_now = 0.0;

double Sigma_ext_prev;
				//plastic strain variance, total strain variance
double dSigma2, dEps2, dEps2_total, PR1, PR2, Eps_tot_avg;



//---------------------- save parameters -----------------------------------------------//
int nthreads=1;										//
#pragma omp parallel if(sizeX*sizeY>parallel_threshold)					//
{											//
nthreads = omp_get_num_threads();							//
}											//
log("==============================================================================\n");//
log("\n" + currentDateTime()+"\t run started; number of threads:", (long)nthreads);	//
//--------------------------------------------------------------------------------------//


//remove previous files
if (ens_id==1){
Sigma.save(data_folder+"stress_map/Sigma_init"+file_suff);      //
SigmaY.save(data_folder+"stress_map/SigmaY_init"+file_suff);      //

Epsilon.save(data_folder+"strain_map/Epsilon_plastic"+file_suff);
}
//Ux.save(data_folder+"displacements/Ux"+file_suff);
//Uy.save(data_folder+"displacements/Uy"+file_suff);
//Sigma.save(data_folder+"stress_map/Sigma"+file_suff);


//======================================================================================//
//======================= Simulation starts here =======================================//
//======================================================================================//


//------------------------------- transient --------------------------------------------//
	std::ofstream stress_strain_file; 
	stress_strain_file.open((data_folder+"stress_strain/stress_strain_curve_"+file_suff).c_str()); 

	std::ofstream stress_strain_file_old; 
	stress_strain_file_old.open((data_folder+"stress_strain/stress_strain_curve_old_"+file_suff).c_str()); 

	
	std::ofstream stress_variance_file; 
	stress_variance_file.open((data_folder+"stress_variance/stress_variance_"+file_suff).c_str()); 

	std::ofstream slipping_site_file; 
	slipping_site_file.open((data_folder+"slipping_site/slipping_site_"+file_suff).c_str()); 
	
//for (i=0; i<0*transient; i++){
//from Mehdi's paper, rule of thumb: steady state reached when <ep>/d^0.5 = 10.0; put 50.0 to be sure
//while (ep<50.0*sqrt(d)){
	//add 0,0 first
stress_strain_file << std::setprecision(10)<< 0.0 <<"\t"<<0.0<<"\t"<<0.0<<"\t"<< 0.0 <<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<std::endl;
while (eps_load<10.0){
i++;
//	sigma = propag(xplnew, yplnew, Sigma, Ux, Uy, SigmaY, G, G_Ux, G_Uy, Epsilon, d, eta, sigmaY_max);

Sigma_prev = Sigma;
Epsilon_prev = Epsilon;
SigmaY_prev = SigmaY;
Sigma_ext_prev = Sigma_ext_final;
eps_load_prev = eps_load;

propag_sync(ep, Sigma_ext_start, Sigma_ext_final, deps_load, unstable_indices, Sigma, Ux, Uy, SigmaY, G, G_Ux, G_Uy, Epsilon, d, eta, late_dist_y, dx_late,eps_load, SigmaY_type, alpha, mu_shear, c0, c1);//, xmin);


stress_strain_file_old << std::setprecision(10)<< eps_load <<"\t"<<Sigma_ext_start<<std::endl;
stress_strain_file_old << std::setprecision(10)<< eps_load <<"\t"<<Sigma_ext_final<<std::endl;


//stress_strain_file << std::setprecision(10)<< eps_load <<"\t"<<Sigma_ext_final<<std::endl;

/*Sigma_ext+=2.0*mu_shear*deps_load_incr;
eps_load+=deps_load_incr;*/


//---------------------------------------- save data -------------------------------------------// 
	//here should save at constant strain window, not steps...
	//should compute strain variance?...; ehh, will compute after from the maps
	//however, quantities saved here should be saved before the propagator!

//if the loading exceeded the next window limit
if (eps_load>max_sampled_up_to_now){
	//it is possible to exceed several window limits
	while (max_sampled_up_to_now < eps_load){
		//!! ONLY QUANTITIES THAT DO NOT CHANGE ON THE ELASTIC BRANCH SHOULD BE SAVED HERE
		if (ens_id<=250){
		Epsilon_prev.append(data_folder+"strain_map/Epsilon_plastic"+file_suff);
		Sigma_prev.append(data_folder+"stress_map/Sigma"+file_suff);
		SigmaY_prev.append(data_folder+"stress_map/SigmaY"+file_suff);
		}
		Sigma2=Sigma_prev;
		Sigma2*=Sigma_prev;
		dSigma2 = Sigma2.avg() - Sigma_prev.avg()*Sigma_prev.avg();

		Epsilon2 = Epsilon_prev;
		Epsilon2*=Epsilon_prev;
		dEps2 = Epsilon2.avg() - Epsilon_prev.avg()*Epsilon_prev.avg();

		//first type of PR (including homogeneous strain)
		Epsilon_tmp = Epsilon_prev;
		Sigma_tmp = Sigma_prev;
		Sigma_tmp*=1.0/(2.0*mu_shear);
		Epsilon_tmp +=Sigma_tmp;
		Epsilon_tmp +=(Sigma_ext_prev + 2.0*mu_shear*(max_sampled_up_to_now-eps_load_prev)) / (2.0*mu_shear);
		Epsilon_tmp*=Epsilon_tmp; //e^2

		PR1 = Epsilon_tmp.sum()*Epsilon_tmp.sum() / (sizeX*sizeY);
		Epsilon_tmp*=Epsilon_tmp; //e^4
		PR1 /= Epsilon_tmp.sum();		


		//second type of PR (not including homogeneous strain)
		Epsilon_tmp = Epsilon_prev;
		Sigma_tmp = Sigma_prev;
		Sigma_tmp*=1.0/(2.0*mu_shear);
		Epsilon_tmp +=Sigma_tmp;
		Epsilon_tmp*=Epsilon_tmp; //e^2

		PR2 = Epsilon_tmp.sum()*Epsilon_tmp.sum() / (sizeX*sizeY);
		Epsilon_tmp*=Epsilon_tmp; //e^4
		PR2 /= Epsilon_tmp.sum();	


		//total strain variance
		Epsilon_tmp = Epsilon_prev;
		Sigma_tmp = Sigma_prev;
		Sigma_tmp*=1.0/(2.0*mu_shear);
		Epsilon_tmp +=Sigma_tmp;
		Eps_tot_avg = Epsilon_tmp.avg();
		Epsilon_tmp*=Epsilon_tmp;
		dEps2_total = Epsilon_tmp.avg() - Eps_tot_avg*Eps_tot_avg;


		stress_strain_file << std::setprecision(10)<< max_sampled_up_to_now <<"\t"<<Sigma_ext_prev + 2.0*mu_shear*(max_sampled_up_to_now-eps_load_prev)<<"\t"<<dEps2<<"\t"<<dSigma2<<"\t"<<PR1<<"\t"<<PR2<<"\t"<<dEps2_total<<std::endl;

		max_sampled_up_to_now+=sample_window;
	}

}

/*
	if (i%100==0){
		Epsilon.append(data_folder+"strain_map/Epsilon_plastic"+file_suff);
		//Ux.append(data_folder+"displacements/Ux"+file_suff);
		//Uy.append(data_folder+"displacements/Uy"+file_suff);
		Sigma.append(data_folder+"stress_map/Sigma"+file_suff);

		SigmaY.append(data_folder+"stress_map/SigmaY"+file_suff);



		Sigma2=Sigma;
		Sigma2*=Sigma;

		//strain_variance_file << std::setprecision(10)<< eps <<"\t"<< Epsilon2.avg() - eps*eps <<std::endl;	//
		stress_variance_file << std::setprecision(10)<< eps_load <<"\t"<< Sigma2.avg() - Sigma.avg()*Sigma.avg() <<std::endl;	

	}	
*/	
												//
	if(i%1000==0){									//
		log("\n" + currentDateTime()+"\t transition? step #", i);   			//
		log("\n plastic strain (eps): ", ep);						//
		log("\n applied loading strain (eps_load): ", eps_load);			
		//
	}											//
//----------------------------------------------------------------------------------------------//			
}

//---------------------------- end of transient ------------------------------------//


stress_variance_file.close();
slipping_site_file.close();
stress_strain_file.close();
log("\n=================== DONE =========================");
}
