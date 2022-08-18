/*
Function for performance test and logging. 
Usage:

	timestamp_t t0 = get_timestamp();

	//piece of code the runtime we are interested in

	timestamp_t t1 = get_timestamp();
	printf("%ld\n", t1-t0);

*/


#include <sys/time.h>
#include <ctime>
#include<fstream>

#ifndef PERFORMANCE_H
#define PERFORMANCE_H


typedef unsigned long long timestamp_t;

//returns timestamp in microseconds
static timestamp_t get_timestamp();


typedef unsigned long long timestamp_t;

    static timestamp_t
    get_timestamp ()
    {
      struct timeval now;
      gettimeofday (&now, NULL);
      return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
    }

//global file for logging
std::string logfile;

//function for timestamp in human readable format
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d \t %X", &tstruct);

    return buf;
}



//function for logging
void log(){
	std::ofstream file2; 
	file2.open(logfile.c_str());
	file2.close();

}

void log(std::string text){
	std::ofstream file2; 
	file2.open(logfile.c_str(),std::ios::app); 
	file2 << text <<" " << std::setprecision(10);
	file2 << std::endl;
	file2.close();
}



void log(std::string text,  double data){
	std::ofstream file2; 
	file2.open(logfile.c_str(),std::ios::app); 
	file2 << text <<" " << std::setprecision(10)<< data <<" ";
	file2 << std::endl;
	file2.close();
}

void log(std::string text,  std::string data){
	std::ofstream file2; 
	file2.open(logfile.c_str(),std::ios::app); 
	file2 << text <<" " << data <<" ";
	file2 << std::endl;
	file2.close();
}

//function for logging
void log(std::string text,  long data){
	std::ofstream file2; 
	file2.open(logfile.c_str(),std::ios::app); 
	file2 << text <<" " << std::setprecision(10)<< data <<" ";
	file2 << std::endl;
	file2.close();
}

void log(std::string text,  int data){
	std::ofstream file2; 
	file2.open(logfile.c_str(),std::ios::app); 
	file2 << text <<" " << std::setprecision(10)<< data <<" ";
	file2 << std::endl;
	file2.close();
}



#endif
