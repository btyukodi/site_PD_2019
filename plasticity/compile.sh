#!/bin/bash

g++ -c -O3 utils/source/fft.cpp -o utils/obj/fft.o
g++ -c -O3 utils/source/fft2D.c -o utils/obj/fft2D.o
g++ -c -fopenmp -O3 utils/source/histogram.cpp -o utils/obj/histogram.o
g++ -c -O3 utils/source/random.cpp -o utils/obj/random.o
g++ -c -O3 utils/source/avalanches.cpp -o utils/obj/avalanches.o
g++ -fopenmp utils/obj/fft.o utils/obj/fft2D.o utils/obj/histogram.o  utils/obj/random.o utils/obj/avalanches.o  plasticity_parall.cpp -O3 -o plasticity_parall

