// ConsoleApplication1.cpp
// use monte carlo to calculate pi
/** simple program to test the MPI stuff to see if it works **/
/** includes **/
#include "stdafx.h"
#include <iostream>
#include <mpi.h>
#include <stdlib.h>    
#include <random>
#include <chrono>

int world_rank;

double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}


int calculateHits(int iterations)
{
	double x = 0;
	double y = 0; 
	double z = 0;
	int i = 0;
	int hits = 0;

	// random 
	std::mt19937_64 rng;
	// init the random generator with seed
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count()*world_rank;
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	// initialize a uniform distribution between 0 and 1
	std::uniform_real_distribution<double> unif(-1, 1);

	for (i = 0; i<iterations; i++) 
	{
		
		x = unif(rng);
		y = unif(rng);
		z = x*x + y*y;
		if (z <= 1)
		{
			hits++;
		}
	}
	return hits;
}

/* master node method */
void coordinator(int iterations, int world_size)
{
	double pi;
	int hits = 0;
	int total_hits;

	total_hits = calculateHits(iterations);

	// gather data from other nodes
	for (int i = 1; i < world_size; i++) 
	{
		MPI_Recv(&hits, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		total_hits += hits;
	}

	// calculate pi
	pi = (double) total_hits * 4 / (iterations * world_size);
	
	// output result
	std::cout.precision(15);
	std::cout.setf(std::ios::fixed, std::ios::floatfield); // floatfield set to fixed
	std::cout << "PI ~ "<< pi << std::endl; 
	
}

/* slave node method */
void participant(int iterations)
{
	int hits;
	hits = calculateHits(iterations);

	// send data to master node
	MPI_Send(&hits, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
}

int main(int argc, char** argv) 
{
	MPI_Init(NULL, NULL);
	int iterations = 0;
	int world_size;


	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	/* get the number of iterations*/
	sscanf_s(argv[1], "%d", &iterations); 

	
	// master node 
	if(world_rank == 0) 
	{
		auto start_time = MPI_Wtime();
		coordinator(iterations/world_size, world_size);
		auto end_time = MPI_Wtime();
		// output the time
		std::cout << "time in s:" << end_time-start_time << std::endl;

	}
	// other nodes
	else
	{
		participant(iterations/world_size);
	}

	MPI_Finalize();
	return 0;
}





