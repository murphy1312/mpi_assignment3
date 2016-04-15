// ConsoleApplication1.cpp
// parallel matrix multiplier using many nodes

/** includes **/
#include "stdafx.h" // visual studio only
#include <iostream>
#include <mpi.h>
#include <stdlib.h>   
#include <cstdlib>
#include <cmath>



int world_rank;

double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

template <size_t sizeRow, size_t sizeColumn>
void initIntMatrix(int(&arr)[sizeRow][sizeColumn])
{
	for (int i = 0; i < sizeRow; i++)
	{
		for (int j = 0; j < sizeColumn; j++)
		{
			arr[i][j] = (int) fRand(0,50);
		}
	}
}

#pragma warning (disable : 4996) // visual studio does not compile fopen because of security  
void printMatrixFromFile() 
{
	double matrix[64];
	FILE *input = fopen("matrix.dat", "rd");
	if (!input) {
		return;
	}
	fread(matrix, sizeof(double), 64, input);
	fclose(input);
	for (int i = 0; i < 64; ++i) {
		if (i % 8 == 0)
			std::cout << std::endl;
		else
			std::cout << matrix[i] << " ";
	}
}

// print a 2 dimensional array
template <typename Type, size_t sizeRow, size_t sizeColumn>
void printMatrix(Type(&arr)[sizeRow][sizeColumn])
{
	for (int i = 0; i < sizeRow; i++) 
	{
		for (int j = 0; j < sizeColumn; j++)
		{
			std::cout << arr[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

// returns the dot product of two int matricies
template <size_t rowsA, size_t colsA, size_t rowsB, size_t colsB>
int dotProduct(int(&matrixA)[rowsA][colsA], int row, int(&matrixB)[rowsB][colsB], int column)
{
	int result = 0;
	for (size_t i = 0; i < colsA; i++)
	{
		result += matrixA[row][i] * matrixB[i][column];
	}
	return result;
}

// function that generates a matrix in binary form on disk
void generateMatrixFile() {
	// generate a matrix of values that need to be written to disk in the form of a one dimensional array this will write out an 8x8 matrix
	double matrix[64] = {
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1
	};

	// open up a file for writing
	/*
	FILE *output;
	fopen_s(&output, "matrix.dat", "wb");
	*/

	FILE *output = fopen("matrix.dat", "wb");
	if (!output) {
		return;
	}

	// do a simple fwrite to write the matrix to file
	fwrite(matrix, sizeof(double), 64, output);


	// close the file when we are finished writing
	fclose(output);

	printMatrixFromFile();

}


/* master node method */
void coordinator(int world_size)
{
	int arr[8][8];
	int arr2[8][8];
	initIntMatrix(arr2);
	initIntMatrix(arr);

	printMatrix(arr);
	printMatrix(arr2);

	int dotproduct = dotProduct(arr, 0, arr2, 0);
	std::cout << dotproduct << " ";




	


	
}

/* slave node method */
void participant()
{
	
}

int main(int argc, char** argv) 
{
	MPI_Init(NULL, NULL);
	int world_size;


	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	/*  get the matricies */
	// sscanf_s(argv[1], "%d", &iterations); 

	
	// master node 
	if(world_rank == 0) 
	{
		auto start_time = MPI_Wtime();
		coordinator(world_size);
		auto end_time = MPI_Wtime();
		// output the time
		std::cout << "time in s:" << end_time-start_time << std::endl;

	}
	// other nodes
	else
	{
		participant();
	}

	MPI_Finalize();
	return 0;
}





