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
			arr[i][j] = 1;
				//(int) fRand(0,50);
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
template <size_t sizeRowA, size_t sizeColumnA, size_t sizeRowB, size_t sizeColumnB>
int dotProduct_two_dimensional(int(&matrixA)[sizeRowA][sizeColumnA], int row, int(&matrixB)[sizeRowB][sizeColumnB], int column)
{
	// only works if the size of the row of A is equal to size of column of B
	if(sizeRowA != sizeColumnB)
	{
		return EXIT_FAILURE;
	}
	int result = 0;
	for (size_t i = 0; i < sizeColumnA; i++)
	{
		result += matrixA[row][i] * matrixB[i][column];
	}
	return result;
}

// returns the dot product of two int matricies
template <size_t sizeRowA, size_t sizeRowB, size_t sizeColumnB>
int dotProduct(int(&arrayA)[sizeRowA], int(&arrayB)[sizeRowB][sizeColumnB], int column)
{
	// only works if the size of the row of A is equal to size of column of B
	if (sizeRowA != sizeColumnB)
	{
		return EXIT_FAILURE;
	}

	int result = 0;

	for (int i = 0; i < sizeRowA; i++)
	{
		result += arrayA[i] * arrayB[i][column];
	}
	return result;
}



template <size_t sizeRowA, size_t sizeRowB, size_t sizeColumnB>
int* multiplyStripe(int(&arrayA)[sizeRowA], int(&arrayB)[sizeRowB][sizeColumnB])
{
	int* array_stripe_C = new int[sizeColumnB];

	for(int i = 0; i < sizeColumnB; i++)
	{
		array_stripe_C[i] = dotProduct(arrayA, arrayB, i);
	}


	return array_stripe_C;
}

 int(&fillarr(int(&arr)[5]))[5]{ // no decay; argument must be size 5
	 return arr;
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

	arr2[0][0] = 5; // testing

	printMatrix(arr);
	printMatrix(arr2);

	//int dotproduct = dotProduct_two_dimensional(arr, 0, arr2, 0);
	// std::cout << dotproduct << " ";


	int arr3[8][8];
	int* result = new int[8];

	for (int i = 0; i < 8; i++)
	{
		result = multiplyStripe(arr[i], arr2);	
		for (int j = 0; j < sizeof(result); j++)
		{
			arr3[i][j] = result[j];
		}
	}

	printMatrix(arr3);
	



	


	
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





