// ConsoleApplication1.cpp
// parallel matrix multiplier using many nodes

/** includes **/
#include "stdafx.h" // visual studio only
#include <iostream>
#include <mpi.h>
#include <stdlib.h>   
#include <cstdlib>
#include <cmath>
#include <cstdio>


int world_rank;
int **arr;


double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

// template <size_t sizeRow, size_t sizeColumn>
void initIntMatrix(int** arr, int matrix_size)
{
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			arr[i][j] = (int) fRand(0,10);
				
		}
	}
}

#pragma warning (disable : 4996) // visual studio does not compile fopen because of security  
template <typename Type>
void printMatrixFromFile(Type** arr, int size)
{
	double* buffer;
	long lSize;
	size_t result;

	// FILE *input;
	// fopen_s(&input,"matrix.dat", "rb");

	FILE *input = fopen("matrix.txt", "rb");
	if (!input) 
	{
		return;
	}

	// obtain file size:
	fseek(input, 0, SEEK_END);
	lSize = ftell(input);
	rewind(input);

	// allocate memory to contain the whole file:
	buffer = (double*) malloc(sizeof(double)*lSize);

	std::cout << "printMatrixFromFile sizeof(double) " << sizeof(double) << " ";
	std::cout << "printMatrixFromFile size*size " << (size*size) << " " << std::endl;


	result = fread_s(buffer, lSize, sizeof(double), size*size, input);
	
	fclose(input);
	for (int i = 0; i < size*size; i++)
	{
		std::cout << buffer[i] << " ";
		if( (i+1) % size == 0)
		{
			std::cout <<  std::endl;
		}	
	}

}

// print a 2 dimensional array
template <typename Type>
void printMatrix(Type** arr, int n, int m)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			std::cout << arr[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

// returns the dot product of two int matricies
int dotProduct(int* arrayA, int** arrayB, int column, int matrix_size)
{

	int result = 0;

	for (int i = 0; i < matrix_size; i++)
	{
		result += arrayA[i] * arrayB[i][column];
	}
	return result;
}


int* multiplyStripe(int* arrayA, int** arrayB, int matrix_size)
{
	int* array_stripe_C = new int[matrix_size];

	for(int i = 0; i < matrix_size; i++)
	{
		array_stripe_C[i] = dotProduct(arrayA, arrayB, i, matrix_size);
	}


	return array_stripe_C;
}

int malloc2dint(int ***array, int n, int m) {

	/* allocate the n*m contiguous items */
	int *p = (int *)malloc(n*m*sizeof(int));
	if (!p) return -1;

	/* allocate the row pointers into the memory */
	(*array) = (int **)malloc(n*sizeof(int*));
	if (!(*array)) {
		free(p);
		return -1;
	}

	/* set up the pointers into the contiguous memory */
	for (int i = 0; i<n; i++)
		(*array)[i] = &(p[i*m]);

	return 0;

	
}

int free2dint(int ***array) 
{
	/* free the memory - the first element of the array is at the start */
	free(&((*array)[0][0]));

	/* free the pointers into the memory */
	free(*array);

	return 0;
}


// function that generates a matrix in binary form on disk
template <typename Type>
void generateMatrixFile(Type** arr, int matrix_size) 
{
	// generate a matrix of values that need to be written to disk in the form of a one dimensional array this will write out an 8x8 matrix
	/*
	double matrix[64] = 
	{
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1
	};
	*/

	// convert to 1d array
	double* p = new double[matrix_size*matrix_size];

	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			p[i * matrix_size + matrix_size] = arr[i][j];
			std::cout << "p[i * matrix_size + matrix_size] " << p[i * matrix_size + matrix_size] << " " << std::endl;
		}
	}

	// open up a file for writing
	FILE *output;
	fopen_s(&output, "matrix.txt", "wb");
	
	/*
	FILE *output = fopen("matrix.dat", "wb");
	*/

	if (!output) 
	{
		return;
	}


	std::cout << "generateMatrixFile sizeof(Type) " << sizeof(double) << " ";
	std::cout << "generateMatrixFile matrix_size*matrix_size " << matrix_size*matrix_size << " " << std::endl;
	// do a simple fwrite to write the matrix to file
	fwrite(p, sizeof(double), matrix_size*matrix_size, output);


	// close the file when we are finished writing
	fclose(output);

	// printMatrix(p);

}


/* master node method */
void coordinator(int world_size, int matrix_size)
{

	// A complete
	arr = new int*[matrix_size];
	for (int i = 0; i < matrix_size; i++)
	{
		arr[i] = new int[matrix_size];
	}

	initIntMatrix(arr, matrix_size);
	/*printMatrix(arr, matrix_size, matrix_size);
*/

	// broadcast size of an individual row 
	MPI_Bcast(&matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	

	int **arrB;
	malloc2dint(&arrB, matrix_size, matrix_size);
	initIntMatrix(arrB, matrix_size);
	// Broadcast B
	MPI_Bcast(&(arrB[0][0]), matrix_size*matrix_size, MPI_INT, 0, MPI_COMM_WORLD);

	/*std::cout << "node: " << world_rank << std::endl;
	printMatrix(arrB, matrix_size, matrix_size);
	//printMatrix(arr);
	std::cout << std::endl;*/


	int size_for_each_node = matrix_size / world_size;
	// broadcast the partition size / stripe size
	MPI_Bcast(&size_for_each_node, 1, MPI_INT, 0, MPI_COMM_WORLD);

	/*int **arrA;
	malloc2dint(&arrA, size_for_each_node, matrix_size);*/

	int ** arrC;
	malloc2dint(&arrC, size_for_each_node, matrix_size);


	int** matrix = new int *[matrix_size];
	int *data = new int[matrix_size*matrix_size];
	for (int i = 0; i < matrix_size; i++)
		matrix[i] = &(data[matrix_size*i]);

	initIntMatrix(matrix, matrix_size);

	printMatrix(matrix, matrix_size, matrix_size);

	/*int** arrA = new int *[matrix_size];
	int *data2 = new int[matrix_size*matrix_size];
	for (int i = 0; i < matrix_size; i++)
		matrix[i] = &(data[matrix_size*i]);*/

	int *partition = new int[size_for_each_node*matrix_size];

	// scatter A
	// scatter the partition to each node
	 MPI_Scatter(matrix, matrix_size * size_for_each_node, MPI_INT, partition, matrix_size * size_for_each_node, MPI_INT, 0, MPI_COMM_WORLD);

	// printMatrix(arrA, size_for_each_node, matrix_size);
	
	
	
/*

	// print A and B
	printMatrix(arr);
	printMatrix(arr2);
	
	// results from stripes
	int* result = new int[matrix_size];
	 
	for (int i = 0; i < matrix_size; i++)
	{
		result = multiplyStripe(arr[i], arr2);	
		for (int j = 0; j < matrix_size; j++)
		{
			arr3[i][j] = result[j];
		}
	}
	
	// print the result
	 printMatrix(arr3);*/
	// generateMatrixFile(arr);

	// printMatrixFromFile(arr, matrix_size);
/*
	free2dint(&arrA);
	free2dint(&arrB);
	free2dint(&arrC);*/

	
}

void print1darrayAs2d(int n, int m, int* arr)
{
	for (int i = 0; i < m*n; i++)
	{

		std::cout << " " << arr[i];
		if((i+1)% n == 0)
		{
			std::cout << std::endl;
		}
	}
}

/* slave node method */
void participant()
{	

	int matrix_size = 0;
	
	// broadcast the matrix size
	MPI_Bcast(&matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
/*	std::cout << "matrix_size: " << matrix_size << std::endl;*/

	int **arrB;
	malloc2dint(&arrB, matrix_size, matrix_size);

	// broadcast B
	MPI_Bcast(&(arrB[0][0]), matrix_size*matrix_size, MPI_INT, 0, MPI_COMM_WORLD);

	/*
	std::cout << "node: " << world_rank << std::endl;
	printMatrix(arrB, matrix_size, matrix_size);
	//printMatrix(arr);
	std::cout << std::endl;
	*/


	int size_for_each_node = 0;
	// broadcast the partition size / stripe size
	MPI_Bcast(&size_for_each_node, 1, MPI_INT, 0, MPI_COMM_WORLD);
	/*std::cout << " " << size_for_each_node << std::endl;
*/
	int **arrA;
	malloc2dint(&arrA, size_for_each_node, matrix_size);

	int ** arrC;
	malloc2dint(&arrC, size_for_each_node, matrix_size);

	int ** arr;
	malloc2dint(&arr, matrix_size, matrix_size);


	int** matrix = new int *[matrix_size];
	int *data = new int[matrix_size*matrix_size];
	for (int i = 0; i < matrix_size; i++)
		matrix[i] = &(data[matrix_size*i]);

/*
	int** arrA = new int *[matrix_size];
	int *data2 = new int[matrix_size*matrix_size];
	for (int i = 0; i < matrix_size; i++)
		matrix[i] = &(data[matrix_size*i]);*/

	int *partition = new int[size_for_each_node*matrix_size];

	MPI_Scatter(matrix, matrix_size * size_for_each_node, MPI_INT, partition, matrix_size * size_for_each_node, MPI_INT, 0, MPI_COMM_WORLD);

	print1darrayAs2d(matrix_size, size_for_each_node, partition);

	// printMatrix(arrA, matrix_size, matrix_size);


	// Matrix test1 = new Matrix(arrA);





	/*free2dint(&arrA);
	free2dint(&arrB);
	free2dint(&arrC);*/

	


}

int main(int argc, char** argv) 
{
	MPI_Init(NULL, NULL);
	int world_size;
	int matrix_size = 16;


	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	/*  get the matricies */
	sscanf_s(argv[1], "%d", &matrix_size);

	
	// master node 
	if(world_rank == 0) 
	{
		std::cout << "matrix size: " << matrix_size << std::endl;
		auto start_time = MPI_Wtime();
		coordinator(world_size, matrix_size);
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


