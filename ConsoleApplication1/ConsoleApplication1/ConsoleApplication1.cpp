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

template <typename Type>
void printMatrix(Type** arr, int n, int m);
int free2dint(int ***array);
double fRand(double fMin, double fMax);
void fillIntMatrixRnd(int** arr, int matrix_size);
int** readMatrixFromFile(const std::string file, int matrix_size);
int dotProduct(int* arrayA, int** arrayB, int column, int matrix_size);
void coordinator(int world_size, int matrix_size);
int dotProduct(int* arrayA, int** arrayB, int column, int matrix_size);
int* multiplyStripe(int* arrayA, int** arrayB, int matrix_size);
int malloc2dint(int ***array, int n, int m);
void print1darrayAs2d(int n, int m, int* arr);
void participant();



int world_rank;

double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

// template <size_t sizeRow, size_t sizeColumn>
void fillIntMatrixRnd(int** arr, int matrix_size)
{
	for (int i = 0; i < matrix_size; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			arr[i][j] = 1;
				// (int) fRand(0,10);
				
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


int** readMatrixFromFile(const std::string file, int matrix_size)
{
	int** arr = new int *[matrix_size];
	int *dataC = new int[matrix_size*matrix_size];
	for (int i = 0; i < matrix_size; i++)
		arr[i] = &(dataC[matrix_size*i]);

	FILE *input;
	fopen_s(&input, file.c_str(), "rb");
	if (!input) 
	{
		return NULL;
	}
	fread(&(arr[0][0]), sizeof(int), matrix_size*matrix_size, input);
	fclose(input);

	return arr;

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

int malloc2dint(int ***array, int n, int m) 
{

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

void print1darrayAs2d(int n, int m, int* arr)
{
	for (int i = 0; i < m*n; i++)
	{

		std::cout << arr[i] << " ";
		if ((i + 1) % n == 0)
		{
			std::cout << std::endl;
		}
	}
}

/* master node method */
void coordinator(int world_size, int matrix_size, std::string filenameA, std::string filenameB)
{
	// get matrixA and matrixB from files
	int **arrA = readMatrixFromFile(filenameA, matrix_size);
	int **arrB = readMatrixFromFile(filenameB, matrix_size);

	printMatrix(arrA, matrix_size, matrix_size);

	printMatrix(arrB, matrix_size, matrix_size);

	// broadcast size of an individual row 
	MPI_Bcast(&matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Broadcast B
	MPI_Bcast(&(arrB[0][0]), matrix_size*matrix_size, MPI_INT, 0, MPI_COMM_WORLD);

	int size_for_each_node = matrix_size / world_size;
	// broadcast the partition size / stripe size
	MPI_Bcast(&size_for_each_node, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int** arrResultForThisNode = new int *[size_for_each_node];
	int *dataC = new int[size_for_each_node*matrix_size];
	for (int i = 0; i < size_for_each_node; i++)
		arrResultForThisNode[i] = &(dataC[matrix_size*i]);

	int** arrStripeForThisNode = new int *[size_for_each_node];
	int *data2 = new int[size_for_each_node*matrix_size];
	for (int i = 0; i < size_for_each_node; i++)
		arrStripeForThisNode[i] = &(data2[matrix_size*i]);

	// scatter A
	// scatter the partition to each node
	 MPI_Scatter(&(arrA[0][0]), matrix_size * size_for_each_node, MPI_INT, data2, matrix_size * size_for_each_node, MPI_INT, 0, MPI_COMM_WORLD);
	
	 int** overall_result = new int *[matrix_size];
	 int *overall_result_1d = new int[matrix_size*matrix_size];
	 for (int i = 0; i < matrix_size; i++)
		 overall_result[i] = &(overall_result_1d[matrix_size*i]);

	 // caluclation
	 int *result = new int[matrix_size];
	 for (int i = 0; i < size_for_each_node; i++)
	 {
		 result = multiplyStripe(arrStripeForThisNode[i], arrB, matrix_size);
		 for (int j = 0; j < matrix_size; j++)
		 {
			 arrResultForThisNode[i][j] = result[j];
		 }
	 }

	 MPI_Gather(&(arrResultForThisNode[0][0]), matrix_size*size_for_each_node, MPI_INT, overall_result_1d, size_for_each_node*matrix_size, MPI_INT, 0, MPI_COMM_WORLD);
	 printMatrix(overall_result, matrix_size, matrix_size);

	// cleanup
	delete[] overall_result, arrStripeForThisNode, data2, arrResultForThisNode, dataC, overall_result_1d, arrA, arrB, result;
	
}



/* slave node method */
void participant()
{	
	int matrix_size = 0;
	// broadcast the matrix size
	MPI_Bcast(&matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int **arrB;
	malloc2dint(&arrB, matrix_size, matrix_size);

	// broadcast B
	MPI_Bcast(&(arrB[0][0]), matrix_size*matrix_size, MPI_INT, 0, MPI_COMM_WORLD);

	int size_for_each_node = 0;
	// broadcast the partition size (how many rows a node should compute)
	MPI_Bcast(&size_for_each_node, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int** arrStripeForThisNode = new int *[size_for_each_node];
	int *data2 = new int[size_for_each_node*matrix_size];
	for (int i = 0; i < size_for_each_node; i++)
		arrStripeForThisNode[i] = &(data2[matrix_size*i]);

	int** arrResultForThisNode = new int *[size_for_each_node];
	int *dataC = new int[size_for_each_node*matrix_size];
	for (int i = 0; i < size_for_each_node; i++)
		arrResultForThisNode[i] = &(dataC[matrix_size*i]);

	// scatter the value of the rows this node should calculate to partition array
	MPI_Scatter(NULL, matrix_size * size_for_each_node, MPI_INT, data2, matrix_size * size_for_each_node, MPI_INT, 0, MPI_COMM_WORLD);

	// actual caluculation
	int* result = new int[matrix_size];
	for (int i = 0; i < size_for_each_node; i++)
	{
		result = multiplyStripe(arrStripeForThisNode[i], arrB, matrix_size);
		for (int j = 0; j < matrix_size; j++)
		{
			arrResultForThisNode[i][j] = result[j];
		}
	}

	MPI_Gather(&(arrResultForThisNode[0][0]), matrix_size*size_for_each_node, MPI_INT, NULL, size_for_each_node*matrix_size, MPI_INT, 0, MPI_COMM_WORLD);

	// cleanup
	delete[] arrB, arrStripeForThisNode, data2, arrResultForThisNode, dataC, arrB, result;

}


int main(int argc, char** argv) 
{
	MPI_Init(NULL, NULL);
	int world_size;
	int matrix_size = 16;

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	bool shutdown = false;

	// master node 
	if(world_rank == 0) 
	{
		// timer
		auto start_time = MPI_Wtime();

		if (argc <4)
		{
			shutdown = true;
			std::cout << "not all arguments provided [MatrixA filepath] [MatrixB filepath] [Matrixsize]" << std::endl;
		}
		if (argv[1] == NULL)
		{
			shutdown = true;
			std::cout << "Please enter the path to the first matrix" << std::endl;
		}
		if (argv[2] == NULL)
		{
			shutdown = true;
			std::cout << "Please enter the path to the second matrix" << std::endl;			
		}
		if (argv[3] == NULL)
		{
			shutdown = true;
			std::cout << "Please enter the matrix size" << std::endl;
		}
		MPI_Bcast(&shutdown, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
		if (!shutdown)
		{

			//  get the matrix size
			sscanf_s(argv[3], "%d", &matrix_size);
			std::cout << "matrix size: " << matrix_size << std::endl;

			// filenames
			std::string maA(argv[1]);
			std::string maB(argv[2]);

			coordinator(world_size, matrix_size, maA, maB);
			auto end_time = MPI_Wtime();
			// output the time
			std::cout << "time in s:" << end_time - start_time << std::endl;
		}
	}
	// other nodes
	else
	{
		MPI_Bcast(&shutdown, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
		if (!shutdown)
		{
			participant();
		}	
	}

	MPI_Finalize();
	return 0;
}



