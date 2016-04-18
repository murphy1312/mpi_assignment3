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
int **arr;


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
	
	printMatrix(arr, matrix_size, matrix_size);

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

		std::cout << " " << arr[i];
		if ((i + 1) % n == 0)
		{
			std::cout << std::endl;
		}
	}
}

/* master node method */
void coordinator(int world_size, int matrix_size)
{
	// get matrixA and matrixB from files
	int **arrA = readMatrixFromFile("matrix.dat", matrix_size);
	int **arrB = readMatrixFromFile("matrix2.dat", matrix_size);

	printMatrix(arrA, matrix_size, matrix_size);

	printMatrix(arrB, matrix_size, matrix_size);

	// broadcast size of an individual row 
	MPI_Bcast(&matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Broadcast B
	MPI_Bcast(&(arrB[0][0]), matrix_size*matrix_size, MPI_INT, 0, MPI_COMM_WORLD);


	int size_for_each_node = matrix_size / world_size;
	// broadcast the partition size / stripe size
	MPI_Bcast(&size_for_each_node, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int** arrC = new int *[size_for_each_node];
	int *dataC = new int[size_for_each_node*matrix_size];
	for (int i = 0; i < size_for_each_node; i++)
		arrC[i] = &(dataC[matrix_size*i]);


	// matrix[7][7] = 5;

	// original matrix A
	


	// result array for this nodes A
	int** arrResult = new int *[size_for_each_node];
	int *data2 = new int[size_for_each_node*matrix_size];
	for (int i = 0; i < size_for_each_node; i++)
		arrResult[i] = &(data2[matrix_size*i]);


	int *partition = new int[size_for_each_node*matrix_size];


	// scatter A
	// scatter the partition to each node
	 MPI_Scatter(&(arrA[0][0]), matrix_size * size_for_each_node, MPI_INT, partition, matrix_size * size_for_each_node, MPI_INT, 0, MPI_COMM_WORLD);

	 // printMatrix(arrA, size_for_each_node, matrix_size);
	
	// print1darrayAs2d(matrix_size, size_for_each_node, partition);
	
	// 1d array partition into 2d array arrA
	 for (int i = 0; i < size_for_each_node; i++)
	 {
		 for (int j = 0; j < matrix_size; j++)
		 {
			 arrResult[i][j] = partition[j*size_for_each_node + i];
		 }
	 }

	 // print1darrayAs2d(matrix_size, size_for_each_node, partition);
	 // printMatrix(arrA, size_for_each_node, matrix_size);
	

	 int *result = new int[matrix_size];

	 for (int i = 0; i < size_for_each_node; i++)
	 {
		 result = multiplyStripe(arrResult[i], arrB, matrix_size);
		 for (int j = 0; j < matrix_size; j++)
		 {
			 arrC[i][j] = result[j];
		 }
	 }

	// printMatrix(arrC, size_for_each_node, matrix_size);


	 

	 int *matrixC2 = new int[matrix_size*matrix_size];

	 MPI_Gather(&(arrC[0][0]), matrix_size*size_for_each_node, MPI_INT, matrixC2, size_for_each_node*matrix_size, MPI_INT, 0, MPI_COMM_WORLD);


	print1darrayAs2d(matrix_size, matrix_size, matrixC2);

/*
	free2dint(&arrA);
	free2dint(&arrB);
	free2dint(&arrC);*/

	
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
	// std::cout << "size_for_each_node " << size_for_each_node << std::endl;

	// result array for this nodes A
	int** arrA = new int *[size_for_each_node];
	int *data2 = new int[size_for_each_node*matrix_size];
	for (int i = 0; i < size_for_each_node; i++)
		arrA[i] = &(data2[matrix_size*i]);



	int** arrC = new int *[size_for_each_node];
	int *dataC = new int[size_for_each_node*matrix_size];
	for (int i = 0; i < size_for_each_node; i++)
		arrC[i] = &(dataC[matrix_size*i]);


	int **arr;
	malloc2dint(&arr, matrix_size, matrix_size);


	int** matrix = new int *[matrix_size];
	int *data = new int[matrix_size*matrix_size];
	for (int i = 0; i < matrix_size; i++)
		matrix[i] = &(data[matrix_size*i]);


	int *partition = new int[size_for_each_node*matrix_size];

	// scatter the value of the rows this node should calculate to partition array
	MPI_Scatter(&(matrix[0][0]), matrix_size * size_for_each_node, MPI_INT, partition, matrix_size * size_for_each_node, MPI_INT, 0, MPI_COMM_WORLD);
	
	// 1d array partition into 2d array arrA
	for (int i = 0; i < size_for_each_node; i++)
	{
		for (int j = 0; j < matrix_size; j++)
		{
			arrA[i][j] = partition[j*size_for_each_node + i];
		}
	}

	// print1darrayAs2d(matrix_size, size_for_each_node, partition);
	// printMatrix(arrA, size_for_each_node, matrix_size);

	// actual caluculation
	int* result = new int[matrix_size];
	for (int i = 0; i < size_for_each_node; i++)
	{
		result = multiplyStripe(arrA[i], arrB, matrix_size);
		for (int j = 0; j < matrix_size; j++)
		{
			arrC[i][j] = result[j];
		}
	}

	// print result
	// printMatrix(arrC, size_for_each_node, matrix_size);
	

	MPI_Gather(&(arrC[0][0]), matrix_size*size_for_each_node, MPI_INT, NULL, size_for_each_node*matrix_size, MPI_INT, 0, MPI_COMM_WORLD);



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

	//  get the matricies 
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



