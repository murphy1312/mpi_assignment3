#include "stdafx.h" // visual studio only
#include <iostream>
#include <cstdio>
using namespace std;


void printMatrix() 
{
	int matrix[64];

	FILE *input;
	fopen_s(&input, "matrix.dat", "rb");
	if (!input) {
		return;
	}
	fread(matrix, sizeof(int), 64, input);
	fclose(input);
	for (int i = 0; i < 64; ++i) {
		if (i % 8 == 0)
			cout << endl;
		else
			cout << matrix[i] << " ";
	}
}


// function that generates a matrix in binary form on disk
void generateMatrixFile() {
	// generate a matrix of values that need to be written to disk in the form of a one dimensional array this will write out an 8x8 matrix
	int matrix[64] = {
		5, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1
	};

	// open up a file for writing
	FILE *output;
	fopen_s(&output, "matrix2.dat", "wb");

	// FILE *output = fopen("matrix.dat", "wb");
	if (!output) {
		return;
	}

	// do a simple fwrite to write the matrix to file
	fwrite(matrix, sizeof(int), 64, output);


	// close the file when we are finished writing
	fclose(output);

	printMatrix();

}

// the main function of the program
/*
int main(int argc, char** argv) 
{
	generateMatrixFile();
}
*/
