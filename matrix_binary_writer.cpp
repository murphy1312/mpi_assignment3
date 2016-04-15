// simple program to take a matrix and write it to disk in the form of a binary file

// includes that are necessary to make the program work
#include <iostream>
#include <cstdio>

#include <stdio.h>
using namespace std;


void printMatrix() {
	double matrix[64];
	FILE *input = fopen("matrix.dat", "rd");
	if (!input) {
		return;
	}
	fread(matrix, sizeof(double), 64, input);
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

	printMatrix();

}

// the main function of the program
int main(int argc, char** argv) {
	generateMatrixFile();
}
