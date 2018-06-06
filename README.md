# eigen-jacobi
Eigenvalues of a cauchy matrix using Jacobi rotation

Program to find the Eigenvalues of the Cauchy matrix using multiple Jacobi rotations.

Language : Fortran
Operating System : Ubuntu
Compiler : ifort
Install the ifort compiler for respective operating system.
Precision : Double Precision
To obtain double precision link the library –r16
Comments : In the .tar file, submitted, the number of iterations in the rotation counter loop are 10. I had kept only those many for testing. It needs to be turned into an infinite loop, and the limit condition will handle its exit. 

Language : C
Operating System : Ubuntu
Compiler : gcc
Libraries : Math library ( -lm ) to be included while compiling
Input : Command Line Input - The dimension of the square matrix( number of rows or columns)
Precision : Double precision
Concept :  The Cauchy matrix is made to undergo Jacobi rotations till the off-diagonal elements become insignificant.
Observe that the sum of squares of the elements remain constant.
The diagonal elements will tend to the eigenvalues of the matrix.

Author : Anant Shah
Email : anantshah200@gmail.com
