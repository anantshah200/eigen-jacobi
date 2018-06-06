/*Anant Shah
 * Date : 19-12-2017
 **/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

void mat_gen(double*  mat,int DIM ); /* Generate the Cauchy matrix */
double mat_sqr_norm( double* mat,int DIM ); /* Return the square of the elements of the matrix */
double mat_diag_sqr_norm( double* mat,int DIM ); /* Return the sum of squares of the diagonal elements */
double mat_offdiag_sqr_norm( double* mat,int DIM ); /* Returns the sum of squares of the off diagonal elements */
void jacobi( double* mat,double* mat_rot,int DIM ); /* Creates the rotation matrix */
void ROTATE( double* mat,double* mat_rot,double* mat_mul_temp,double* mat_mul,int DIM); /* Gives the matrix after one jacobi rotation */

int main(int argc,char **argv)
{
        /* Variable Initialization */
        double*         mat; /* Cauchy matrix : A(i,j) =1/i+j */
        int             DIM; /* Dimension of the square matrix */
        double*         mat_mul; /* Multiplication matrix : mat_rot(T)*mat*mat_rot */
        double*         mat_rot; /* Rotation matrix */
        double*         mat_mul_temp; /* Stores the temporary multiplication */
        int             rot_counter = 1; /* Jacobi rotation counter */
        double          sqr_norm; /* Square norm of all elements in the matrix */
        double          diag_sqr_norm; /* Square norm of the diagonal elements */
        double          offdiag_sqr_norm; /* Square norm of the off diagonal elements */
        int             i,j; /* Loop control variables */
        double          sqr_norm_solve; /* Square norm obtained as the sum of diagonal square norm and off diagonal square norm */
        double*         eigenvalues; /* Will store the diagonal elements of the rotated matrix i.e the eigenvalues */


        if(argc!=2)
        {
                printf("error : invalid number of arguments \n");
                exit(0);
        }

        /* Taking arguments from the command line */
        DIM =atoi( argv[1] );

       /* Allocate memory to the matrices */
        mat = (double *)malloc( DIM*DIM*sizeof(double) );
        mat_mul = (double *)malloc( DIM*DIM*sizeof(double) );
        mat_rot = (double *)malloc( DIM*DIM*sizeof(double) );
        mat_mul_temp = (double *)malloc( DIM*DIM*sizeof(double) );
        eigenvalues = (double *)malloc( DIM*sizeof(double) );

        /* Generate the Cauchy matric */
        mat_gen( mat,DIM );

        /* Initialize the temporary multiplication matrices */
        for(i=0;i<DIM;i++)
        {
                for(j=0;j<DIM;j++)
                {
                        mat_mul_temp[i*DIM+j] = 0.0;
                        mat_mul[i*DIM+j] = 0.0;
                }
        }

        sqr_norm = mat_sqr_norm( mat,DIM );
        printf("Square norm of all elements :  \t %.20f \n",sqr_norm);

        while( rot_counter!=0 )
        {
                jacobi( mat,mat_rot,DIM );
                ROTATE( mat,mat_rot,mat_mul_temp,mat_mul,DIM );

                /* Matrix after rotation obtained */
                diag_sqr_norm = mat_diag_sqr_norm( mat,DIM );
                offdiag_sqr_norm = mat_offdiag_sqr_norm( mat,DIM );

                sqr_norm_solve = diag_sqr_norm + offdiag_sqr_norm;
                /* Print the rotation matrix */
                for(i=0;i<DIM;i++)
                {
                        for(j=0;j<DIM;j++)
                        {
                                printf("%.20f ",mat_rot[i*DIM+j]);
                        }
                        printf("\n");
                }

                /* Print the matrix after each rotation */

                if( offdiag_sqr_norm<pow(10,-10) ) {


                        printf("Diagonal sqaure norm : %.20f \n",diag_sqr_norm);
                        printf("Off diagonal square norm : %.20f \n",offdiag_sqr_norm);

                        printf("The obtained square norm is : %.20f \n",sqr_norm_solve);

                        for(i=0;i<DIM;i++)
                        {
                                for(j=0;j<DIM;j++)
                                {
                                        printf("%.20f ",mat[i*DIM+j]);
                                }
                                printf("\n");
                        }
                        break;
                }
                else{
                        printf("Diagonal square norm : \t %.20f \n",diag_sqr_norm);
                        printf("Off Diagonal square norm : \t %.20f \n",offdiag_sqr_norm);

                        printf("The obtained square norm is : \t %.20f \n",sqr_norm_solve);
                        /* Print the matrix */
                        for(i=0;i<DIM;i++){
                                for(j=0;j<DIM;j++){
                                        printf("%.14f ",mat[i*DIM+j]);
                                }
                                printf("\n");
                        }
                }
                rot_counter++;
        }
       /* Assigning the diagonal elements to the <eigenvalues> array */
        for(i=0;i<DIM;i++)
                eigenvalues[i] = mat[i*DIM+i];

        printf("\n \n The eigenvalues of the Cauchy matrix are : \n");
        for(i=0;i<DIM;i++)
                printf("%d. %.14f \n",(i+1),eigenvalues[i]);
}

double mat_sqr_norm(double* mat,int DIM)
{
        /* Variable Initialization */
        int             i,j; /* Loop control variables */
        double          sqr_norm = 0.0; /* Square norm of all the elements */

        for(i=0;i<DIM;i++)
                for(j=0;j<DIM;j++)
                        sqr_norm+=mat[i*DIM+j]*mat[i*DIM+j];

        return sqr_norm;
}

double mat_diag_sqr_norm(double* mat,int DIM)
{
        /* Variable Initialization */
        int             i,j; /* Loop control variables */
        double          diag_sqr_norm = 0.0; /* Square norm of all diagonal elements */

        for(i=0;i<DIM;i++)
                        diag_sqr_norm+=mat[i*DIM+i]*mat[i*DIM+i];

        return diag_sqr_norm;
}

double mat_offdiag_sqr_norm(double* mat,int DIM)
{
        /* Variable Initialization */
        int             i,j; /* Loop control variables */
        double          offdiag_sqr_norm = 0.0; /* Square norm of off diagonal elements */
        for(i=0;i<DIM;i++)
                for(j=0;j<DIM;j++)
                {
                        if(i!=j)
                        {
                                offdiag_sqr_norm+=mat[i*DIM+j]*mat[i*DIM+j];
                        }
                }

        return offdiag_sqr_norm;
}

void ROTATE(double* mat,double* mat_rot,double* mat_mul_temp,double* mat_mul,int DIM)
{
        /* Variable Initialization */
        int             i,j,k; /* Loop initialization variables */

        for(i=0;i<DIM;i++)
                for(j=0;j<DIM;j++)
                        for(k=0;k<DIM;k++)
                                mat_mul_temp[i*DIM+j] += mat_rot[k*DIM+i]*mat[k*DIM+j];

        for(i=0;i<DIM;i++)
                for(j=0;j<DIM;j++)
                        for(k=0;k<DIM;k++)
                                mat_mul[i*DIM+j] += mat_mul_temp[i*DIM+k]*mat_rot[k*DIM+j];

        /* Copy the result back into mul */
        for(i=0;i<DIM;i++)
                for(j=0;j<DIM;j++)
                        mat[i*DIM+j] = mat_mul[i*DIM+j];

        /* Initializa the temporary matrices back to 0 */
        for(i=0;i<DIM;i++)
                for(j=0;j<DIM;j++)
                {
                        mat_mul[i*DIM+j] = 0.0;
                        mat_mul_temp[i*DIM+j] = 0.0;
                }
}

void jacobi(double* mat,double* mat_rot,int DIM)
{
        double          s,c; /* Sin and Cosine to be added to the rotation matrix */
        double          max; /* Maximum absolute values in the cauchy matrix */
        double          tau; /* Value of tau */
        int             i,j; /* Loop control variables */
        int             row_num,col_num; /* Row and column of element with maximum absolute values */
        double          t_one,t_two; /* Two roots of the quadratic; Choose which one is smaller */

        max = fabs( mat[1] );
        row_num = 0;
        col_num = 1;
        /* Find the maximum absolute value in the matrix */
        for(i=0;i<DIM;i++)
                for(j=0;j<DIM;j++)
                {
                        if( ( fabs(mat[i*DIM+j]) > max ) && (i!=j) )
                        {
                                max = fabs(mat[i*DIM + j]);
                                row_num = i;
                                col_num = j;
                        }
                }

        tau = (mat[col_num*DIM+col_num] - mat[row_num*DIM+row_num])/(2.0*mat[row_num*DIM+col_num]);

        t_one = (-tau) + sqrt( 1+tau*tau );
        t_two = (-tau) - sqrt( 1+tau*tau );

        if( fabs( t_one ) < fabs( t_two ) )
        {
                c = 1.0/sqrt( 1+t_one*t_one );
                s = c*t_one;
        }
        else
        {
                c = 1.0/sqrt( 1+t_two*t_two );
                s =  c*t_two;
                printf("%.14f c \n",c);
        }


        /* Rotation mtrix is such that mat(row_num,row_num) = c , mat(col_num,col_num = c) , mat(row_num,col_num) = s, mat(col_num,row_ 
           num) = (-s) , all other diagonal elements = 1, all other off diagonal elements = 0 */
        /* Assign these values into the rotation matrix */
        for(i=0;i<DIM;i++)
        {
                for(j=0;j<DIM;j++)
                {
                        if( (i==row_num) && (j==row_num) )
                                mat_rot[i*DIM+j] = c;
                        else if( (i==col_num) && (j==col_num) )
                                mat_rot[i*DIM+j] = c;
                        else if( (i==row_num) && (j==col_num) )
                                mat_rot[i*DIM+j] = s;
                        else if( (i==col_num) && (j==row_num) )
                                mat_rot[i*DIM+j] = (-s);

                        else if( i==j )
                                mat_rot[i*DIM+j] = 1.0;
                        else
                                mat_rot[i*DIM+j] = 0.0;
                }
        }
}

void mat_gen(double* mat,int DIM)
{
        /* Variable Initialization */
        int             i,j; /*Loop control variables */

        /* Fill in the values of the Cauchy matrix */
        for(i=0;i<DIM;i++)
                for(j=0;j<DIM;j++)
                        mat[i*DIM+j] = (double) 1.0/( (i+1) + (j+1) );
}

