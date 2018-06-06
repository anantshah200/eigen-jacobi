! Anant Shah
! Jacobi Rotation
! Date : 27-12-2017

	program JacobiRotation

        implicit none

        integer :: dim_mat ! Matrix dimension
        real(kind = 16),dimension(:,:),allocatable               :: mat ! Cauchy matrix 
        real(kind = 16),dimension(:,:),allocatable               :: mat_rot ! Roation matrix
        real(kind = 16),dimension(:),allocatable                 ::eigenvalues
        real :: square_norm
        real :: diag_square_norm
        real :: offdiag_square_norm
        integer :: i,j ! Loop control variables
        real,parameter  :: limit = 1.0e-10
        integer :: rot_count

        print *,"Enter the dimension of the matrix : "
        read *,dim_mat ! Read the dimension of the matrix; users input

        ! Allocate memory to the rotation matrix and the Cauchy matrix
        allocate(mat(dim_mat,dim_mat),mat_rot(dim_mat,dim_mat) )
        allocate( eigenvalues(dim_mat) )

        mat = mat_gen(dim_mat)

        square_norm = mat_sqr_norm(mat,dim_mat)
        print *,"Square =",square_norm
        do rot_count = 1,10
                mat_rot = jacobi(mat,dim_mat)
                mat = rotate(mat,mat_rot,dim_mat)

                ! Printing the rotation matrix
                print *,mat_rot

                diag_square_norm = mat_diag_sqr_norm(mat,dim_mat)
                offdiag_square_norm = mat_offdiag_sqr_norm(mat,dim_mat)
                if( offdiag_square_norm .LT. limit ) then
                        print *,"Diagoanl square =",diag_square_norm
                        print *,"Off Diagonal Square =",offdiag_square_norm
                        print *,mat
                        print *,'Exiting the program'
                        exit
                else
                        print *,"Diagonal Square =",diag_square_norm
                        print *,"Off diagonal Square =",offdiag_square_norm
                        print *,mat
                end if
        end do

        do i=1,dim_mat
                eigenvalues(i) = mat(i,i)
        end do
        print *,"Eigenvalues : "
        print *,eigenvalues

        contains
        function        mat_gen(dim_mat) result(mat)

                integer                    :: i,j ! Loop control variables
                real(kind = 16),dimension(:,:),allocatable      :: mat
                integer ,intent(in)        :: dim_mat ! Matrix dimension                

                ! Allocate the memory
                allocate( mat(dim_mat,dim_mat) )
                outer : do i=1,dim_mat
                        inner : do j=1,dim_mat
                                mat(i,j) = 1.0/( i+j )
                        end do inner
                end do outer
        end function    mat_gen

        real*8 function        mat_sqr_norm(mat,dim_mat) result(elements_sum)
                integer                 :: i,j ! Loop control variables
                real(kind = 16),dimension(:,:),allocatable    :: mat ! Cauchy matrix input
                integer    :: dim_mat ! Dimension
                elements_sum = 0.0

                outer1   : do i=1,dim_mat
                        inner1 : do j=1,dim_mat
                                elements_sum = elements_sum + mat(i,j)**2
                        end do inner1
                end do outer1
        end function    mat_sqr_norm

        real*8 function        mat_diag_sqr_norm(mat,dim_mat) result(dia_sum)
                integer         :: i,j ! Loop control variables
                real(kind = 16),dimension(:,:),allocatable   :: mat ! Cauchy matrix input
                integer,intent(in)      :: dim_mat

                dia_sum = 0.0
                do i=1,dim_mat
                        dia_sum = dia_sum + mat(i,i)**2
                end do
        end function mat_diag_sqr_norm

        real*8 function        mat_offdiag_sqr_norm(mat,dim_mat) result(off_sum)
                integer                 :: i,j ! Loop control variables
                real(kind = 16),dimension(:,:),allocatable   :: mat ! Cauchy matrix
                integer,intent(in)      :: dim_mat ! Dimension of the matrix

                off_sum = 0.0

                outer2 : do i=1,dim_mat
                        inner2 : do j=1,dim_mat
                                if ( i .NE. j ) then
                                        off_sum = &
                                        off_sum + mat(i,j)**2
                                end if
                        end do inner2
                end do outer2
        end function mat_offdiag_sqr_norm

        function        jacobi(mat,dim_mat) result(mat_rot)
                integer,intent(in)                      :: dim_mat ! Dimension of matrix
                real(kind = 16),dimension(:,:),allocatable       :: mat ! Cauchy matrix
                real(kind = 16),dimension(:,:),allocatable       :: mat_rot ! Rotation

                real :: mat_max
                real :: s,c,tau
                integer :: row_num,col_num
                real :: t_one,t_two

                allocate( mat_rot(dim_mat,dim_mat) )


                row_num = 1
                col_num = 2
                mat_max = abs( mat(row_num,col_num) )

                outer3 : do i=1,dim_mat
                        inner3 : do j=1,dim_mat
                                if( i .NE. j) then
                                        if( abs( mat(i,j) ) .GT. mat_max ) then
                                                mat_max = abs( mat(i,j) )
                                                row_num = i
                                                col_num = j
                                        end if
                                end if
                        end do inner3
                end do outer3

                tau = ( mat(col_num,col_num) - mat(row_num,row_num) )&
                      /( 2.0 * mat(row_num,col_num) )
                t_one = (-tau) + sqrt(1 + tau**2)
                t_two = (-tau) - sqrt(1 + tau**2)

                if( abs(t_one) .LT. abs(t_two) ) then
                        c = 1.0/sqrt(1+t_one**2)
                        s = c*t_one
                else
                        c = 1.0/sqrt(1+t_two**2)
                        s = c*t_two
                end if

                print *,row_num
               print *,col_num

                outer4 : do i=1,dim_mat
                        inner4 : do j=1,dim_mat
                        if( (i .EQ. row_num) .AND. (j .EQ. row_num) )  then
                                mat_rot(i,j) = c
                        else if( (i .EQ. col_num) .AND. (j .EQ. col_num) ) then
                                mat_rot(i,j) = c
                        else if( (i .EQ. row_num) .AND. (j .EQ. col_num) ) then
                                mat_rot(i,j) = s
                        else if( (i .EQ. col_num) .AND. (j .EQ. row_num) ) then
                                mat_rot(i,j) = (-s)
                        else if( i .EQ. j) then
                                mat_rot(i,j) = 1.0
                        else
                                mat_rot(i,j) = 0.0
                        end if
                        end do inner4
                end do outer4
        end function    jacobi
        function rotate(mat,mat_rot,dim_mat) result(mat_temp)
                integer,intent(in)         :: dim_mat ! Dimension of the matrix
                real(kind = 16),dimension(:,:),allocatable :: mat
                real(kind = 16),dimension(:,:),allocatable :: mat_rot
                real(kind = 16),dimension(:,:),allocatable :: mat_temp

                allocate( mat_temp(dim_mat,dim_mat) )

                mat_temp = matmul(matmul(transpose(mat_rot),mat),mat_rot)


        end function rotate
    	end
