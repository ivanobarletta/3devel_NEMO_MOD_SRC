MODULE solpetsc
   !!======================================================================
   !!                     ***  MODULE  solfet
   !! Ocean solver :  preconditionned conjugate gradient solver
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   sol_petsc    : petsc krylov subspace solver
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE sol_oce         ! ocean solver variables
   USE lib_mpp         ! distributed memory computing
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager
   USE lib_fortran     ! Fortran routines library
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sol_petsc_init, sol_petsc_finalize    ! 
   PUBLIC   petsc_set_mat , petsc_set_rhs, petsc_solve    ! 

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !! Petsc include files     
#  include "petsc/finclude/petscsys.h"
#  include "petsc/finclude/petscvec.h"
#  include "petsc/finclude/petscmat.h"
#  include "petsc/finclude/petscksp.h"
#  include "petsc/finclude/petscpc.h"
#  include "petsc/finclude/petscviewer.h"
#  include "petsc/finclude/petscvec.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: solpcg.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   PetscInt       ::   m, N
   PetscInt       ::   one, two, three, four, five
   PetscInt       ::   ji,jj, ishift, jshift 
   PetscInt       ::   LI, col5(5), col4(4), col3(3), row
   PetscErrorCode ::   ierr
   PetscBool      ::   flg, print_out
   PetscScalar    ::   dot,ione, values5(5), values4(4), values3(3)
   Vec            ::   x,rhs,residual, mv
   Mat            ::   Amat
   KSP            ::   ksp
   PC             ::   pc
   PetscScalar, pointer ::  xx_v(:) 
   PetscMPIInt    ::   rank, comm_size
   PetscDraw      ::   draw           
   PetscViewer    ::   vview, resview
   PetscViewer    ::   mview  ! DO NOT USE matview. It is the name of a PETSc routine 
   INTEGER, ALLOCATABLE, DIMENSION (:) :: south_coo, west_coo, east_coo, north_coo     
   INTEGER, ALLOCATABLE, DIMENSION (:) :: nfd_coo, nfd_coo2, nfd_coo_glo, nfd_coo_glo2
   INTEGER, ALLOCATABLE, DIMENSION (:) :: grid_points 
   INTEGER, ALLOCATABLE, DIMENSION (:) :: d_nnz, o_nnz  
   INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: COORD  
   INTEGER  :: nldi_pet, nlei_pet, nldj_pet, nlej_pet, prev_gp
   INTEGER  :: nx, ny , IDX        
   CHARACTER( LEN=32 ) :: filename     

CONTAINS
   
   SUBROUTINE sol_petsc_init
      CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)     

      five  = 5
      four  = 4
      three = 3
      two   = 2
      one   = 1

      ! global dimension of matrix
      N = ( jpiglo - 2 ) * ( jpjglo - 2 )

      ! define new subdomain edges  
      nldi_pet = nldi
      nlei_pet = nlei
      nldj_pet = nldj
      nlej_pet = nlej

      IF((MOD( narea, jpni ) == 1) .OR. (jpni == 1)) nldi_pet = nldi_pet + 1  ! west process
      IF((MOD( narea, jpni ) == 0) .OR. (jpni == 1)) nlei_pet = nlei_pet - 1  ! east process

      IF( (narea <=  jpni) .OR. (jpnj == 1) ) nldj_pet = nldj_pet + 1         ! first row of processor
      IF( (narea > (jpnj-1)*jpni ) .OR. (jpnj == 1)) nlej_pet = nlej_pet -1   ! last  row of processor

      ! # of local rows of matrix
      nx = nlei_pet - nldi_pet + 1
      ny = nlej_pet - nldj_pet + 1
      m = nx * ny

      ALLOCATE (  grid_points( mppsize ) )  
      ALLOCATE (  west_coo ( ny ),  east_coo ( ny ) )
      ALLOCATE ( south_coo ( nx ), north_coo ( nx ) )

      !! Only processes involved in North-Fold allocate 
      !! and manage these arrays             
      IF ( npolj >= 3 ) THEN
        ALLOCATE ( nfd_coo     ( nx )       , nfd_coo2     ( nx )       )
        ALLOCATE ( nfd_coo_glo ( jpiglo-2 ) , nfd_coo_glo2 ( jpiglo-2 ) ) 
        nfd_coo     = -99; nfd_coo2     = -99
        nfd_coo_glo = -99; nfd_coo_glo2 = -99
      END IF  

      ALLOCATE ( COORD ( jpi, jpj , 5 ) )  

      ALLOCATE ( d_nnz ( m ), o_nnz( m ) )  

      !! Initialize coordinates at boudaries  
      west_coo  = -99  ;  east_coo = -99
      south_coo = -99  ; north_coo = -99
      COORD = -99  

      ishift = nldi_pet - 1
      jshift = nldj_pet - 1

      !! Gather number of grid points of each process  
      grid_points = 0  
      grid_points(narea) = m  
      CALL mppsum_int_array(  jpnij, grid_points )
      prev_gp = 0
      IF ( narea > 1 ) prev_gp = SUM ( grid_points(1:narea-1) )

      !! Build coordinates at boundaries  
      DO jj = nldj_pet, nlej_pet
        DO ji = nldi_pet, nlei_pet
          IDX = ( jj-1-jshift) * nx + ji - ishift - 1     
          IF ( bmask(ji,jj ) /= 0 ) THEN
            IF ( jj == nlej_pet-1 .AND. npolj >=3 ) THEN 
              nfd_coo2(ji-ishift) = prev_gp + IDX
            END IF    
            IF ( jj == nldj_pet   ) south_coo(ji-ishift) = prev_gp + IDX
            IF ( jj == nlej_pet   ) north_coo(ji-ishift) = prev_gp + IDX
            IF ( ji == nldi_pet   ) west_coo(jj-jshift)  = prev_gp + IDX
            IF ( ji == nlei_pet   ) east_coo(jj-jshift)  = prev_gp + IDX
          END IF
        END DO
      END DO
        
      !! correction for edge processes in case of non-cyclic 
      !! boundary conditions  
      IF ( nbondi == -1 )  west_coo = -99
      IF ( nbondi ==  1 )  east_coo = -99
      IF ( nbondj == -1 ) south_coo = -99
      IF ( nbondj ==  1 ) north_coo = -99

      ! communicate the boundary conditions  
      CALL lbc_coo  

      !! Initialize arrays for preallocation  
      d_nnz = 1 ;  o_nnz = 0  

      !! Build COORD, d_nnz, o_nnz arrays  
      DO jj = nldj_pet, nlej_pet
        DO ji = nldi_pet, nlei_pet
          IDX = ( jj-1-jshift) * nx + ji - ishift
          row = prev_gp + IDX - 1
          COORD(ji,jj,3) = row        ! DIAG 
          IF ( bmask(ji,jj) /= 0 ) THEN
            !! SOUTH COEFFICIENT   
            IF ( gcp2(ji,jj,1) /= 0.d0 ) THEN 
              IF ( jj == nldj_pet ) THEN   
                COORD(ji,jj,1) = south_coo(ji-ishift); o_nnz( IDX ) = o_nnz( IDX ) + 1
              ELSE
                COORD(ji,jj,1) = row - nx            ; d_nnz( IDX ) = d_nnz( IDX ) + 1
              END IF              
            END IF 
            !! WEST COEFFICIENT (different treatment in case of cyclic condition)      
            IF ( gcp2(ji,jj,2) /= 0.d0 ) THEN 
              IF ( ji == nldi_pet ) THEN
                COORD(ji,jj,2) = west_coo(jj-jshift)
                IF ( nbondi == 2 ) THEN 
                  d_nnz( IDX ) = d_nnz( IDX ) + 1
                ELSE
                  o_nnz( IDX ) = o_nnz( IDX ) + 1
                END IF
              ELSE
                COORD(ji,jj,2) = row - 1             ; d_nnz( IDX ) = d_nnz( IDX ) + 1
              END IF
            END IF    
            !! EAST COEFFICIENT (different treatment in case of cyclic condition)  
            IF ( gcp2(ji,jj,3) /= 0.d0 ) THEN 
              IF ( ji == nlei_pet ) THEN
                COORD(ji,jj,4) = east_coo(jj-jshift)
                IF ( nbondi == 2 ) THEN
                  d_nnz( IDX ) = d_nnz( IDX ) + 1
                ELSE  
                  o_nnz( IDX ) = o_nnz( IDX ) + 1
                END IF
              ELSE
                COORD(ji,jj,4) = row + 1             ; d_nnz( IDX ) = d_nnz( IDX ) + 1
              END IF
            END IF
            !! NORTH COEFFICIENT       
            IF ( gcp2(ji,jj,4) /= 0.d0 ) THEN 
              IF ( jj == nlej_pet ) THEN   
                COORD(ji,jj,5) = north_coo(ji-ishift)
                IF ( npolj < 3 ) o_nnz( IDX ) = o_nnz( IDX ) + 1
              ELSE
                COORD(ji,jj,5) = row + nx            ; d_nnz( IDX ) = d_nnz( IDX ) + 1
              END IF              
            END IF                    
          END IF     
        END DO
      END DO

      9876 FORMAT(2i6,4i9)  
      9877 FORMAT(2i6,2a,2i9)  

      IF ( .FALSE. ) THEN  
        SELECT CASE ( narea )
          CASE ( 1  ) ; filename = "ij_coo_0.txt"
          CASE ( 2  ) ; filename = "ij_coo_1.txt"
          CASE ( 3  ) ; filename = "ij_coo_2.txt"
          CASE ( 4  ) ; filename = "ij_coo_3.txt"
          CASE ( 5  ) ; filename = "ij_coo_4.txt"
          CASE ( 6  ) ; filename = "ij_coo_5.txt"
          CASE ( 7  ) ; filename = "ij_coo_6.txt"
          CASE ( 8  ) ; filename = "ij_coo_7.txt"
          CASE ( 9  ) ; filename = "ij_coo_8.txt"
          CASE ( 10 ) ; filename = "ij_coo_9.txt"
          CASE ( 11 ) ; filename = "ij_coo_10.txt"
          CASE ( 12 ) ; filename = "ij_coo_11.txt"
          CASE ( 13 ) ; filename = "ij_coo_12.txt"
          CASE ( 14 ) ; filename = "ij_coo_13.txt"
          CASE ( 15 ) ; filename = "ij_coo_14.txt"
          CASE ( 16 ) ; filename = "ij_coo_15.txt"
          CASE ( 17 ) ; filename = "ij_coo_16.txt"
          CASE ( 18 ) ; filename = "ij_coo_17.txt"
          CASE ( 19 ) ; filename = "ij_coo_18.txt"
          CASE ( 20 ) ; filename = "ij_coo_19.txt"
          CASE ( 21 ) ; filename = "ij_coo_20.txt"
          CASE ( 22 ) ; filename = "ij_coo_21.txt"
          CASE ( 23 ) ; filename = "ij_coo_22.txt"
          CASE ( 24 ) ; filename = "ij_coo_23.txt"
        END SELECT  


        OPEN ( 811+narea, file=filename, action='write')  
        DO jj = nldj_pet, nlej_pet
          DO ji = nldi_pet, nlei_pet    
            write(811+narea,9876) ji,jj,d_nnz(row) ,o_nnz(row) ,COORD(ji,jj,3), COORD(ji,jj,1) 
            write(811+narea,9877) ji,jj,'         ','         ',COORD(ji,jj,3), COORD(ji,jj,2)
            write(811+narea,9877) ji,jj,'         ','         ',COORD(ji,jj,3), COORD(ji,jj,3)
            write(811+narea,9877) ji,jj,'         ','         ',COORD(ji,jj,3), COORD(ji,jj,4)
            write(811+narea,9877) ji,jj,'         ','         ',COORD(ji,jj,3), COORD(ji,jj,5)
          END DO
        END DO    

        CLOSE(811+narea)  

      END IF  

      CALL MatCreate( mpi_comm_opa,Amat,ierr)
      !CALL MatCreate( PETSC_COMM_WORLD,A,ierr)
      CALL MatSetSizes(Amat,m,m,N,N,ierr)
      CALL MatSetType(Amat, MATMPIAIJ,ierr)
      CALL MatSetFromOptions(Amat,ierr)

      !CALL MatMPIAIJSetPreallocation(Amat,PETSC_NULL_INTEGER,d_nnz,PETSC_NULL_INTEGER,o_nnz,ierr)
      CALL MatMPIAIJSetPreallocation(Amat,five,PETSC_NULL_INTEGER,two,PETSC_NULL_INTEGER,ierr)

      !! < RIGHT HAND SIDE >
      CALL VecCreate( mpi_comm_opa, rhs, ierr)
      !CALL VecCreate( PETSC_COMM_WORLD, rhs, ierr)
      CALL VecSetSizes(rhs, m, N, ierr)
      CALL VecSetFromOptions(rhs,ierr)

      !! Solution Vector  
      CALL VecDuplicate(rhs, x       , ierr)

      !! Krylov SubSpace         
      CALL KSPCreate( mpi_comm_opa, ksp, ierr)
      !CALL KSPCreate( PETSC_COMM_WORLD, ksp, ierr)
      CALL KSPSetOperators(ksp, Amat, Amat, ierr)
      CALL KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)
      !CALL KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED, ierr)
      CALL KSPSetFromOptions(ksp, ierr)

      !! --  
   END SUBROUTINE sol_petsc_init               

   SUBROUTINE sol_petsc_finalize
      !! Finalize Petsc Package  

      DEALLOCATE ( south_coo )  
      DEALLOCATE (  west_coo )  
      DEALLOCATE (  east_coo )  
      DEALLOCATE ( north_coo )  

      DEALLOCATE ( COORD )          

      DEALLOCATE ( d_nnz , o_nnz )  

      CALL MatDestroy(Amat, ierr)
      CALL VecDestroy(rhs, ierr)
      CALL VecDestroy(x, ierr)
      CALL KSPDestroy(ksp, ierr)
      CALL PetscFinalize(ierr)     

   END SUBROUTINE sol_petsc_finalize               

   !!=====================================================================

   SUBROUTINE petsc_set_mat  
      !! Set Matrix Values   
      IF( nn_timing == 1 )  CALL timing_start('petsc_set_mat')

      DO jj = nldj_pet, nlej_pet
        DO ji = nldi_pet, nlei_pet
          values5(1:5) = (/gcp2(ji,jj,1),gcp2(ji,jj,2),gcdmat2(ji,jj),gcp2(ji,jj,3),gcp2(ji,jj,4)/)
          col5(1:5)    = COORD(ji,jj,1:5)      
          CALL MatSetValues(Amat,one,COORD(ji,jj,3),five,col5,values5(1:5), INSERT_VALUES,ierr)
        END DO
      END DO

      CALL MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr)
      CALL MatAssemblyEnd  (Amat,MAT_FINAL_ASSEMBLY,ierr)

      IF( nn_timing == 1 )  CALL timing_stop('petsc_set_mat')

   END SUBROUTINE petsc_set_mat

   SUBROUTINE petsc_set_rhs
      !! Set RHS   
      IF( nn_timing == 1 )  CALL timing_start('petsc_set_rhs')

      DO jj = nldj_pet, nlej_pet
        DO ji = nldi_pet, nlei_pet
          CALL VecSetValue(rhs,COORD(ji,jj,3),gcb2(ji,jj),INSERT_VALUES,ierr)
        END DO
      END DO

      CALL VecAssemblyBegin(rhs,ierr)
      CALL VecAssemblyEnd(rhs,ierr)

      IF( nn_timing == 1 )  CALL timing_stop('petsc_set_rhs')

   END SUBROUTINE petsc_set_rhs          

   SUBROUTINE petsc_solve( kindic )
      INTEGER, INTENT(inout) ::   kindic   ! solver indicator, < 0 if the conver-
      !                                    ! gence is not reached: the model is stopped in step
      !                                    ! set to zero before the call of solver
  
      IF( nn_timing == 1 )  CALL timing_start('petsc_solve')

      ! < SET FIRST GUESS >
      CALL VecGetArrayF90(x,xx_v,ierr)
      DO jj = nldj_pet, nlej_pet
         DO ji = nldi_pet, nlei_pet
            xx_v(COORD(ji,jj,3) - prev_gp + 1) = gcx(ji,jj)
         END DO
      END DO
      CALL VecRestoreArrayF90(x,xx_v,ierr)

      ! < SOLVE SYSTEM >
      CALL KSPSolve(ksp, rhs, x, ierr)
      CALL KSPGetIterationNumber(ksp, niter, ierr)

      IF( niter == 0 ) kindic = -2

      ! < RETURN SOLUTION TO NEMO >
      CALL VecGetArrayReadF90(x,xx_v,ierr)
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            gcx(ji,jj) = xx_v(COORD(ji,jj,3) - prev_gp + 1) 
         END DO
      END DO
      CALL VecRestoreArrayReadF90(x,xx_v,ierr)

      CALL lbc_lnk( gcx, c_solver_pt, 1. )      ! Output in gcx with lateral b.c. applied
 
      IF( nn_timing == 1 )  CALL timing_stop('petsc_solve')

   END SUBROUTINE petsc_solve          

   SUBROUTINE mppsum_int_array( arr_size , ktab )
      !!----------------------------------------------------------------------
      !!                 ***  routine mppsum_int_array  ***
      !!
      !! ** Purpose :   Global integer sum or an array
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT (in) :: arr_size
      INTEGER, DIMENSION(arr_size), INTENT(inout) ::   ktab
      !!
      INTEGER :: ierror
      INTEGER, DIMENSION(arr_size) ::  iwork
      !!----------------------------------------------------------------------
      !
      CALL mpi_allreduce( ktab, iwork, arr_size, mpi_integer, mpi_sum, mpi_comm_opa, ierror )
      !
      ktab = iwork
      !
   END SUBROUTINE mppsum_int_array

   SUBROUTINE lbc_coo
      INTEGER :: ierr        
      INTEGER :: ijm1  
      INTEGER, ALLOCATABLE, DIMENSION (:) :: displ, nx_nfd    
      INTEGER :: north_comm_size  

      nx = nlei_pet - nldi_pet + 1  
      ny = nlej_pet - nldj_pet + 1 

      ishift = nldi_pet - 1  

      !! number of processes below 
      !! north belt    
      ijm1 = (jpnj-1)*jpni  

      !! Tags of messages:
      !! 1 - Send to West  / Receive from East    
      !! 2 - Send to East  / Receive from West    
      !! 3 - Send to South / Receive from North    
      !! 4 - Send to North / Receive from South   

      !! East-West Boundaries Migrations 
      SELECT CASE ( nbondi ) 
      CASE ( -1 )     
         CALL MPI_Send(  east_coo, ny, MPI_INT, noea,2, mpi_comm_opa, ierr)
         CALL MPI_Recv(  east_coo, ny, MPI_INT, noea,1, mpi_comm_opa, MPI_STATUS_IGNORE, ierr)         
      CASE ( 0  )  
         CALL MPI_Send(  east_coo, ny, MPI_INT, noea,2, mpi_comm_opa, ierr)
         CALL MPI_Send(  west_coo, ny, MPI_INT, nowe,1, mpi_comm_opa, ierr)
         CALL MPI_Recv(  west_coo, ny, MPI_INT, nowe,2, mpi_comm_opa, MPI_STATUS_IGNORE, ierr) 
         CALL MPI_Recv(  east_coo, ny, MPI_INT, noea,1, mpi_comm_opa, MPI_STATUS_IGNORE, ierr)       
      CASE ( 1 )  
         CALL MPI_Send(  west_coo, ny, MPI_INT, nowe,1, mpi_comm_opa, ierr) 
         CALL MPI_Recv(  west_coo, ny, MPI_INT, nowe,2, mpi_comm_opa, MPI_STATUS_IGNORE, ierr)
      END SELECT
      
      !! Only 1 Process in x direction  
      IF ( nbondi == 2 ) THEN   
        IF ( jperio == 1 .OR. jperio >= 3 ) THEN
          !! Swap east-west coordinates 
          east_coo = east_coo + west_coo            
          west_coo = east_coo - west_coo            
          east_coo = east_coo - west_coo            
        END IF     
      END IF 

      !! North-South Boundaries Migrations  
      SELECT CASE ( nbondj )
      CASE ( -1 )
         CALL MPI_Send( north_coo, nx, MPI_INT, nono,4, mpi_comm_opa, ierr)
         CALL MPI_Recv( north_coo, nx, MPI_INT, nono,3, mpi_comm_opa, MPI_STATUS_IGNORE, ierr)  
      CASE ( 0 ) 
         CALL MPI_Send( north_coo, nx, MPI_INT, nono,4, mpi_comm_opa, ierr)
         CALL MPI_Send( south_coo, nx, MPI_INT, noso,3, mpi_comm_opa, ierr)
         CALL MPI_Recv( south_coo, nx, MPI_INT, noso,4, mpi_comm_opa, MPI_STATUS_IGNORE, ierr)             
         CALL MPI_Recv( north_coo, nx, MPI_INT, nono,3, mpi_comm_opa, MPI_STATUS_IGNORE, ierr)             
      CASE ( 1 )
         CALL MPI_Send( south_coo, nx, MPI_INT, noso,3, mpi_comm_opa, ierr)        
         CALL MPI_Recv( south_coo, nx, MPI_INT, noso,4, mpi_comm_opa, MPI_STATUS_IGNORE, ierr)        
      END SELECT 

      north_comm_size = 1  

      !! Arrays for North-Fold Treatment 
      IF ( npolj >= 3 ) THEN
        IF ( jpni > 1 ) CALL MPI_Comm_size( ncomm_north, north_comm_size, ierr ) 
        ALLOCATE ( displ ( north_comm_size ) , nx_nfd( north_comm_size ) )
        !! grid points along x for north_comm processes 
        nx_nfd = grid_points(mppsize-north_comm_size+1:mppsize) 
        nx_nfd = nx_nfd / ny
        !! displacement array for ALL_GATHERV
        displ = 0
        IF ( ndim_rank_north >= 2 ) THEN
          DO ji=2,north_comm_size
            displ(ji) = sum(nx_nfd(1:ji-1)) 
          END DO
        END IF
        !WRITE(6,*) 'displacement    ', displ  
        IF ( jpni > 1 ) THEN     
          !! jj = jpj-2        
          CALL MPI_ALLGATHERV(nfd_coo2,nx,MPI_INT,nfd_coo_glo2,nx_nfd,displ,MPI_INT,ncomm_north, ierr)
        ELSE
          nfd_coo_glo2 = nfd_coo2      
        END IF

        DO ji=3,jpiglo-1
          nfd_coo_glo(ji-ishift) = nfd_coo_glo2(jpiglo-ji+ishift)       
        END DO
        IF ( jpni > 1 ) THEN
          !! scatter back values      
          DO ji=1,nx
            north_coo(ji) = nfd_coo_glo(ji+displ(narea-ijm1))  
          END DO
        ELSE
          north_coo = nfd_coo_glo
        END IF    
      END IF  

      IF ( npolj >= 3 ) THEN   
        DEALLOCATE( displ  )
        DEALLOCATE( nx_nfd )    
      END IF  

   END SUBROUTINE lbc_coo     
        
END MODULE solpetsc
