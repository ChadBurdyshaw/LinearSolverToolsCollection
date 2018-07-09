program linearSolv
  use kinddefs, only : dp,solver_params_type
  !use mpi
  use pcomm
  use linalg
  use blas
  use armslib
  use arms2
  !use hb_file_module
  implicit none

  integer :: maxnodes,maxnnz,maxnsub,maxsearch,iwriteres
  integer :: nsubiters,nnodes,nsub,nnz
  !integer :: k,icycle,mypre,method
  parameter(maxnodes=15000,maxnnz=18*maxnodes,maxnsub=6)
  parameter(maxsearch=500)

  real(dp) :: convgtol,resnorm,droptol
  integer :: i,j,ierror,istat,bwidth
  !integer::mreord,rscal,cscal,rscal_S,cscal_S,pq_S
  !integer :: fill_L,fill_U,fill_EU,fill_LF,fill_S,fill_LS,fill_US
  !integer :: mformat

  real(dp), pointer, dimension(:,:,:) ::  A!(maxnsub*maxnsub*maxnnz)   ! Matrix for Ax+b = 0
  real(dp), pointer, dimension(:,:) ::    q!(maxnsub*maxnodes)     ! Solution vector
  real(dp), pointer, dimension(:,:) ::    rhs!(maxnsub*maxnodes)    ! Right-hand-side
  real(dp), pointer, dimension(:,:) ::    res!(maxnsub*maxnodes)    ! Right-hand-side

  real(dp), pointer, dimension(:,:)::x
  real(dp), pointer, dimension(:,:)::b

  integer, pointer, dimension(:) :: ia!(maxnodes+1)
  integer, pointer, dimension(:) :: iau!(maxnodes)
  integer, pointer, dimension(:) :: ja!(maxnnz) ! Compressed row storage

  integer, pointer, dimension(:) :: ianew!(maxnodes+1)
  integer, pointer, dimension(:) :: iaunew!(maxnodes)
  integer, pointer, dimension(:) :: janew!(maxnnz) ! Compressed row storage
  integer, pointer, dimension(:) :: jnewtojold!(maxnnz) ! Compressed row storage
  integer, pointer, dimension(:) :: oldtonew
  integer, pointer, dimension(:) :: list

  integer, pointer, dimension(:) :: iafill!(maxnodes+1)
  integer, pointer, dimension(:) :: iaufill!(maxnodes)
  integer, pointer, dimension(:) :: jafill!(maxnnz) ! Compressed row storage
  real(dp), pointer, dimension(:,:,:) :: Afill!(maxnsub*maxnsub*maxnnz)   ! Matrix for Ax+b = 0

  integer:: nnzfill

  integer :: my_nnz
  integer,   pointer, dimension(:) :: my_ia!(maxnodes+1)
  integer,   pointer, dimension(:) :: my_iau!(maxnodes)
  integer,   pointer, dimension(:) :: my_ja!(maxnnz) ! Compressed row storage
  real(dp),  pointer, dimension(:,:,:) :: my_A!(maxnsub*maxnsub*maxnnz)   ! Matrix for Ax+b = 0

  real(dp)::sum

  !arms data types
  type(parms_Operator),pointer ::op
  type(parms_Mat) :: self
  type(parms_FactParam),target :: param
  type(cs_type),pointer :: mat
  type(p4_type),pointer::levmat
  integer :: bsize !min block size for independent sets
  integer :: nlev

  type(solver_params_type)::sp

  integer :: iretval,ierr,rank,np

  !character(len=130)::matrix_file
  nullify(op)
  nullify(mat)


       !call MPI_INIT(iretval)
       !call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
       !call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

write(6,*)"MPI results np=",np," rank=",rank


  !call READ_PARAMFILE(matrix_file,mformat,method,mypre,nsubiters,k,icycle,tol,droptol, &
  !                    fill_L,fill_U,fill_EU,fill_LF,fill_S,fill_LS,fill_US, &
  !                    bsize,nlev,mreord,rscal,cscal,rscal_S,cscal_S,pq_S)

  call readparams(sp)

  write(6,*) "from input matrix file ",sp%matrix_file
  write(6,*) 'Solving equations'


  if(sp%mformat.eq.1) then
     call READERf(sp%matrix_file,nnodes,nsub,nnz,A,ia,iau,ja,ianew,iaunew,janew,list,oldtonew,jnewtojold,bwidth,rhs)!matrix.txt
  else if(sp%mformat.eq.2) then
     call READER_HB(sp%matrix_file,nnodes,nsub,nnz,A,ia,iau,ja,ianew,iaunew,janew,list,oldtonew,jnewtojold,bwidth,rhs)!matrix.txt
  end if



  call READRESTART(nnodes,nsub,q)

  ! Solve for q
  allocate(res(nsub,nnodes),STAT = istat)
  !=====compute initial residual and display l2norm =======
  call compute_residual(nnodes,nsub,nnz,ia,ja,A,q,rhs,res,resnorm)
  write(6,*) "initial resnorm= ",resnorm
  !=============
  
  iwriteres = 1
  if(sp%method.eq.1)then !Gauss Seidel solver
write(6,*)"method is Gauss Seidel"
     call SOLVE(nnodes,nsub,nnz,ia,ja,iau,A,q,rhs,sp%nsubiters,iwriteres)
  else if(sp%method.eq.2)then !Block Kaczmarz solve
write(6,*)"method is Block Kaczmarz"
     call BKACZMARZ(nnodes,nsub,nnz,ia,ja,iau,A,q,rhs,sp%nsubiters,iwriteres)
  else if(sp%method.eq.3)then !Kaczmarz solver
write(6,*)"method is Kaczmarz"
     call KACZMARZ(nnodes,nsub,nnz,ia,ja,iau,A,q,rhs,sp%nsubiters,iwriteres)

  else if(sp%method.eq.4) then !GMRES Solve ======================================

     write(6,*)"method is GMRES"
     my_nnz=nnz
     my_ia=>ia
     my_ja=>ja
     my_iau=>iau
     my_A=>A     

     select case(sp%precon)
     case(5)!iluk preconditioner symbolic factorization
        write(6,*)"preconditioner is ILUK(",sp%fill_BL,")"
        ! if iluk, create a new set of ia,iau,ja,A with indexing for additional columns and copy original A into new Afill
        ! actual ILU factoring will occur inside gmres as an ilu(0) with the additional columns included
        allocate(iafill(nnodes+1),STAT = istat)
        allocate(iaufill(nnodes),STAT = istat)
        ! currently requires storing both A and Afill
        call ILUK_factorization(sp%fill_BL,nnodes,nsub,nnz,A,ianew,iaunew,janew, &
             jnewtojold,bwidth,nnzfill,iafill,iaufill,jafill,Afill)
        
        my_nnz=nnzfill
        my_ia=>iafill
        my_ja=>jafill
        my_iau=>iaufill
        my_A=>Afill
        
     case(6)!function for ilut (threshold based fill in)
           write(6,*)"preconditioner is ILUT"
     case(7)!ARMS
        write(6,*)"preconditioner is VARMS"
        !ceb much of this will have to be reworked to be efficient in memory and time
        call create_parms_Mat(self)!calls create_parms_Map
        call setupCS(mat,nnodes,nsub,nsub,1)!creates mat object
        call copy_to_CS(A,ia,iau,ja,mat)!allocates rows and copies A to mat object
        call create_FactParam(param, &
             sp%fill_BL,sp%fill_BU,sp%fill_EU,sp%fill_LF,sp%fill_S,sp%fill_SL,sp%fill_SU, &
             sp%droptol,nnodes,nsub,nsub,sp%nlev,sp%block_size,sp%krylovdirs, &
             sp%mreord,sp%rscale_B,sp%cscale_B,sp%rscale_S,sp%cscale_S,sp%pq_S)
        call parms_OperatorCreate(op)
        call arms_factorization(self,param,mat,op)
        !x=>q
        !b=>rhs
     end select




     write(6,*) "Calling gmres solver"
     call GMRES(nnodes,sp%krylovdirs,sp%restarts,sp%convgtol,ierror,sp%precon,sp%nsubiters,my_nnz,nsub,&
                !my_ia,my_ja,my_iau,my_A,rhs,q,sp%fill_BL,sp%droptol,& 
                my_ia,my_ja,my_iau,my_A,res,q,sp%fill_BL,sp%droptol,& 
                op) 



  else if(sp%method.eq.5) then! VARMS =====================================================

     write(6,*)"method is VARMS"
     
     call create_parms_Mat(self)!calls create_parms_Map
     call setupCS(mat,nnodes,nsub,nsub,1)!creates mat object
     call copy_to_CS(A,ia,iau,ja,mat)!allocates rows and copies A to mat object
     call create_FactParam(param, &
          sp%fill_BL,sp%fill_BU,sp%fill_EU,sp%fill_LF,sp%fill_S,sp%fill_SL,sp%fill_SU, &
          sp%droptol,nnodes,nsub,nsub,sp%nlev,sp%block_size,sp%krylovdirs, &
          sp%mreord,sp%rscale_B,sp%cscale_B,sp%rscale_S,sp%cscale_S,sp%pq_S)
     call parms_OperatorCreate(op)
     call arms_factorization(self,param,mat,op)
     write(6,*)"post arms_factorization"
     write(6,*)""
     write(6,*)"pre arms_sol"
     x=>q
     b=>rhs
     call parms_arms_sol_vcsr(op, b, x)
     write(6,*)"post arms_sol"
     
     !do i=1,mat%n
     !   do j=1,mat%nsubrow
     !      write(6,*)"x(",j,i,")=",x(j,i)
     !   end do
     !end do
     
     
     !ceb temp convergence check for 1's vector solution
     sum=0.
     do i=1,nnodes
        do j=1,nsub
           sum=sum+abs(i-x(j,i))
        end do
     end do
     sum=sum/(nnodes*nsub)
     write(6,*) "ARMS sol resid=",sum
     
     
  end if
  
  !close(5)!input.dat

  !=====compute final residual and display l2norm =======
  call compute_residual(nnodes,nsub,nnz,ia,ja,A,q,rhs,res,resnorm)
  write(6,*) "final resnorm= ",resnorm
  !=============
  !
  !write(6,*)"checkpt end mygmres1"
  call WRITEOUT(nnodes,nsub,q)
  !write(6,*)"checkpt end mygmres2"
  deallocate(ianew,STAT=istat)
  deallocate(iaunew,STAT=istat)
  deallocate(janew,STAT=istat)
  !write(6,*)"checkpt end mygmres3"
  deallocate(jnewtojold,STAT=istat)
  deallocate(oldtonew,STAT=istat)
  !deallocate(list,STAT=istat)
  !write(6,*)"checkpt end mygmres4"
  deallocate(iafill,STAT=istat)
  deallocate(iaufill,STAT=istat)
  deallocate(jafill,STAT=istat)
  deallocate(Afill,STAT=istat)
  !write(6,*)"checkpt end mygmres5"
  !if(ALLOCATED(A))deallocate(A,STAT=istat)
  if(associated(A))deallocate(A,STAT=istat)
  !write(6,*)"checkpt end mygmres6"
  deallocate(ia,STAT=istat)
  deallocate(iau,STAT=istat)
  deallocate(ja,STAT=istat)
  !write(6,*)"checkpt end mygmres6.5"
  !deallocate(ka,STAT=istat)
  !deallocate(kstart,STAT=istat)
  !write(6,*)"checkpt end mygmres6.6"
  !  deallocate(rhs,STAT=istat)
  deallocate(res,STAT=istat)
  !deallocate(q,STAT=istat) !this kills the code for some reason
  ! my be causing illegal memory access somewhere 
  write(6,*)"checkpt end mygmres final"

   !call MPI_COMM_FREE(MPI_COMM_WORLD)
   !call MPI_FINALIZE(ierr)

  stop
end program linearSolv

























































