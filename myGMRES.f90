program linearSolv
  use kinddefs, only : dp
  use linalg
  use blas
  implicit none

  integer maxnodes,maxnnz,maxnsub,maxsearch,method,iwriteres,nsubiters,ntot,k,icycle,mypre
  parameter(maxnodes=15000,maxnnz=18*maxnodes,maxnsub=6)
  parameter(maxsearch=500)

  real(dp) :: tol,resnorm,droptol
  integer :: i,ierror,istat,fill_level,bwidth
  integer :: nnodes,nsub,nnz

  real(dp), dimension(:,:,:), allocatable ::  A!(maxnsub*maxnsub*maxnnz)   ! Matrix for Ax+b = 0
  real(dp), dimension(:,:), allocatable ::    q!(maxnsub*maxnodes)     ! Solution vector
  real(dp), dimension(:,:), allocatable ::    rhs!(maxnsub*maxnodes)    ! Right-hand-side
  real(dp), dimension(:,:), allocatable ::    res!(maxnsub*maxnodes)    ! Right-hand-side

  integer, dimension(:), allocatable :: ia!(maxnodes+1)
  integer, dimension(:), allocatable :: iau!(maxnodes)
  integer, dimension(:), allocatable :: ja!(maxnnz) ! Compressed row storage
  !integer, dimension(:), allocatable :: ka!(maxnnz) ! Compressed row storage
  !integer, dimension(:), allocatable :: kstart!(maxnodes) ! Arrays we need to try the least squares linear solver

  integer, dimension(:), allocatable :: ianew!(maxnodes+1)
  integer, dimension(:), allocatable :: iaunew!(maxnodes)
  integer, dimension(:), allocatable :: janew!(maxnnz) ! Compressed row storage
  integer, dimension(:), allocatable :: jnewtojold!(maxnnz) ! Compressed row storage
  integer, dimension(:), allocatable :: oldtonew
  integer, dimension(:), allocatable :: list

  integer, dimension(:), allocatable :: iafill!(maxnodes+1)
  integer, dimension(:), allocatable :: iaufill!(maxnodes)
  integer, dimension(:), allocatable :: jafill!(maxnnz) ! Compressed row storage
  real(dp), dimension(:,:,:), allocatable :: Afill!(maxnsub*maxnsub*maxnnz)   ! Matrix for Ax+b = 0
  integer:: nnzfill

  character(len=130)::matrix_file
  character(len=130)::param_file
  param_file="input.dat"

  open(UNIT=5,FILE=param_file,form='formatted',STATUS='old')

  read(5,*)
  read(5,*) matrix_file
  open(UNIT=299,FILE=matrix_file,form='formatted',STATUS='old')
  write(6,*) "from input matrix file ",matrix_file

  write(6,*) 'Solving equations'
  read(299,*)
  read(299,*)nnodes
  !write(6,*) "myGmres nnodes=",nnodes
  if(nnodes.gt.maxnodes)then
     write(6,*)'Increase maxnodes to at least',nnodes
     stop
  end if
  !read(299)nsub
  read(299,*)
  read(299,*)nsub
  !write(6,*) "myGmres nsub=",nsub
  !read(299)nnz
  read(299,*)
  read(299,*)nnz
  !write(6,*) "myGmres nnz=",nnz
  if(nnz.gt.maxnnz)then
     write(6,*)'Increase maxnnz to at least',nnz
     stop
  end if
  close(299)!matrix.bin
  
  !write(6,*)'1 nnodes=',nnodes
  !write(6,*)'1 nsub=',nsub
  !write(6,*)'1 nnz=',nnz

  allocate(ia(nnodes+1),STAT = istat)
  allocate(iau(nnodes),STAT = istat)
  allocate(ja(nnz),STAT = istat)

  allocate(ianew(nnodes+1),STAT = istat)
  allocate(iaunew(nnodes),STAT = istat)
  allocate(janew(nnz),STAT = istat)
  allocate(jnewtojold(nnz),STAT = istat)
  allocate(oldtonew(nnodes),STAT = istat)
  allocate(list(nnodes),STAT = istat)

  !allocate(ka(nnz),STAT = istat)
  !allocate(kstart(nnodes),STAT = istat)

  allocate(A(nsub,nsub,nnz),STAT = istat)
  allocate(q(nsub,nnodes),STAT = istat)
  allocate(rhs(nsub,nnodes),STAT = istat)
  allocate(res(nsub,nnodes),STAT = istat)

  ! Read in binary data file
  !call READA(nnodes,nnz,ia,ja,iau,A,rhs,nsub)!matrix.bin
  !also reorders using cuthill-mckee
  call READERf(matrix_file,nnodes,nsub,nnz,A,ia,iau,ja,ianew,iaunew,janew,list,oldtonew,jnewtojold,bwidth,rhs)!matrix.txt

  call READRESTART(nnodes,nsub,q)

  ! Choose point iterative or gmres
  !write(6,*)"Enter 1 for point iterative 2 for Kaczmarz (anything else will be gmres)"
  read(5,*)
  read(5,*) method
  write(6,*) "method :",method

  ! Solve for q

  !=====compute initial residual and display l2norm =======
  call compute_residual(nnodes,nsub,nnz,ia,ja,A,q,rhs,res,resnorm)
  write(6,*) "initial resnorm= ",resnorm
  !=============
  
  iwriteres = 1
  if(method.eq.1)then !Gauss Seidel solver
     !write(6,*)"Enter number of subiterations"
     read(5,*)
     read(5,*)nsubiters
     call SOLVE(nnodes,nsub,nnz,ia,ja,iau,A,q,rhs,nsubiters,iwriteres)
  else if(method.eq.2)then !Block Kaczmarz solve
     !write(6,*)"Enter number of subiterations"
     read(5,*)
     read(5,*)nsubiters
     call BKACZMARZ(nnodes,nsub,nnz,ia,ja,iau,A,q,rhs,nsubiters,iwriteres)
  else if(method.eq.3)then !Kaczmarz solver
     !call getKA(nnodes,nnz,ia,iau,ja,ka,kstart)
     !write(6,*)"Enter number of subiterations"
     read(5,*)
     read(5,*)nsubiters
     call KACZMARZ(nnodes,nsub,nnz,ia,ja,iau,A,q,rhs,nsubiters,iwriteres)
  else if(method.eq.4) then !GMRES Solve
     ntot = nsub*nnodes
     !write(6,*)"Enter nsearch icycle tol mypre nsubiters(only used for mypre 1 (Gauss Seidel) or mypre 2 (Kaczmarz))"
     read(5,*)
     read(5,*)nsubiters
     read(5,*)
     read(5,*)k
     read(5,*)
     read(5,*)icycle
     read(5,*)
     read(5,*)tol
     read(5,*)
     read(5,*)mypre
     read(5,*)
     read(5,*)fill_level
     write(6,*) "ILU fill level :",fill_level
     read(5,*)
     read(5,*)droptol

     if (mypre.eq.4) then!iluk preconditioner symbolic factorization
        ! if iluk, create a new set of ia,iau,ja,A with indexing for additional columns and copy original A into new Afill
        ! actual ILU factoring will occur inside gmres as an ilu(0) with the additional columns included
        allocate(iafill(nnodes+1),STAT = istat)
        allocate(iaufill(nnodes),STAT = istat)
        ! currently requires storing both A and Afill
 call ILUK_factorization(fill_level,nnodes,nsub,nnz,A,ianew,iaunew,janew,jnewtojold,bwidth,nnzfill,iafill,iaufill,jafill,Afill)
        ! GMRES is called to solve Afill only
 call GMRES(ntot,k,icycle,tol,ierror,mypre,nsubiters,nnzfill,nsub,nnodes,iafill,jafill,iaufill,Afill,rhs,q,fill_level,droptol) 
     else !if(mypre.eq.5)!function for ilut (threshold based fill in)
        write(6,*) "Calling gmres solver"
        call GMRES(ntot,k,icycle,tol,ierror,mypre,nsubiters,nnz,nsub,nnodes,ianew,janew,iaunew,A,rhs,q,fill_level,droptol)
     end if
  else if(method.eq.5) then! VARMS
     !call ARMS()
     ! to be implemented
  end if

  close(5)!input.dat

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
  if(ALLOCATED(A))deallocate(A,STAT=istat)
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

  stop
end program linearSolv

























































