module linalg
  use kinddefs, only : dp
  use blas
  !use armslib
  use ilut
  use pilu_lib
  !use hb_file_module
  !use pcomm
implicit none
contains

!=================================== READA ===================================
!
! Read in A matrix
!
!=============================================================================
  subroutine READA(nnodes,nnz,ia,ja,iau,A,rhs,nbsize)
    integer,intent(inout):: nnodes,nnz,nbsize
    real(dp),intent(inout) :: A(nbsize,nbsize,nnz),rhs(nbsize,nnodes)
    integer, intent(inout) ::ia(nnodes+1),iau(nnodes),ja(nnz)
    integer i,j,k,m,jstart,jend    

    character(len=130)::matrix_file
    matrix_file="matrix.bin"
    open(UNIT=299,FILE=matrix_file,form='unformatted',STATUS='unknown')
   
    !read(299) !Dimension
    read(299)nnodes
    !read(299) !block size
    read(299)nbsize
    !read(299)!number of non zero entries in matrix
    read(299)nnz
  
    !write(6,*)'2 nnodes=',nnodes
    !write(6,*)'2 nsub=',nbsize 
    !write(6,*)'2 nnz=',nnz

    !read(299)!"IA dimension=",nnodes+1
    do  i = 1,nnodes+1
       read(299)ia(i)
    end do 

    !read(299)!"IAU dimension=",nnodes
    do  i = 1,nnodes
       read(299)iau(i)
    end do  

    !read(299)!"JA nnz=",nnz
    do i = 1,nnz
       read(299)ja(i)
    end do
   
    !read(299)!"A dimension=",nsub*nnz," (nsub*nnz)"
    do  i = 1,nnodes
       jstart = ia(i)
       jend   = ia(i+1) - 1
       do  j = jstart,jend
          do m = 1,nbsize
             read(299)(A(m,k,j),k=1,nbsize)
          end do
       end do
    end do

    !read(299)!"Right hand side dimension=",nnodes
    do  i = 1,nnodes
       read(299)(rhs(m,i),m=1,nbsize)
       ! Reverse the sign because we solve Ax+b=0 but this was generated
       ! (in converter) as Ax=b
       do m = 1,nbsize
          rhs(m,i) = -rhs(m,i)
       end do
    end do
    
    close(299)
    return
  end subroutine READA
!===================================================================





!===================================================================
subroutine READ_PARAMFILE(matrix_file,format,method,mypre,nsubiters,k,icycle,tol,droptol, &
                          fill_L,fill_U,fill_EU,fill_LF,fill_S,fill_LS,fill_US, &
                          bsize,nlev,mreord,rscal,cscal,rscal_S,cscal_S,pq_S)

  character(len=130),intent(inout) :: matrix_file
  integer,intent(inout) :: format,method,mypre,nsubiters,k,icycle
  integer,intent(inout) :: fill_L,fill_U,fill_EU,fill_LF,fill_S,fill_LS,fill_US
  integer,intent(inout) :: bsize,nlev,mreord,rscal,cscal,rscal_S,cscal_S,pq_S
  real(dp),intent(inout):: tol,droptol

  character(len=130)::param_file
  param_file="input.dat"
  open(UNIT=5,FILE=param_file,form='formatted',STATUS='old')

  read(5,*)
  read(5,*) matrix_file
  read(5,*)
  read(5,*) format
  read(5,*)
  read(5,*) method
  write(6,*) "method :",method
  read(5,*)
  read(5,*)mypre
  write(6,*) "preconditioner :",mypre
  read(5,*)
  read(5,*)nsubiters
  read(5,*)
  read(5,*)k
  write(6,*) "krylov search dirs :",k
  read(5,*)
  read(5,*)icycle
  read(5,*)
  read(5,*)tol
  write(6,*) "convg tol :",tol
  read(5,*)
  read(5,*)droptol
  write(6,*) "ilut drop tol :",droptol
  read(5,*)
  read(5,*)fill_L
  write(6,*) "B ILU fill level L:",fill_L
  read(5,*)
  read(5,*)fill_U
  write(6,*) "B ILU fill level U:",fill_U
  read(5,*)
  read(5,*)fill_EU
  write(6,*) "EU fill level :",fill_EU
  read(5,*)
  read(5,*)fill_LF
  write(6,*) "LF fill level :",fill_LF
  read(5,*)
  read(5,*)fill_S
  write(6,*) "S fill level :",fill_S
  read(5,*)
  read(5,*)fill_LS
  write(6,*) "LS fill level :",fill_LS
  read(5,*)
  read(5,*)fill_US
  write(6,*) "US fill level :",fill_US
  read(5,*)
  read(5,*)bsize
  write(6,*)" min block size:",bsize
  read(5,*)
  read(5,*)nlev
  write(6,*)"Number of ARMS levels:",nlev
  read(5,*)
  read(5,*)mreord
  write(6,*)"reorder type:",mreord
  read(5,*)
  read(5,*)rscal !row scaling
  write(6,*)"B scale rows:",rscal
  read(5,*)
  read(5,*)cscal !column scaling
  write(6,*)"B scale colss:",cscal
  read(5,*)
  read(5,*)rscal_S !row scaling
  write(6,*)"S scale rows:",rscal_S
  read(5,*)
  read(5,*)cscal_S !column scaling
  write(6,*)"S scale cols:",cscal_S
  read(5,*)
  read(5,*)pq_S !pq permutation on S
  write(6,*)"PQ permute S:",pq_S
  
  close(5) 
  return
end subroutine READ_PARAMFILE
!===================================================================



!============================ READERf ================================
!
! Reads the matrix in Chad's format and converts to mine
! The only difference is that he reads in nnodes+1 iau's
! This routine also gives the option to transpose 
! or divide by the diagonal
!
!===================================================================
!subroutine READERf(nnodes,nsub,nnz,A,ia,iau,ja,list,tag,ianew,iaunew,janew,oldtonew,rhs)
subroutine READERf(matrix_file,nnodes,nsub,nnz,A,ia,iau,ja,ianew,iaunew,janew,list,oldtonew,jnewtojold,bwidth,rhs)
  integer,intent(inout) :: nnodes,nsub,nnz
  integer,intent(inout) :: bwidth

  !real(dp),intent(inout) :: A(nsub,nsub,nnz),rhs(nsub,nnodes)
  !real(dp),pointer,intent(inout),dimension(:,:,:),allocatable :: A
  !real(dp),pointer,intent(inout),dimension(:,:),allocatable ::rhs
  real(dp),pointer,intent(inout),dimension(:,:,:):: A
  real(dp),pointer,intent(inout),dimension(:,:) ::rhs

  !integer,intent(inout) :: ia(nnodes + 1),iau(nnodes),ja(nnz)
  !integer,intent(inout) :: ianew(nnodes + 1),iaunew(nnodes),janew(nnz),jnewtojold(nnz),oldtonew(nnodes),list(nnodes)
  !integer,pointer,intent(inout),dimension(:),allocatable :: ia,iau,ja,ianew,iaunew,janew,jnewtojold,oldtonew,list
  integer,pointer,intent(inout),dimension(:) :: ia,iau,ja,ianew,iaunew,janew,jnewtojold,oldtonew,list

  integer :: nsubprev!ceb test
  integer :: testsize!ceb test
  real(dp)::tmp(5)!ceb test

  integer i,j,k,m,ii,jj,iread,jstart,jend,icol,istat
  real(dp),dimension(:),allocatable:: s!(nsub)
  integer ireorder
  character(len=130),intent(inout) :: matrix_file

  write(6,*) "reading matrix file ",matrix_file
  ireorder=0

  open(UNIT=55,FILE=matrix_file,form='formatted',STATUS='old')
    
  read(55,*)!Dimension
  read(55,*) nnodes
  write(6,*)'2 nnodes=',nnodes
  read(55,*)!block size
  read(55,*) nsub
  read(55,*)!number of non zero entries in matrix
  read(55,*) nnz
  
  !testsize=5!ceb test
  !nsubprev=nsub!ceb test
  !if(nsub>1)nsub=testsize!ceb test
  write(6,*)'2 nsub=',nsub
  write(6,*)'2 nnz=',nnz

  !ceb
  allocate(A(nsub,nsub,nnz),STAT=istat)
  allocate(rhs(nsub,nnodes),STAT=istat)
  allocate(ia(nnodes+1),STAT=istat)
  allocate(iau(nnodes),STAT=istat)
  allocate(ja(nnz),STAT=istat)
  allocate(ianew(nnodes+1),STAT=istat)
  allocate(iaunew(nnodes),STAT=istat)
  allocate(janew(nnz),STAT=istat)
  allocate(jnewtojold(nnz),STAT=istat)
  allocate(oldtonew(nnodes),STAT=istat)
  allocate(list(nnodes),STAT=istat)
  allocate(s(nsub),STAT=istat)
  !ceb

  read(55,*)!IA dimension=nnodes+1
  do  i = 1,nnodes+1
     read(55,*) ia(i)
  end do

  read(55,*)!IAU dimension=nnodes
  do  i = 1,nnodes
     read(55,*) iau(i)
  end do

  read(55,*)!JA nnz=nnz
  do  i = 1,nnz
     read(55,*) ja(i)
!write(6,*)"rf: ja(",i,")=",ja(i)
  end do

  read(55,*)!
  do  i = 1,nnz
     !read in column major
     do j=1,nsub
     !do j=1,nsubprev
        read(55,*) (A(j,k,i),k=1,nsub)
        !read(55,*) (tmp(k),k=1,nsubprev)
        !do k=1,nsubprev
        !   if(j.le.testsize .and. k.le.testsize)A(j,k,i)=tmp(k)
        !end do
     end do
  end do

 close(55)
  
  !ceb hardwired for test case where solution is 1's vector
  ! Write out the matrix and the rhs so the answer is all 1's
  !
  ! First transpose A
  !
  !     call transposeA(nnodes,nnz,ia,ja,iau,A)
  !     call scaleAD(nnodes,nnz,ia,ja,iau,A)
  !
  ! Generate RHS. Use qnode for scratch
  !
  iread = 0
  if(iread.eq.0)then

     do  i = 1,nnodes
        do j=1,nsub
           rhs(j,i) = 0.
        end do
     end do

     do  i = 1,nnodes
        jstart = ia(i) ! Beginning of row
        jend   = ia(i+1) - 1 ! Ending of row
        do k=1,nsub
           s(k)=0.
        end do

        do  j=jstart,jend
           icol   = ja(j) ! We don't care about the column because we have all 1's as the solution
           do k=1,nsub
              do m=1,nsub
                 !s(k) =  s(k) + A(k,m,j)*1.0  
           !rhs(k,i) = rhs(k,i) + A(k,m,j)*1.0
            !rhs(k,i) = rhs(k,i) + A(k,m,j)*1.0
            rhs(k,i) = rhs(k,i) + A(k,m,j)*(icol) !set up so that x(i) == i
!write(6,*) "A[",j,",",k,",",m,"]=", A(k,m,j)
             end do
!write(6,*) "rhs:B(",k,",",i,")=",rhs(k,i) 

           end do
        end do

        !do k=1,nsub
        !   rhs(k,i) = rhs(k,i) + s(k)
        !end do
     end do

     !else
     !read(5,*)
     !do i = 1,nnodes
     !   do k=1,nsub
     ! ! Reverse the sign because we solve Ax+b=0 but this was generated
     ! ! (in converter) as Ax=b
     !    rhs(k,i) = -rhs(k,i)
     !   end do
     !end do


  !do ii=1,nnodes
  !   do jj=1,nsub
  !      write(6,*) "rhs:B(",ii,",",jj,")=",rhs(jj,ii) 
  !   end do
  !end do

  end if
  ! Re-order the matrix using Cuthill-McKee
  call BANDWIDTH(ia,iau,ja,nnodes,nnz,bwidth)
  call reOrder(nnodes,nsub,nnz,ia,ja,iau,ianew,iaunew,janew,list,oldtonew,jnewtojold,ireorder)
  call BANDWIDTH(ianew,iaunew,janew,nnodes,nnz,bwidth)
   
  !call ILUK_factorization(nnodes,nnz,ndim,A,ianew,iaunew,janew)
  deallocate(s,STAT=istat)
  return
end subroutine READERf






!============================ READER_HB ================================
!
! Reads the matrix in Harwell Boeing format and converts to mine
! The only difference is that he reads in nnodes+1 iau's
! This routine also gives the option to transpose 
! or divide by the diagonal
!
!===================================================================
!subroutine READERf(nnodes,nsub,nnz,A,ia,iau,ja,list,tag,ianew,iaunew,janew,oldtonew,rhs)
subroutine READER_HB(matrix_file,nnodes,nsub,nnz,A,ia,iau,ja,ianew,iaunew,janew,list,oldtonew,jnewtojold,bwidth,rhs)
  use hb_file_module
  integer,intent(inout) :: nnodes,nsub,nnz
  integer,intent(inout) :: bwidth

  !real(dp),intent(inout) :: A(nsub,nsub,nnz),rhs(nsub,nnodes)
  !integer,intent(inout) :: ia(nnodes + 1),iau(nnodes),ja(nnz)
  !integer,intent(inout) :: ianew(nnodes + 1),iaunew(nnodes),janew(nnz),jnewtojold(nnz),oldtonew(nnodes),list(nnodes)

  real(dp),intent(inout),pointer,dimension(:,:,:) :: A
  real(dp),intent(inout),pointer,dimension(:,:) :: rhs!(nsub,nnodes)
  integer,intent(inout),pointer,dimension(:) :: ia,iau,ja
  integer,intent(inout),pointer,dimension(:) :: ianew,iaunew,janew,jnewtojold,oldtonew,list

  character(len=130),intent(inout) :: matrix_file
  integer ::i,j,k,m,ii,jj,iread,jrow,colptr_cnt,jstart,jend,jcol,istat,ireorder

  real(dp),allocatable,dimension(:) ::s
  integer,allocatable,dimension(:) ::count,ca

  write(6,*) "reading matrix file ",matrix_file
  ireorder=0
  open(UNIT=55,FILE=matrix_file,form='formatted',STATUS='old')
    
  !read(55,*)!Dimension
  !read(55,*) !nnodes
  !read(55,*)!block size
  !read(55,*) !nsub
  !read(55,*)!number of non zero entries in matrix
  !read(55,*) !nnz
  !  write(6,*)'2 nnodes=',nnodes
  !write(6,*)'2 nsub=',nsub
  !write(6,*)'2 nnz=',nnz


  !read(55,*)!IA dimension=nnodes+1
  !do  i = 1,nnodes+1
  !   read(55,*) ia(i)
  !end do

  !read(55,*)!IAU dimension=nnodes
  !do  i = 1,nnodes
  !   read(55,*) iau(i)
  !end do

  !read(55,*)!JA nnz=nnz
  !do  i = 1,nnz
  !   read(55,*) ja(i)
!!write(6,*)"rf: ja(",i,")=",ja(i)
  !end do

  !read(55,*)!
  !do  i = 1,nnz
  !   !read in column major
  !   do j=1,nsub
  !      read(55,*) (A(j,k,i),k=1,nsub)
  !   end do
  !end do

  !submatrix size
  nsub=1!hardwired for HB matrix unless altered after read

  !HB matrices are in compressed column format.
  !we want compressed row format, so well have to convert.
  !this can be done by transposing the matrix
  !but since the matrix can be nonsymmetric we
  !cant just do a simple element swap
  call hb_file_read(55)
  close(55)


  !Dimension
  nnodes=nrow !ceb we currently assume square matrices
  !number of non zero entries in matrix
  nnz=nnzero
  write(6,*)'HB_read: nnodes=',nnodes
  write(6,*)'HB_read: nsub=',nsub
  write(6,*)'HB_read: nnz=',nnz

  allocate(ia(nnodes+1),STAT=istat)
  allocate(iau(nnodes),STAT=istat)
  allocate(ja(nnz),STAT=istat)
  allocate(A(nsub,nsub,nnz),STAT=istat)

  allocate(count(nnodes),STAT=istat)
  allocate(ca(nnz),STAT=istat)

  !initialize ia
  ia(1)=1
  do i=1,nnodes
     ia(i+1)=0
     iau(i)=-1
     count(i)=0
  end do
  colptr_cnt=1
  do j=1,nnz
     jrow=rowind(j)
     ia(jrow+1)=ia(jrow+1)+1
     if(j.ge.colptr(colptr_cnt+1)) colptr_cnt=colptr_cnt+1
     ca(j)=colptr_cnt
!write(6,*)"ca(",j,")=",ca(j)
  end do
  do i=1,nnodes
     ia(i+1)=ia(i+1)+ia(i)
!write(6,*)"READER_HB: ia(",i+1,")=",ia(i+1)
  end do
  do j=1,nnz
     jrow=rowind(j)
     do ii=1,nsub
        do jj=1,nsub
           A(ii,jj,ia(jrow)+count(jrow))=values(j)
        end do
     end do
     ja(ia(jrow)+count(jrow))=ca(j)
!write(6,*)"ja(",ia(jrow)+count(jrow),")=",ca(j)
     if(jrow==ca(j)) iau(jrow)=ia(jrow)+count(jrow)
     count(jrow)=count(jrow)+1
  end do

!  do j=1,nnz
!write(6,*)"ja(",j,")=",ja(j)
!  end do

  do i=1,nnodes
     do j=ia(i),ia(i+1)-1
        write(31,*)i,ja(j)
     end do
     write(31,*)"iau(",i,")=",iau(i)
  end do

  deallocate(count)
  deallocate(ca)
  !need to create iau
  allocate(rhs(nsub,nnodes),STAT=istat)
  allocate(s(nsub),STAT=istat)



  !ceb hardwired for test case where solution is 1's vector
  ! Write out the matrix and the rhs so the answer is all 1's
  !
  ! First transpose A
  !
  !     call transposeA(nnodes,nnz,ia,ja,iau,A)
  !     call scaleAD(nnodes,nnz,ia,ja,iau,A)
  !
  ! Generate RHS. Use qnode for scratch
  !
  iread = 0
  if(iread.eq.0)then

     do  i = 1,nnodes
        do j=1,nsub
           rhs(j,i) = 0.
        end do
     end do

     do  i = 1,nnodes
        jstart = ia(i) ! Beginning of row
        jend   = ia(i+1) - 1 ! Ending of row
        do k=1,nsub
           s(k)=0.
        end do

        do  j=jstart,jend
           jcol   = ja(j) ! We don't care about the column because we have all 1's as the solution
           do k=1,nsub
              do m=1,nsub
                 !s(k) =  s(k) + A(k,m,j)*1.0  
           !rhs(k,i) = rhs(k,i) + A(k,m,j)*1.0
            !rhs(k,i) = rhs(k,i) + A(k,m,j)*1.0
            rhs(k,i) = rhs(k,i) + A(k,m,j)*(jcol) !set up so that x(i) == i
!write(6,*) "A[",j,",",k,",",m,"]=", A(k,m,j)
             end do
!write(6,*) "rhs:B(",k,",",i,")=",rhs(k,i)
           end do
        end do

        !do k=1,nsub
        !   rhs(k,i) = rhs(k,i) + s(k)
        !end do
     end do

     !else
     !read(5,*)
     !do i = 1,nnodes
     !   do k=1,nsub
     ! ! Reverse the sign because we solve Ax+b=0 but this was generated
     ! ! (in converter) as Ax=b
     !    rhs(k,i) = -rhs(k,i)
     !   end do
     !end do


  !do ii=1,nnodes
  !   do jj=1,nsub
  !      write(6,*) "rhs:B(",ii,",",jj,")=",rhs(jj,ii) 
  !   end do
  !end do
  deallocate(s)


  end if
  ! Re-order the matrix using Cuthill-McKee
  call BANDWIDTH(ia,iau,ja,nnodes,nnz,bwidth)
write(6,*)"READER_HB:checkpt 5"

  allocate(ianew(nnodes+1),STAT=istat)
  allocate(iaunew(nnodes),STAT=istat)
  allocate(oldtonew(nnodes),STAT=istat)
  allocate(list(nnodes),STAT=istat)
  allocate(janew(nnz),STAT=istat)
  allocate(jnewtojold(nnz),STAT=istat)

!write(6,*)"READER_HB:checkpt 5.5"
  call reOrder(nnodes,nsub,nnz,ia,ja,iau,ianew,iaunew,janew,list,oldtonew,jnewtojold,ireorder)
!write(6,*)"READER_HB:checkpt 6"
  call BANDWIDTH(ianew,iaunew,janew,nnodes,nnz,bwidth)
!write(6,*)"READER_HB:checkpt 7"
   
  !call ILUK_factorization(nnodes,nnz,ndim,A,ianew,iaunew,janew)

  return
end subroutine READER_HB



!=================================== GETKA  ==================================
!
! Get ka array so we can quickly go through the columns
!
!=============================================================================
subroutine getKA(nnodes,nnz,ia,iau,ja,ka,kstart)
  integer,intent(in) :: nnodes,nnz
  integer:: i,index
  integer :: ia(nnodes+1),iau(nnodes),ja(nnz),ka(nnz),kstart(nnodes)
  !
  ! First fill kstart, which holds the initial positions in the ka array where
  ! we will put things
  !
  do i = 1,nnodes
     kstart(i) = ia(i)
  end do
  !
  ! Now go through the ja array and 
  !
  do i = 1,nnz
     index = kstart(ja(i))
     ka(index) = i
     kstart(ja(i)) = kstart(ja(i)) + 1
  end do
  return
end subroutine getKA



!=================================== READRESTART ================================
!
! Write out the solution
!
!=============================================================================
subroutine READRESTART(nnodes,nsub,x)
  integer :: nnodes,nsub,i,j,itmp,ierror,istat
  real(dp),intent(inout),pointer,dimension(:,:) :: x!(nsub,nnodes)
  character(len=130)::input_file

  allocate(x(nsub,nnodes),STAT=istat)

  input_file="soln.dat"
  open(UNIT=444,FILE=input_file,form='formatted',STATUS='old',IOSTAT=ierror)

!ceb write logic here to kill old file if format is incorrect ie neqn!=nsub nnodes!=nrows

  if (ierror==0) then
     write(6,*)"reading input from soln.dat"
     do i=1,nnodes
        read(444,*)itmp,(x(j,i),j=1,nsub)
     end do
  else
     write(6,*)"error: failed reading input from soln.dat"
     ! initialize to zero
     do i=1,nnodes
        do j=1,nsub 
           x(j,i)=0.0
        end do
     end do
  end if
  close(444)
write(6,*)"exit readrestart"
  return
end subroutine READRESTART

!=================================== WRITEOUT ================================
!
! Write out the solution
!
!=============================================================================
subroutine WRITEOUT(nnodes,nsub,x)
  integer :: nnodes,nsub,i,j
  real(dp),intent(in) :: x(nsub,nnodes)
  character(len=130)::output_file
  output_file="soln.dat"
  open(UNIT=444,FILE=output_file,form='formatted',STATUS='replace')
  do i=1,nnodes
     write(444,*)i,(x(j,i),j=1,nsub)
  end do
  close(444)
  return
end subroutine WRITEOUT



!=================================== DQX =====================================
!
! Copy dq into the 'X' array for GMRES
! The reason we have to go to this trouble is because we have to renumber 
!
!=============================================================================
subroutine DQX(nnodes,dq,X,nsub)
  integer, intent(in) :: nnodes,nsub
  real(dp),intent(in) :: dq(nsub,nnodes)
  real(dp),intent(out) :: X(nnodes*nsub)
  integer :: i,m,index
  do i = 1,nnodes
     index = nsub*i - (nsub - 1)
     do m = 1,nsub
        X(index) = dq(m,i)
        index = index + 1
     end do
  end do
  return
end subroutine DQX


!=================================== XDQ =====================================
!
! Copy 'X' back into dq 
! The reason we have to go to this trouble is because we have to renumber 
!
!=============================================================================
subroutine XDQ(nnodes,dq,X,nsub)
  integer, intent(in) :: nnodes,nsub
  real(dp) dq(nsub,nnodes),X(nnodes*nsub)
  integer :: i,m,index
  do  i = 1,nnodes
     !index = nsub*i - (nsub - 1)
     index = nsub*(i-1)
     do m = 1,nsub
        index = index + 1
        dq(m,i) = X(index) 
     end do
  end do
  return
end subroutine XDQ



!=================================== QTOQ =====================================
!
! Simply copy X into dq where both arrays are 2D (nsub,nnodes)
!
!=============================================================================
subroutine QTOQ(nnodes,dq,X,nsub)
  integer, intent(in) :: nnodes,nsub
  integer :: i,m
  real(dp) dq(nsub,nnodes),X(nsub,nnodes)
  do  i = 1,nnodes
     do m = 1,nsub
        dq(m,i) = X(m,i) 
     end do
  end do
  return
end subroutine QTOQ


!================================ ILUPRE ====================================
!
! Preconditioning using ILU
!
!============================================================================
!subroutine ILUPRE(nnodes,nnz,nsub,idolu,dq,A,ALU,ia,ja,iau,phi,iprecon,mypre,nnz_pc,ia_pc,ja_pc,iau_pc,fill,droptol)
subroutine ILUPRE(nnodes,nnz,nsub,idolu,A,ALU,ia,ja,iau,phi,mypre,nnz_pc,ia_pc,ja_pc,iau_pc,fill,droptol)
  integer :: nnodes,nnz,nsub,idolu,mypre,istat
  real(dp), intent(inout), allocatable, dimension (:,:,:):: ALU
  real(dp), intent(in), dimension(nsub,nsub,nnz) :: A
  integer,  intent(in), dimension(nnodes+1) :: ia
  integer,  intent(in), dimension(nnz) :: ja
  integer,  intent(in), dimension(nnodes) :: iau

  integer, intent(inout) :: nnz_pc
  integer, intent(inout), allocatable, dimension(:) :: ia_pc
  integer, intent(inout), allocatable, dimension(:) :: ja_pc
  integer, intent(inout), allocatable, dimension(:) :: iau_pc

  integer,  allocatable, dimension(:) :: ia_re
  integer,  allocatable, dimension(:) :: iau_re
  integer,  allocatable, dimension(:) :: ja_re
  integer ::ind(nnodes)

  real(dp),intent(inout), dimension(nsub,nnodes) :: phi
  !real(dp),intent(inout), dimension(nsub,nnodes) :: dq
  integer :: i,j,k

  !ilut and schur complement vars
  integer :: schur_start
  integer :: fill ! max fill level
  real(dp), intent(in) :: droptol !dropping tolerance

  integer jstart,jend,jcol,pos,tmpj

  integer bsize !min size of indep block
  integer nbnd !num interior vars
  integer ind_nnodes !number of elements in independent set
  integer iord(nnodes) !map old row index to new reordered row index (permutation array)
  integer riord(nnodes) !map old row index to new reordered row index (permutation array)
  real(dp) indtol !tolerance for excluding row from independent set (diagonal dominance value)
 
 
  ! First we must do the ILU decompositions so lets assemble the A matrix
  select case(mypre)
  case(1)
     write(6,*)"preconditioner is diag(A)"
  case(2)
     write(6,*)"preconditioner is GS"
  case(3)
     write(6,*)"preconditioner is bkaczmarz"
  case(4:5)
     if(mypre.eq.4) fill=0
     write(6,*)"preconditioner is ILU(",fill,")" 
     allocate(ALU(nsub,nsub,nnz),STAT=istat)
     ! First copy A into ALU 
     do  i = 1,nnz
        do  j = 1,nsub
           do  k = 1,nsub
              ALU(j,k,i) = A(j,k,i)
           end do
        end do
     end do
     call BLKILU_0(nnodes, nsub, nnz, ia, ja, iau, ALU)
  case (6)
     write(6,*)"preconditioner is ILUT" 
     write(6,*) "pre blkilut"
     call BLKILUT(nnodes,nsub,nnz,ia,ja,iau,A,ALU,droptol,fill,nnz_pc,ia_pc,ja_pc,iau_pc)
     write(6,*) "pst blkilut"
  case(7)!ARMS
     write(6,*)"preconditioner is ARMS"
  case default
     write(6,*)"no preconditioner"       
  end select
  
  return
end subroutine ILUPRE


!======================================================
subroutine diag_invpc(iau,A,x,b,nnodes,ndim,nnz)
  integer, intent(in) :: nnodes,ndim,nnz
  integer ::  i,k,idiag
  integer, intent(in) :: iau(nnodes)
  real(dp), intent(in), dimension (ndim,ndim,nnz) :: A
  real(dp), intent(inout) :: x(ndim,nnodes)
  real(dp), intent(in) :: b(ndim,nnodes)
  real(dp) :: LU(ndim,ndim)
  real(dp) :: rhs(ndim)
  do i=1,nnodes
     idiag=iau(i)
     do k=1,ndim
        rhs(k)=b(k,i)
     end do
     call get_submat(A,LU,idiag,nnz,ndim)
     !ceb test by setting LU to Identity
     !call blank_2D_mat(LU,ndim,ndim)
     !do k=1,ndim
     !   LU(k,k)=1.0
     !enddo
     !ceb end test
     call LUgeneral_rm(LU,ndim)
     !solving  (LU)x = b   array for x
     call BACKSUBgeneral(LU,rhs,ndim)
     do k=1,ndim
        x(k,i)=rhs(k)
     end do
  end do
  return
end subroutine diag_invpc
!======================================================


!=======================================================================
subroutine BLK_LU_INV_SOLVE(dim, neqn, nnz, y, x, alu, ia, ja, iau)
  !//!! ja[k] must not be a phantom or ghost node!! 
  !// actually this function and blkilu need to properly handle phantom indices!!
  integer, intent(in) :: dim,neqn,nnz
  integer, intent(in), dimension(dim+1) :: ia
  integer, intent(in), dimension(nnz) ::   ja
  integer, intent(in), dimension(dim) ::   iau
  real(dp), intent(in), dimension (neqn,neqn,nnz) :: ALU
  real(dp), intent(inout),dimension(neqn,dim) :: x
  real(dp), intent(in), dimension(neqn,dim) :: y
  integer :: i,j,jj,ii,m1,m2,p,jrow,diag
  real(dp):: t(neqn)
  !//c----------------------------------------------------------
  !//c solve LUx=y for x
  !//c----------------------------------------------------------
  !L[i,j]*y[j]=b[i](off diagonals only)
  !U[i,j]*x[j]=y[i](off diagonals only)
  !D[i]*x[i]= 

  !//compute y[i]= rhs[i]-L[p,j]*y[j]
  do p=1,dim
     m1=ia(p)
     m2=iau(p)-1
     do ii=1,neqn
        x(ii,p)=y(ii,p)!//x=rhs
        do j=m1,m2!//Lower Matrix
           do jj=1,neqn
              x(ii,p) = x(ii,p) - alu(ii,jj,j)*x(jj,ja(j))
           end do
        end do
     end do
  end do

  !//compute x[p]= y[p]-U[p,j]*x[j]
  do p=dim,1,-1
     diag=iau(p)
     m1=diag+1
     m2=ia(p+1)-1
     do ii=1,neqn
        do j=m1,m2!//Upper Matrix
           do jj=1,neqn
              x(ii,p) = x(ii,p) - alu(ii,jj,j)*x(jj,ja(j))
           end do
        end do
     end do

     !//x[p]=D[p]*x[p]
     do ii=1,neqn
        t(ii)=0.0
        do jj=1,neqn		
           t(ii) = t(ii) + alu(ii,jj,diag)*x(jj,p)
        end do
     end do
     do ii=1,neqn
        x(ii,p)=t(ii)
     end do
  end do
  return
end subroutine BLK_LU_INV_SOLVE
!===================================================================

!=======================================================================
!subroutine BLK_LU_INV_SOLVE_form2(dim, neqn, nnz, y, x, alu, ia, ja, iau)
subroutine BLK_LU_INV_SOLVE_form2(dim, neqn, y, x, L,U)
  !//!! ja[k] must not be a phantom or ghost node!! 
  !// actually this function and blkilu need to properly handle phantom indices!!
  integer, intent(in) :: dim,neqn
  !integer, intent(in), dimension(dim+1) :: ia
  !integer, intent(in), dimension(nnz) ::   ja
  !integer, intent(in), dimension(dim) ::   iau
  !real(dp), intent(in), dimension (neqn,neqn,nnz) :: ALU

  !integer, intent(in), dimension(:) :: U_rownnz
  !type(ja_type), intent(in), dimension(:) :: U_ja
  !type(submat_type), intent(in), dimension(:) :: U
  !integer, intent(in), dimension(:) :: L_rownnz
  !type(ja_type), intent(in), dimension(:) :: L_ja
  !type(submat_type), intent(in), dimension(:) :: L
  type(matrix_type),intent(in) :: U
  type(matrix_type),intent(in) :: L


  real(dp), intent(inout),dimension(neqn,dim) :: x
  real(dp), intent(in), dimension(neqn,dim) :: y
  integer :: i,j,jj,ii,m1,m2,p,jrow,diag,jnode
  real(dp):: t(neqn)
  !//c----------------------------------------------------------
  !//c solve LUx=y for x
  !//c----------------------------------------------------------
  !L[i,j]*y[j]=b[i](off diagonals only)
  !U[i,j]*x[j]=y[i](off diagonals only)
  !D[i]*x[i]= 

  !//compute y[i]= rhs[i]-L[p,j]*y[j]
  do i=1,dim 
     do j=1,L%rownnz(i)
        jnode=L%ja(i)%cols(j)
        do ii=1,neqn
           x(ii,i)=y(ii,i)!//x=rhs
           do jj=1,neqn
              x(ii,i)= x(ii,i) - L%row(i)%submat(ii,jj,j)*x(jj,jnode)
           end do
        end do
     end do
  end do

  !//compute x[p]= y[p]-U[p,j]*x[j]
  do i=dim,1,-1
     do j=1,U%rownnz(i)
        jnode=U%ja(i)%cols(j)
        do ii=1,neqn
           do jj=1,neqn
              x(ii,i)= x(ii,i) - U%row(i)%submat(ii,jj,j)*x(jj,jnode)
           end do
        end do
     end do

     do ii=1,neqn
        t(ii)=0.0
        do jj=1,neqn
           t(ii)= t(ii) - U%row(i)%submat(ii,jj,1)*x(jj,i)
        end do
     end do

     do ii=1,neqn
        x(ii,i)=t(ii)
     end do
  end do

  return
end subroutine BLK_LU_INV_SOLVE_form2
!===================================================================




!===================================================================
subroutine BLKILU_0(dim, neqn, nnz, ia, ja, iau, alu)
  integer :: dim,neqn,nnz
  integer, intent(in), dimension(dim+1) :: ia
  integer, intent(in), dimension(nnz) :: ja
  integer, intent(in), dimension(dim) :: iau
  real(dp), intent(inout), dimension(neqn,neqn,nnz):: ALU
  integer :: j,k,L1,L2,ii,jj,kk,jcol,j1,j2,jw,p,diagj,diagp
  real(dp) :: pivot(neqn,neqn)
  real(dp) :: Lij(neqn,neqn)
  real(dp) :: Uij(neqn,neqn)
  real(dp) :: Aij(neqn,neqn)
  real(dp) :: D(neqn)
  real(dp) :: temp  
  integer :: iw(dim+1)!//this vector points to the column in alu corresponding to a given node index
  !write(6,*) "dim+1=",dim+1
  !//initialize column indexing array entries to -1
  do jj=1,dim
     iw(jj)=-1
     !write(6,*)"iw(",jj,")=",iw(jj)
  end do
  
  do p=1,dim!//loop over block rows
     !{
!write(6,*)"iw checkpt line 692 p=",p
     j1=ia(p)
     diagp=iau(p)!diagonal index of rowp in A
     do jj=j1,ia(p+1)-1
        iw(ja(jj))=jj
     end do!//for current row p, record j column indices into alu

     !==============================================================================
     do j=j1,diagp-1!//loop over column block entries in row p, left of diagonal 
        !{
        jcol = ja(j)!get node index for column j in row p   
        ! get index of diagonal for node corresponding to column j in row p
        diagj=iau(jcol)
        ! set current block entry A(p,j)= A(p,j)*A(diag(rowj))

        !   Compute the block multiplier (pivot) for jrow.
        !   Multiply alu(rowp,colj) by alu(diag(colj)) (mult row p col j entry (submat) by diagonal entry of col j above row p) !left of diagonal only
        !call blank_2D_mat(pivot,neqn,neqn)
        do L1=1,neqn
           do L2=1,neqn
              pivot(L1,L2) = 0
              do kk=1,neqn
                 pivot(L1,L2) = pivot(L1,L2) + alu(L1,kk,j)*alu(kk,L2,diagj) ! T = A(p,j)*A(j,j) = L(p,j)*D(j)
              end do
              !write(6,*)"[",p,jcol,"] pivot(",L1,L2,")=",pivot(L1,L2)
           end do
        end do 

        !ceb test
        !if (L1Norm(pivot,neqn) <= 1.0e-1) then !skip to end of loop over rows in L
        !   write(6,*)"skip block col ",jcol," (dropping) and go to next"
        !   cycle!continue
        !endif
        !ceb test

        !write(6,*)"L1(pivot([",p,jcol,"]=",L1Norm(pivot,neqn)

        !//perform linear combinations on block columns (jrow) right of diagonal 
        ! row operation: jrow = jrow - jrow*pivot !performed for lower matrix only ??
        do jj= diagj+1,ia(jcol+1)-1! //block rows right of diagonal and above row p
           jw=iw(ja(jj))!//get block column index for nonzero node in this jrow
           if(jw .ne. -1) then
              !//loop over submatrix entry
              do L1=1,neqn
                 do L2=1,neqn
                    do kk=1,neqn
                       alu(L1,L2,jw) = alu(L1,L2,jw) - pivot(L1,kk)*alu(kk,L2,jj) ! U(j,jw)=U(j,jw)-Pvt*U(j,jj) where jj>j
                    end do
                    !write(6,*)"U[",jcol,ja(jw),"](",L1,L2,")=",alu(L1,L2,jw)
                 end do
              end do
           end if
        end do

        ! Replace L(p,j) with A(p,j)*A(j,j) (set col entry to (col entry mult diag))
        do jj=1,neqn
           do kk=1,neqn
              alu(jj,kk,j)=pivot(jj,kk) ! L(p,j)= A(p,j)*A(j,j)=Pvt = L(p,j)*D(j) where j<p
           end do
        end do
        !}
     end do !end j loop over block columns (left of diag) of block row p
     !==============================================================================

     if (ja(diagp) .ne. p) then!//if (jrow .ne. k) goto 600
        write(6,*)"Zero pivot found on row ",p
        return
     end if






     
     !//define and store lower part ====================================
     !//compute L,U for pivot (diagonal of combined LU matrix)
     ! Lij,Uij submatrices refer only to the diagonal block element
     temp=0.0
     do jj=1,neqn
        do kk=1,neqn
           if (jj.eq.kk) then!//set diagonal (sub diagonal)
              temp=alu(jj,jj,diagp)!A[p,p](jj,jj)
              do ii=1,jj-1
                 temp = temp - (Lij(ii,jj)*Uij(jj,ii)) ! A[p,p](jj,jj)=A[p,p](jj,jj)- L[pj](ii,jj)*U[pj](jj,ii)
              end do
              Lij(kk,kk)= 1.0/temp 
           else if(kk<jj) then!//set lower
              temp = alu(jj,kk,diagp)!A[p,p](jj,kk) L 
              do ii=1,kk-1
                 temp = temp - Lij(ii,jj)*Uij(kk,ii)
              end do
              Lij(kk,jj)=temp
           else if(kk>jj) then !//set upper
              temp = alu(jj,kk,diagp)!A[p,p](jj,kk) U
              do ii=1,jj-1
                 temp = temp - Lij(ii,jj)*Uij(kk,ii)
              end do
              Uij(kk,jj) = temp*Lij(jj,jj)
           end if
        end do
     end do

     !//set alu[diagp] to identity
     do jj=1,neqn
        do kk=1,neqn
           alu(jj,kk,diagp)=0.0
        end do
        alu(jj,jj,diagp)=1.0
     end do

     !//write LU into alu[diagp]
     do kk=1,neqn !loop subcols
        do jj=1,neqn-1 !get D[1,..,neqn-1] ! loop subcols
           D(jj)=alu(jj,kk,diagp) !Identity matrix row jj 
           do ii=1,jj-1
              D(jj) = D(jj) - Lij(ii,jj)*D(ii)
           end do
           D(jj) = D(jj)*Lij(jj,jj)
        end do
              
        do jj=neqn,1,-1 !get alu[diagp,jj,kk] ! loop subcols
           if(jj.eq.neqn) then     
              !D=D-L*D (left of subdiag)
              do ii=1,jj-1
                 alu(jj,kk,diagp) = alu(jj,kk,diagp) - Lij(ii,jj)*D(ii)
              end do
              alu(jj,kk,diagp) = alu(jj,kk,diagp)*Lij(jj,jj)
              !D=
           else      
              !D=D-U*D (right of sub diag)
              alu(jj,kk,diagp)=D(jj)
              do ii=jj+1,neqn
                 alu(jj,kk,diagp) = alu(jj,kk,diagp) - Uij(ii,jj)*alu(ii,kk,diagp)
              end do
           end if
           !write(6,*)"A(",p,",diag)[",jj,kk,"]=",Aij(jj,kk)
           !write(6,*)"U(",p,",diag)[",jj,kk,"]=",alu(jj,kk,diagp)
        end do
     end do

     !reset iw
     !write(6,*) "p=",p," jjstart=",ia(p)," jjend=",ia(p+1)-1
     do jj = ia(p), ia(p+1)-1
        !write(6,*)"jj",jj
        !write(6,*)"ja(",jj,")=",ja(jj)
        iw(ja(jj)) = -1
        !write(6,*)"iw(",ja(jj),")=",iw(ja(jj))
     end do
     !write(6,*)"checkpt iw test"
     !}
  end do!end loop over rows
  
  !write(6,*)"checkpt post iw test"

  !write(6,*)"ilu:"
  !do p=1,dim
  !   !do ii=1,neqn
  !   !   do jj=1,neqn   
  !   !      write(6,*) (alu(ii,jj,j),j=ia(p),ia(p+1)-1)
  !   !   end do
  !   !end do
  !   do ii=1,neqn
  !      write(6,*) (alu(ii,ii,iau(p)))
  !   end do

  !write(6,*)"row:",p
  !     do j=ia(p),ia(p+1)-1
  !        write(6,*)"ja(",j,")=",ja(j)
  !     end do
  !end do

  !do p=1,dim
    !write(6,*)"row:",p
  !     do ii=1,neqn
  !        do jj=1,neqn   
  !           write(6,*) (alu(ii,jj,j),j=ia(p),ia(p+1)-1)
  !        end do
  !     end do
  !end do

  return
end subroutine BLKILU_0












!=================================== MATMUL ==================================
!
! Multiply by A a few times to see if we can start with an eigenvector
!
!=============================================================================
subroutine MATMUL(nnodes,nnz,A,ia,iau,ja,x,B,nsub)
  integer, intent(in) :: nnodes,nnz,nsub
  real(dp),intent(in) :: A(nsub,nsub,nnz)
  real(dp),intent(inout) :: x(nsub,nnodes)
  real(dp),intent(in) :: B(nsub,nnodes)
  real(dp) :: tres(nsub,nnodes)
  integer, intent(in) :: ia(nnodes+1),iau(nnodes),ja(nnz)
  integer i,j,k,m,ii,maxiter,icol,idiag,istart,iend,iter
  
  maxiter = 20
  do  i = 1,nnodes
     do m = 1,nsub
        tres(m,i) = x(m,i)
        x(m,i) = 0.
     end do
  end do
  
  do  iter = 1,maxiter
     do  i = 1,nnodes
        istart = ia(i) ! Row of matrix starts here
        iend   = ia(i+1) - 1 ! Row of matrix ends here
        idiag  = iau(i)      ! Diagonal
        do  ii = istart,iend
           icol = ja(ii)
           do  j=1,nsub
              do  k=1,nsub
                 x(j,i) = x(j,i) + A(j,k,ii)*tres(k,icol)
              end do
           end do
        end do
     end do
     
     do  i = 1,nnodes
        do m = 1,nsub
           tres(m,i) = tres(m,i) + x(m,i) + B(m,i)
           x(m,i) = tres(m,i)
        end do
     end do
  end do
  
  return
end subroutine MATMUL
!==========================================================



!================================= GLU =======================================
!
! LU decomposition of an neqn x neqn matrix
! Here, we loop through the A matrix and compute the LU decomposition
! of the diagonal. 
!
!=============================================================================
subroutine DLU(D,nnodes,nsub)
  integer,intent(in) :: nnodes,nsub
  integer ii
  real(dp),intent(inout) :: D(nsub,nsub,nnodes)
  integer segment,row,column
  real(dp) Dtmp(nsub,nsub)
  ! Now do LU decomposition
  do  ii = 1,nnodes
     call get_submat(D,Dtmp,ii,nnodes,nsub)
     call LUgeneral_rm(Dtmp,nsub)
     call set_submat(D,Dtmp,ii,nnodes,nsub)
  end do
  return
end subroutine DLU

!=============================================================================
subroutine BACKSUBgeneral(a,b,n)
  integer,intent(in) :: n
  real(dp),intent(in) :: a(n,n)
  real(dp),intent(inout) :: b(n)
  integer i,j
  real(dp) sum  
  !// do the forward substitution on Ly = b.  The b then
  !// is formed into "y", which is to be used in the
  !// backward substitution step.
  !// note that the "L" matrix has ones along the diagonal
  b(1) = b(1)/a(1,1)
  do i = 2,n !  // accumulate the rhs term into "sum"
     sum = b(i)
     do j=1,i-1
        sum = sum - a(i,j)*b(j)
     end do
     b(i) = sum/a(i,i)
  end do
  !// now do the backward subsitution on Ux = y.  The
  !// b then is formed into x, the solution to the
  !// system.
  do i=n-1,1,-1
     sum = b(i)
     do j = n,i+1,-1
        sum = sum - a(i,j)*b(j)
     end do
     b(i) = sum
  end do
  return
end subroutine BACKSUBgeneral


!================================= GBACKSUB ==================================
!
! Performs backsubstitution (general neqn x neqn)
!
!=============================================================================
subroutine GBACKSUB(nnodes,A,D,q,rhs,ia,ja,iau,nnz,niter,iwriteres,nsub)
  integer :: nnodes,nnz,niter,iwriteres,nsub
  !real(dp),intent(in) :: A(nsub,nsub,nnodes)
  real(dp),intent(in) :: A(nsub,nsub,nnz)
  real(dp),intent(in) :: D(nsub,nsub,nnodes)
  real(dp),intent(in) :: rhs(nsub,nnodes)
  real(dp),intent(inout) ::q(nsub,nnodes)

  real(dp):: resmonitor(25),rms
  integer :: ia(nnodes+1),ja(nnz),iau(nnodes)
  integer :: i,j,k,ii,istart,iend,idiag,icol

  real(dp) :: Btmp(nsub)
  real(dp) :: Dtmp(nsub,nsub)
  !real(dp) :: res(nsub,nnodes)
  ! loop over the nodes and add the offdiagonal pieces
  rms = 0.0

  do  i = 1,nnodes
     do j = 1,nsub
        Btmp(j)= rhs(j,i)
     end do

     istart = ia(i)       ! Row of matrix starts here
     iend   = ia(i+1) - 1 ! Row of matrix ends here
     idiag  = iau(i)      ! Diagonal
!write(6,*)"ia(",i,")=",istart," ia(",i+1,")-1=",iend
     ! Loop over the columns in the matrix
     do  ii = istart,iend
        icol = ja(ii)
!write(6,*)"GBACKSUB: ja(",ii,")=",icol
        if(ii.ne.idiag)then!skip diagonal
           do j=1,nsub
              do k = 1,nsub
                 Btmp(j) = Btmp(j) - A(j,k,ii)*q(k,icol)
              end do
           end do
        end if
     end do

     call get_submat(D,Dtmp,i,nnodes,nsub)
     call BACKSUBgeneral(Dtmp,Btmp,nsub)
     do j = 1,nsub
        q(j,i) = Btmp(j)
     end do
  end do
  !call findrms_prll(pcomm,ia,ja,A,dq,rhs,res,rms,nnodes,nsub,nnz)
  !call findrms(ia,ja,A,dq,rhs,res,rms,nnodes,nsub,nnz)
  !write(6,*)"In backsub n = ",niter," rms = ",rms
  return
end subroutine GBACKSUB




!=================================== SOLVE ===================================
!
! Solve the matrix problem
! A =   Matrix
! B =   RHS
! D =   Diagonal of A
! X =  Solution
! rhs = Right-hand-side
!
!=============================================================================
subroutine SOLVE(nnodes,nsub,nnz,ia,ja,iau,A,X,B,nsubiters,iwriteres)
  integer, intent(in) :: nnodes,nsub,nnz,nsubiters,iwriteres
  real(dp),intent(in), dimension(nsub,nsub,nnz) ::  A
  real(dp),intent(inout),dimension(nsub,nnodes) :: X
  real(dp), intent(in), dimension(nsub,nnodes) :: B
  integer, intent(in), dimension(nnodes+1) :: ia
  integer, intent(in), dimension(nnodes) :: iau
  integer, intent(in), dimension(nnz) :: ja
  integer :: i,idiagonal,j,k,n
  real(dp), dimension(nsub,nsub,nnodes) :: D
  !real(dp),intent(in), dimension(nsub,nnodes) ::  rhs
  !real(dp), dimension(nsub,nnodes) :: res

  ! Initialize dq to zero
  !call blank_2D_mat(X,nsub,nnodes)! dq is axed in backsub anyway

  ! If doing the LU decomposition save the diagonal and do the lu
  do  i = 1,nnodes
     idiagonal = iau(i)
     do  j = 1,nsub
        do  k = 1,nsub
           D(j,k,i) = A(j,k,idiagonal)
        end do
     end do
  end do
  
  ! LU decomposition on the diagonal
  ! This writes over D, that is why we save the diagonal earlier
  call DLU(D,nnodes,nsub) ! This just inverts the diagonal

  !call findrms_prll(pcomm,ia,ja,A,X,B,res,rms,nnodes,nsub,nnz)
  !write(6,*)"In solve initial rms = ",rms

  do  n = 1,nsubiters!        Update delta Q at this iteration
     call GBACKSUB(nnodes,A,D,X,B,ia,ja,iau,nnz,n,iwriteres,nsub)
  end do

  return
end subroutine SOLVE





!================================= LU =======================================
!
! LU decomposition of a 4 x 4 matrix
! Here, we loop through the A matrix and compute the LU decomposition
! of the diagonal. 
!
!=============================================================================
subroutine LU(nnodes,ndim,D,nnz,iau,neqn)
  integer, intent(in) :: nnodes,ndim,nnz,neqn
  real(dp)  D(neqn,neqn,nnodes)
  integer iau(nnodes)
  integer i,j,k,ii,jj,kk
  real(dp) sum

  ! Now do LU decomposition
  do i = 1,nnodes    
     do jj=1,ndim
        ! procedure one, this solves for betaij (upper)
        do ii=1,jj-1
           sum=D(k,j,i)
           do kk=1,ii-1
              sum = sum - D(ii,kk,i)*D(kk,jj,i)
           end do
           D(ii,jj,i)=sum
        end do      
        !procedure two, this solves for alphaij (lower)
        do ii=jj,ndim
           sum=D(ii,jj,i)
           do k=1,jj-1
              sum = sum - D(ii,kk,i)*D(kk,jj,i)
           end do
           D(ii,jj,i)=sum
        end do
        
        sum=1.0/D(jj,jj,i)
        do ii=jj+1,ndim
           D(ii,jj,i) =  D(ii,jj,i)*sum
        end do
     end do
  end do

  return
end subroutine LU


!=============================================================================
subroutine compute_residual(nnodes,ndim,nnz,ia,ja,A,x,B,res,rms)
  integer,intent(in) :: nnodes,ndim,nnz
  real(dp),intent(in) :: A(ndim,ndim,nnodes),B(ndim,nnodes)
  real(dp),intent(inout) :: res(ndim,nnodes)
  real(dp),intent(in) :: x(ndim,nnodes)
  integer,intent(in) :: ia(nnodes),ja(nnz)

  !real(dp),intent(in),pointer,dimension(:,:,:) :: A
  !real(dp),intent(in),pointer,dimension(:,:) :: B
  !real(dp),intent(inout),pointer,dimension(:,:) :: res
  !real(dp),intent(in),pointer,dimension(:,:) :: x
  !integer,intent(in),pointer,dimension(:) :: ia,ja
  integer i,j,k,m
  real(dp),intent(inout) :: rms

  integer ii,jj,kk,jstart,jend

!write(6,*)"compute_resid checkpt 1"

  !do ii=1,nnodes
  !   do jj=1,nsub
  !      write(6,*) "rhs:B(",ii,",",jj,")=",B(jj,ii) 
  !   end do
  !end do
  !do i=1,nnodes
  !   do jj=1,ndim
  !      write(6,*) "x(",ii,",",jj,")=",x(jj,i) 
  !   end do
  !end do

  !do i=1,nnodes
  !   do jj=1,ndim
  !      write(6,*) "rhs:B(",ii,",",jj,")=",B(jj,i) 
  !   end do
  !end do

!write(6,*)"compute_resid checkpt 5"

  rms=0.
  ! Compute A*dq(current) = A*M(inv)*dq(original)
  call matvect(ia,ja, A, x, res, nnodes, ndim, nnz,1)

!write(6,*)"compute_resid checkpt 5"

!do i=1,nnodes
!   jstart=ia(i)
!   jend=ia(i+1)-1
!   do j=jstart,jend
!      do ii=1,ndim
!         res(ii,i)=0
!         do jj=1,ndim
!            res(ii,i)=res(ii,i)+A(ii,jj,j)*x(ii,i)
!         end do
!      end do
!   end do
!end do

  !do i=1,nnodes
  !   do jj=1,ndim
  !      write(6,*) "res(",jj,",",i,")=",res(jj,i) 
  !   end do
  !end do
  
  do  i = 1,nnodes
     do m = 1,ndim
       res(m,i) = B(m,i) - res(m,i)
       rms=rms+res(m,i)*res(m,i)
     end do
  end do
  rms=sqrt(rms)
  !write(6,*) "residual rms= ",rms 
  return
end subroutine compute_residual

!=============================================================================


!//================================================================
subroutine cond_inf_LU(ia,ja,iau,A,LU,b,dim,neqn,nnz,mypre)

  integer :: neqn,dim,nnz,mypre
  integer ,intent(in), dimension(dim+1) :: ia
  integer ,intent(in), dimension(dim) :: iau
  integer ,intent(in), dimension(nnz) :: ja
  real(dp),intent(in), dimension(neqn,neqn,nnz) :: A
  real(dp),intent(inout), dimension(neqn,neqn,nnz) :: LU

  integer :: nnz_pc
  integer,  allocatable,dimension(:) :: ia_pc
  integer,  allocatable,dimension(:) :: ja_pc
  integer,  allocatable,dimension(:) :: iau_pc

  integer i,j,k,ii,jj,kk
  integer maxind,minind,diag
  real(dp) :: rhs(neqn,dim)  
  real(dp) :: xtemp(neqn,dim)
  real(dp) :: condition_num(neqn)
  real(dp) :: invLU(neqn,neqn)
  real(dp) :: LUtmp(neqn,neqn)
  real(dp) max_cond,maxelem,minelem,norm,maxpivot,pivot

 real(dp) :: res(neqn,dim)
 real(dp),intent(in) :: b(neqn,dim)
 real(dp) rms

  do i=1,dim
     do j=1,neqn
        rhs(j,i)=1.0
     end do
  end do


  !call findrms_prll(pcomm,ia,ja, A, rhs, b, res, rms,dim,neqn,nnz)
  call findrms(ia,ja, A, rhs, b, res, rms,dim,neqn,nnz)
  write(6,*) "cond_inf_LU: rms=",rms

  !//solve for x: LUx=1.0
  !LU should already exist here
  call BLK_LU_INV_SOLVE(dim, neqn, nnz, rhs, xtemp, LU, ia, ja, iau)
  !//compute rms of xtemp (infinity norm)
  do j=1,neqn
     condition_num(j)=0.0
     do k=1,dim   
        if ( abs(xtemp(j,k)) > condition_num(j) ) then
           condition_num(j) = abs(xtemp(j,k))
           !if(isnan(real(xtemp[k*neqn+j]))){condition_num[j]=std::numeric_limits<TheType>::max();k=dim;}
        end if
     end do
  end do

  max_cond=0.0
  do j=1,neqn
     write(6,*) "condition_num[",j,"]=",condition_num(j)
      if(max_cond < condition_num(j)) then
         max_cond=condition_num(j)
      end if
   end do

  write(6,*) "condition number = ",max_cond

  !//find max element in LU
  maxelem=0
  minelem=1.0e32
  do i=1,nnz
     call get_submat(LU,LUtmp,i,nnz,neqn)
     norm=FNorm(LUtmp,neqn)
     if (norm>0) then
        if(norm>maxelem) then
           maxelem=norm
           maxind=i
        end if
        if(norm<minelem) then
           minelem=norm
           minind=i
        end if
     end if
  end do

  write(6,*)"Min LU[",minind,"] elem = ",minelem
  write(6,*)"Min LU[",maxind,"] elem = ",maxelem

  !//find max 1/pivot
  maxpivot=0
  call blank_2D_mat(invLU,neqn,neqn)
  do k=1,dim
     diag=iau(k)
     call get_submat(LU,LUtmp,diag,nnz,neqn)
     call InvertSubMatrx(LUtmp,invLU,neqn)
     pivot = FNorm(invLU,neqn)
!write(6,*) "pivot=",pivot
     if(maxpivot<pivot) maxpivot=pivot
  end do

  write(6,*) "Max 1/pivot =",maxpivot
end subroutine cond_inf_LU







!//================================================================
!/* This function sets up the residual vector and finds the rms, given the ia and ja vectors, the A vector of [BLOCKSIZE][BLOCKSIZE] matrices, the solution vector x of dimension [dim][BLOCKSIZE], the rhs vector b of dimension [dim][BLOCKSIZE] the residual vector of dimension [dim*BLOCKSIZE], the dimension of the matrix dim, and the rms value rms */
subroutine findrms_prll(pcomm,ia,ja, A, x, b, res, rms,dim,neqn,nnz)
  !type(PComm_type),intent(inout)::pcomm
  integer::pcomm
  integer :: dim,neqn,nnz
  real(dp) rms,grms,sum
  integer ,intent(in), dimension(dim+1) :: ia
  integer ,intent(in), dimension(nnz) :: ja
  real(dp),intent(in), dimension(neqn,neqn,nnz) :: A
  real(dp),intent(inout)  :: res(neqn,dim)
  real(dp),intent(inout) :: x(neqn,dim)
  real(dp),intent(in)  :: b(neqn,dim)
  integer :: j, k,ii, jj, jstart, jend,icol,ierr
  rms = 0.0
  do k=1,dim	
     do jj=1,neqn
        res(jj,k)=b(jj,k)
     end do
     !ceb needs a parallel matvec multiply
     jstart = ia(k)
     jend = ia(k+1)-1
     do j=jstart,jend
        icol=ja(j)
        do ii = 1,neqn
           do jj=1,neqn
              res(ii,k) = res(ii,k)-A(ii,jj,j)*x(jj,icol)
           end do
        end do
     end do
!write(6,*)"findrms: res(",1,k,")=",res(1,k)
     do j=1,neqn
        rms=rms+ res(j,k)*res(j,k)
     end do
  end do
  !ierr=GetScalarSum_r(pcomm,rms,grms)
  rms = sqrt(grms)
  return
end subroutine findrms_prll



!//================================================================
!/* This function sets up the residual vector and finds the rms, given the ia and ja vectors, the A vector of [BLOCKSIZE][BLOCKSIZE] matrices, the solution vector x of dimension [dim][BLOCKSIZE], the rhs vector b of dimension [dim][BLOCKSIZE] the residual vector of dimension [dim*BLOCKSIZE], the dimension of the matrix dim, and the rms value rms */
subroutine findrms(ia,ja, A, x, b, res, rms,dim,neqn,nnz)
  !type(PComm_type),intent(inout)::pcomm
  integer :: dim,neqn,nnz
  real(dp) rms,grms,sum
  integer ,intent(in), dimension(dim+1) :: ia
  integer ,intent(in), dimension(nnz) :: ja
  real(dp),intent(in), dimension(neqn,neqn,nnz) :: A
  real(dp),intent(inout)  :: res(neqn,dim)
  real(dp),intent(inout) :: x(neqn,dim)
  real(dp),intent(in)  :: b(neqn,dim)
  integer :: j, k,ii, jj, jstart, jend,icol,ierr
  rms = 0.0
  do k=1,dim	
     do jj=1,neqn
        res(jj,k)=b(jj,k)
     end do
     !ceb needs a parallel matvec multiply
     jstart = ia(k)
     jend = ia(k+1)-1
     do j=jstart,jend
        icol=ja(j)
        do ii = 1,neqn
           do jj=1,neqn
              res(ii,k) = res(ii,k)-A(ii,jj,j)*x(jj,icol)
           end do
        end do
     end do
!write(6,*)"findrms: res(",1,k,")=",res(1,k)
     do j=1,neqn
        rms=rms+ res(j,k)*res(j,k)
     end do
  end do
  grms=rms
  !ierr=GetScalarSum_r(pcomm,rms,grms)
  rms = sqrt(grms)
  return
end subroutine findrms







!======================================================
!applies right preconditioner (preconditioning matrix should already exist in case of ILU)
subroutine precondition(nnodes,nsub,nnz,ia,ja,iau,A,x,b,&
                        nnz_pc,ia_pc,ja_pc,iau_pc,ALU,&
                        nsubiters,mypre,&
                        op)
  !type(PComm_type),intent(inout)::pcomm
  integer :: i,j,nnodes,nnz,nsub,nsubiters,mypre,nnz_pc
  real(dp), intent(in),   dimension (nsub,nsub,nnz)  :: A
  real(dp), intent(in),   pointer,dimension(:,:)     :: b
  integer,  intent(in),   dimension(nnodes+1)        :: ia
  integer,  intent(in),   dimension(nnz)             :: ja
  integer,  intent(in),   dimension(nnodes)          :: iau
  integer,  intent(in),   dimension(nnodes+1)        :: ia_pc
  integer,  intent(in),   dimension(nnz_pc)          :: ja_pc
  integer,  intent(in),   dimension(nnodes)          :: iau_pc
  real(dp), intent(in),   dimension(nsub,nsub,nnz_pc):: ALU

  real(dp),intent(inout), pointer, dimension(:,:)     :: x
  type(parms_Operator),intent(in),pointer ::op!ceb

  real(dp)::rms
  real(dp),dimension(nsub,nnodes)::res

!solves Ax=b


     !do i=1,nnodes
     !   do j=1,nsub
     !      write(6,*)"precond: initial b(",j,i,")=",b(j,i)
     !   end do
     !end do


! ceb Maybe just pass in A or ALU as the same variable A and then 
! ceb use A in all functions since none of these use both A and ALU
! ceb  probably do the same thing with ia,ja,iau,nnz

  ! y[k]=(M^-1)(v[k]/||res||)
  select case(mypre) 
  case(1)!diagonal pc
     call diag_invpc(iau,A,x,b,nnodes,nsub,nnz)
  case(2)!gauss seidel pc
     call SOLVE(nnodes,nsub,nnz,ia,ja,iau,A,x,b,nsubiters,0) ! dq = y[k] 
  case(3)!kaczmarz pc
     call BKACZMARZ(nnodes,nsub,nnz,ia,ja,iau,A,x,b,nsubiters,0)
  case(4:5)!ILU(K) preconditioned matrix ALU should already exist here 
     call BLK_LU_INV_SOLVE(nnodes, nsub, nnz, b, x, ALU, ia, ja, iau)
  case(6)!ILUT 
     call BLK_LU_INV_SOLVE(nnodes, nsub, nnz_pc, b, x, ALU, ia_pc, ja_pc, iau_pc)
  case(7)!arms 
     call parms_arms_sol_vcsr(op, b, x)
  case default  ! no preconditioning
     call QTOQ(nnodes,x,b,nsub) ! copy b into x 
  end select

     !do i=1,nnodes
     !   do j=1,nsub
     !      write(6,*)"precond: final x(",j,i,")=",x(j,i)
     !   end do
     !end do
!ceb
!call findrms(ia,ja,A,x,b,res,rms,nnodes,nsub,nnz)
!write(6,*)"precond rms=",rms
!ceb


  return
end subroutine precondition
!======================================================




!=================================== KACZMARZ ==================================
!
! Solve the matrix problem
! A =   Matrix
! B =   Scratch
! D =   norm of each row of A (nsub,nnodes)
! dq =  Solution
! rhs = Right-hand-side
!
! Solve using Kacsmarz relaxation
!
!=============================================================================
subroutine KACZMARZ(nnodes,nsub,nnz,ia,ja,iau,A,dq,rhs,nsubiters,iwriteres)
  integer, intent(in) :: nnodes,nsub,nnz,nsubiters,iwriteres
  real(dp),intent(in), dimension(nsub,nsub,nnz) ::  A
  real(dp),intent(in), dimension(nsub,nnodes) ::  rhs
  real(dp),intent(inout),dimension(nsub,nnodes) :: dq
  integer, intent(in), dimension(nnodes+1) :: ia
  integer, intent(in), dimension(nnodes) :: iau
  integer, intent(in), dimension(nnz) :: ja

  real(dp), dimension(nsub,nnodes) :: B
  real(dp), dimension(nsub,nnodes) :: D

  real(dp) ::resmonitor(25)
  integer :: i,l,j,k,n,idiag,jstart,jend,icol,istart,iend,m,ii
  real(dp) :: dotProduct,factor,rms

  ! Initialize dq to zero
  !call blank_2D_mat(dq,nsub,nnodes)

  ! Compute the squared norm for each row of A
  nodeLoop1: do i = 1,nnodes    ! Loop over nodes
     do j = 1,nsub
        D(j,i) = 0.
     end do
     !
     do j = ia(i),ia(i+1)-1
        do k = 1,nsub
           do l = 1,nsub
              D(k,i) = D(k,i) + A(k,l,j)*A(k,l,j)
           end do
        end do
     end do
  end do nodeLoop1

  ! Now iterate
  iterationLoop: do n = 1,nsubiters
     nodeLoop: do i = 1,nnodes
        !       call random_number(xx)
        !       i = 1 + xx*(nnodes - 1)
        !
        ! Loop over each row in block (for 4x4 loop over 4)
        ! We need to update all the dq's each time (with blocksize nsub, we will update all the dq's nsub times for each node
        !
        idiag = iau(i)
        kloop: do k = 1,nsub  ! There are nsub rows for each block

           dotProduct = 0.
           jstart = ia(i)
           jend = ia(i+1) - 1
           jloop: do j = jstart,jend
              icol = ja(j)
              do l = 1,nsub
                 dotProduct = dotProduct + A(k,l,j)*dq(l,icol)
              end do
           end do jloop  ! End of j loop over row of A

           ! Now we have the inner product so Update solution 
           ! Note that each time we compute the dot product, we update the guess for all elements of dq
           ! Update dq

           factor = 1.0
           !factor = factor*(-rhs(k,i) - dotProduct)/D(k,i) 
           factor = factor*(rhs(k,i) - dotProduct)/D(k,i) !ceb
           
           ! Now update dq
           jstart = ia(i)
           jend = ia(i+1) - 1
           do j = jstart,jend
              icol = ja(j)
              do l = 1,nsub
                 dq(l,icol) = dq(l,icol) + factor*A(k,l,j)
              end do
           end do ! loop over j
           !GS         end if ! End of if for doing Gauss-Seidel
        end do kloop ! End k loop over rows in block
        
     end do nodeLoop
     
     ! Now compute the residal
     
     ! Set the right-hand side to -rhs
     if(iwriteres.eq.1)then
        do i = 1,nnodes
           do j = 1,nsub
              !B(j,i) =  -rhs(j,i)
              B(j,i) =  rhs(j,i) !ceb
           end do
        end do

        ! Now, loop over the nodes and add the offdiagonal pieces
        rms = 0.0
        do  i = 1,nnodes
           istart = ia(i)       ! Row of matrix starts here
           iend   = ia(i+1) - 1 ! Row of matrix ends here
           idiag  = iau(i)      ! Diagonal
           m = i
           ! Residuals for each row
           do j = 1,nsub
              resmonitor(j) = 0.
           end do

           ! Loop over the columns in the matrix
           do  ii = istart,iend
              icol = ja(ii)
              do j=1,nsub
                 do k = 1,nsub
                    B(j,i) = B(j,i) - A(j,k,ii)*dq(k,icol)
                 end do
              end do
           end do
           do j = 1,nsub
              resmonitor(j) = resmonitor(j) + B(j,i)
           end do

           ! Finish computing the residual
           !
           !         do j = 1,nsub
           !           do k = 1,nsub
           !           resmonitor(j) = resmonitor(j) - A(j,k,idiag)*dq(k,i)
           !           end do
           !         end do
           
           do j = 1,nsub
              rms = rms + resmonitor(j)*resmonitor(j)
           end do
        end do
        rms = sqrt(rms/real(nsub*nnodes))
        write(6,*)"In kacsmarz n = ",n," rms = ",rms
     end if ! End of if(writeres)
     !
  end do iterationLoop
  
  return
end subroutine kaczmarz


!=================================== BKACZMARZ =================================
!
! Solve the matrix problem
! A =   Matrix
! B =   Scratch
! D =   A(transpose)*A for each node (4,4,nnodes)
! dq =  Solution
! rhs = Right-hand-side
!
! Solve using block Kacsmarz relaxation
!
!=============================================================================
!subroutine BKACZMARZ(nnodes,nsub,nnz,ia,ja,iau,A,B,D,dq,rhs,nsubiters,iwriteres)
subroutine BKACZMARZ(nnodes,nsub,nnz,ia,ja,iau,A,dq,rhs,nsubiters,iwriteres)
  integer, intent(in) :: nnodes,nsub,nnz,nsubiters,iwriteres
  real(dp),intent(in), dimension(nsub,nsub,nnz) ::  A
  real(dp),intent(in), dimension(nsub,nnodes) ::  rhs
  real(dp),intent(inout),dimension(nsub,nnodes) :: dq
  integer, intent(in), dimension(nnodes+1) :: ia
  integer, intent(in), dimension(nnodes) :: iau
  integer, intent(in), dimension(nnz) :: ja

  real(dp), dimension(nsub,nnodes) :: B
  real(dp), dimension(nsub,nsub,nnodes) :: D

  real(dp)  resmonitor(25),sum,factor,rms

  integer :: idiag
  integer :: row,column,segment
  integer :: i,j,k,ii,jj,kk,n,iinc,jstart,jend,L,istart,iend,icol,m
  !      integer :: icount(25)

  ! Initialize dq to zero
  !call blank_2D_mat(dq,nsub,nnodes)

  do i = 1,nnodes
     do j=1,nsub
        B(j,i) = 0.0 
        do k = 1,nsub
           D(j,k,i) = 0.
        end do
     end do
  end do

  !
  ! Compute A(transpose)*A for the row correcponsing to each node
  nodeLoop1: do i = 1,nnodes    ! Loop over nodes
     do j = ia(i),ia(i+1)-1
        do ii=1,nsub
           do jj=1,nsub
              do kk=1,nsub
                 D(ii,jj,i) = D(ii,jj,i) +  A(ii,kk,j)*A(jj,kk,j)          
              end do
           end do
        end do
     end do
  end do nodeLoop1

  ! Do the LU decomposition for D
  call DLU(D,nnodes,nsub) 
  
  ! Now iterate
  iterationLoop: do n = 1,nsubiters
     istart = 1
     iend = nnodes
     iinc = 1
     !if((n/2)*2.eq.n)then !if this expression is equal the make nodeloop backwards from nnodes to 1
     !   istart = nnodes
     !   iend = 1
     !   iinc = -1
     !end if
     nodeLoop: do i = istart,iend,iinc

        ! Initialize B to zero
        do k = 1,nsub
           B(k,i) = 0.
        end do

        ! Compute residual for row
        ! First compute row x dq
        jstart = ia(i)
        jend = ia(i+1) - 1
        jloop: do j = jstart,jend
           icol = ja(j)
           do k = 1,nsub
              do L = 1,nsub
                 B(k,i) = B(k,i) + A(k,L,j)*dq(L,icol)
              end do
           end do
        end do jloop  ! End of j loop over row of A
        !
        ! Now add piece to residual from right-hand-side
        do k = 1,nsub
           !B(k,i) =  -rhs(k,i) - B(k,i) 
           B(k,i) =  rhs(k,i) - B(k,i) !ceb
        end do
        !
        ! Solve for y(1-nsub) by solving [D]*y = B
        

        ! solve for y, writing over B
        ! Forward
        m = i
        B(1,m) = B(1,m)/D(1,1,i)
        do j = 2,nsub
           sum = B(j,m)
           do k = 1,j-1
              sum = sum - D(j,k,i)*B(k,m)
              !sum = sum + D(j,k,i)*B(k,m) !ceb
           end do
           B(j,m) = sum/D(j,j,i)
        end do

        ! Backward (note that B(nsub,m) is already obtained so we start at the next row up)
        do j = nsub-1,1,-1
           sum = B(j,m)
           do k = nsub,j+1,-1
              B(j,m) = B(j,m) - D(j,k,i)*B(k,m)
           end do
        end do
        
        factor = 1.0
        
        ! Now update dq
        jstart = ia(i)
        jend = ia(i+1) - 1
        do j = jstart,jend
           icol = ja(j)
           do k = 1,nsub
              do L = 1,nsub
                 dq(k,icol) = dq(k,icol) + factor*A(L,k,j)*B(L,i)
              end do
           end do
        end do ! loop over j
        
     end do nodeLoop



     ! Now compute the residal

     ! Set the right-hand side to -rhs
     monitorRes: if(iwriteres.eq.1)then
        do i = 1,nnodes
           do j = 1,nsub
              !B(j,i) =  -rhs(j,i)
              B(j,i) =  rhs(j,i) !ceb
           end do
        end do

        ! Now, loop over the nodes and add the offdiagonal pieces
        rms = 0.0
        do  i = 1,nnodes
           istart = ia(i)       ! Row of matrix starts here
           iend   = ia(i+1) - 1 ! Row of matrix ends here
           idiag  = iau(i)      ! Diagonal
           m = i
           ! Residuals for each row
           do j = 1,nsub
              resmonitor(j) = 0.
           end do
 
           ! Loop over the columns in the matrix
           do  ii = istart,iend
              icol = ja(ii)
              do j=1,nsub
                 do k = 1,nsub
                    B(j,i) = B(j,i) - A(j,k,ii)*dq(k,icol)
                 end do
              end do
           end do
           do j = 1,nsub
              resmonitor(j) = resmonitor(j) + B(j,i)
           end do

           ! Finish computing the residual         
           do j = 1,nsub
              rms = rms + resmonitor(j)*resmonitor(j)
           end do
        end do
        rms = sqrt(rms/real(nsub*nnodes))
        write(6,*)"In kacsmarz n = ",n," rms = ",rms
     end if monitorRes! End of if(writeres)
     !
  end do iterationLoop
  
  return
end subroutine bkaczmarz


!=================================== KLEAST ==================================
!
! Solve the matrix problem using the least squares method
! A =   Matrix
! B =   Scratch
! D =   norm of each row of A (nsub,nnodes)
! dq =  Solution
! rhs = Right-hand-side
!
!=============================================================================
subroutine KLEAST(nnodes,nnz,ia,ja,iau,ka,A,B,D,dq,rhs,nsub,nsubiters,iwriteres)
  integer, intent(in) :: nnodes,nnz,nsub,nsubiters,iwriteres
  real(dp)  A(nsub,nsub,nnz),D(nsub,nsub),B(nsub,nnodes)
  real(dp)  rhs(nsub,nnodes)
  real(dp)  dq(nsub,nnodes)
  real(dp)  resmonitor(25),rms
  integer ia(nnodes+1),iau(nnodes),ja(nnz),ka(nnz)
  integer :: idiag
  integer :: i,j,k,ii,jj,kk,kstart,kend,length,irowIndex,irow,L,m,n,kcount,icol,istart,iend
  !     integer :: icount(25)
  !
  ! Initialize dq to zero
  !
  do  i = 1,nnodes
     do  j=1,nsub
        !dq(j,i) = 0.0
        D(j,i) = 0.
     end do
  end do
  !
  ! Now iterate
  !
  iterationLoop: do n = 1,nsubiters
     nodeLoop: do i = 1,nnodes
        !
        ! First form D, which is the transpose of the columns time the columns
        !
        do j = 1,nsub
           do k = 1,nsub
              D(j,k) = 0.
           end do
        end do
        
        kstart = ia(i) ! This tells us where in the ja and ka arrays to start
        kend   = ia(i+1) - 1
        length = kend - kstart ! How many elements are in this column
        kloop: do k = 1,length
           irowIndex = ka(kstart + k - 1) ! This is the index into the A array that is on the row we are interested in
           irow = ja(kstart + k - 1)
           do L = 1,nsub
              B(L,k) = -rhs(L,icol) ! Initialize the row sum
              do m = 1,nsub
                 D(l,m) = D(l,m) + A(m,l,irow)*A(l,m,irow)
              end do
           end do
        end do kloop
        !
        ! Form the row sums. Remember that ja contains the columns and the rows because the matrix has symmetric form
        !
        kcount = 0
        kloop2: do k = kstart,kend ! This tells us the rows
           kcount = kcount + 1
           irow = ja(k)
           Lloop: do L = ia(irow),ia(irow+1) - 1 ! Loop over all the elements in row(irow)
              icol = ja(L)
              if(icol.ne.i)then
                 do jj=1,nsub
                    do kk=1,nsub
                       B(jj,kcount) = B(jj,kcount) - A(jj,kk,L)*dq(kk,icol)
                    end do
                 end do
              end if
           end do Lloop
        end do kloop2
        !
        ! Now multiply the B by the transpose of the elements in column i
        !         
        !         do l = 
        !         irow = ia(
        !       
        !       end do kloop2
        !         icount = icount + 1
        !         do l = 1,nsub
        !           R(l,icount) = 0.
        !         end do
        
        !         jstart = 
        !
        ! Now solve for dq
        !
     end do nodeLoop
     !
     ! Now compute the residal
     !
     ! Set the right-hand side to -rhs
     !
     if(iwriteres.eq.1)then
        do i = 1,nnodes
           do j = 1,nsub
              B(j,i) =  -rhs(j,i)
           end do
        end do
        !
        ! Now, loop over the nodes and add the offdiagonal pieces
        !
        rms = 0.0
        do  i = 1,nnodes
           istart = ia(i)       ! Row of matrix starts here
           iend   = ia(i+1) - 1 ! Row of matrix ends here
           idiag  = iau(i)      ! Diagonal
           m = i
           ! Residuals for each row
           do j = 1,nsub
              resmonitor(j) = 0.
           end do
           !
           ! Loop over the columns in the matrix
           !
           do  ii = istart,iend
              icol = ja(ii)
              do j=1,nsub
                 do k = 1,nsub
                    B(j,i) = B(j,i) - A(j,k,ii)*dq(k,icol)
                 end do
              end do
           end do
           do j = 1,nsub
              resmonitor(j) = resmonitor(j) + B(j,i)
           end do
           !
           ! Finish computing the residual
           !
           !         do j = 1,nsub
           !           do k = 1,nsub
           !           resmonitor(j) = resmonitor(j) - A(j,k,idiag)*dq(k,i)
           !           end do
           !         end do
           
           do j = 1,nsub
              rms = rms + resmonitor(j)*resmonitor(j)
           end do
        end do
        rms = sqrt(rms/real(nsub*nnodes))
        write(6,*)"In kacsmarz n = ",n," rms = ",rms
     end if ! End of if(writeres)
     !
  end do iterationLoop
  
  return
end subroutine KLEAST






!=================================== GMRES ===================================
! GMRES with right preconditioning
! Solve Ax + b = zero. 
! To do this solve A*M(inv)*u + b = 0 (or Abar*u + b = 0)
! and then get x by solving Mx = u
!=============================================================================
subroutine GMRES(nnodes,k,icycle,tol,ierror,mypre,nsubiters,nnz,nsub,&
                 ia,ja,iau,A,B,q,fill_level,droptol,&
                 op)
  ! 
  !     PURPOSE 
  !       GMRES COMPUTES THE SOLUTION OF REAL LINEAR
  !       SYSTEM OF EQUATIONS B+AX =0, USING GMRES ALGORITHM
  ! 
  !   INPUT PARAMETERS
  !       NEQ   NUMBER OF EQUATIONS (nsub*nnodes)
  !       X     INITIAL APPROXIMATE SOLUTION
  !       K     NUMBER OF SEARCH DIRECTIONS 
  !       ICYCLE   NUMBER OF GMRES CYCLES (ROUGHLY (K+1)*ICYCLE CALLS TO FCN) 
  !       TOL   CALCULATION WILL TERMINATE IF INITIAL RESIDUAL
  !             ERROR IS REDUCED BY FACTOR TOL. 
  !       NQ    LOGICAL TAPE NUMBER OF FILE (ON SSD) GMRES IS TO USE. 
  !             FILE SHOULD BE OPENED AND CLOSED BY CALLING ROUTINE 
  !             LREC=NEQ,  NUMBER OF RECORDS REQUIRED IS (K+1). 
  !       F     WORK MATRIX USED IN PLACE OF SSD. MATRIX SHOULD BE
  !             DIMENSIONED TO BE (NEQ,K+1)
  !       A     Matrix
  !       ALU   Copy of Matrix for ILU
  !       ICYCLE Number of Restarts
  !       q     Solution vector (nsub,nnodes)
  ! 
  !    INTERNALLY USED ONLY 
  !      H      HESSENBERG MATRIX 
  !      G      RIGHT HAND SIDE OF K DIMENSIONAL LEAST
  !             SQUARES PROBLEM.
  !      C      DIRECTION COSINES 
  !      S      DIRECTION SINES 
  ! 
  !   OUTPUT 
  !      X      FINAL APPROXIMATE SOLUTION
  !      AP     FINAL VALUE OF B+AX 
  !      IERROR  ERROR FLAG 
  ! 
  !            IERROR=-1    MATRIX IS SINGULAR
  !            IERROR=0     RESIDUAL REDUCED BY FACTOR TOL
  !            IERROR=1     RESIDUAL NOT REDUCED BY FACTOR TOL
  !                         IN ICYCLES OF K STEPS EACH. 
  !type(PComm_type),intent(inout)::pcomm
  integer, intent(in) :: k,icycle,nnz,nnodes,nsub,mypre,nsubiters,fill_level
  integer, intent(inout) :: ierror
  integer, intent(in), dimension(nnodes+1) :: ia
  integer, intent(in), dimension(nnz) :: ja
  integer, intent(in), dimension(nnodes) :: iau
  real(dp), intent(in) :: tol
  real(dp), intent(in) :: droptol
  real(dp), intent(inout), dimension(nsub,nnodes) ::   q
  real(dp), intent(inout), dimension(nsub,nnodes) ::   B! rhs
  real(dp), intent(inout), dimension(nsub,nsub,nnz) :: A

  real(dp),target, dimension(:,:), allocatable ::     H!(maxsearch,maxsearch)
  real(dp), dimension(:), allocatable ::       g!(maxsearch) ! Upper Hessenberg matrix
  real(dp), dimension(:), allocatable ::       c!(maxsearch)
  real(dp), dimension(:), allocatable ::       s!(maxsearch)           ! Givens rotations
  !real(dp), dimension(:), allocatable ::       x! Working array inside gmres routine
  real(dp), dimension(:), allocatable ::       x! Working array inside gmres routine
  real(dp), dimension(:), allocatable ::       ap!matrix vector products for gmres
  real(dp),target, dimension(:,:), allocatable ::     f!
  real(dp), dimension(:,:,:), allocatable  ::  ALU! preconditioning matrix
  !real(dp), dimension(:,:), allocatable  ::    phi
  real(dp), dimension(:,:), pointer  ::    phi
  !real(dp), dimension(:,:), allocatable  ::    dq
  real(dp), dimension(:,:), pointer  ::    dq
  real(dp), dimension(:,:), allocatable  ::    BB! copy of B
  real(dp), dimension(:,:,:), allocatable  ::  D! diagonal block matrix

  integer :: nnz_pc
  integer, allocatable, dimension(:) :: ia_pc
  integer, allocatable, dimension(:) :: ja_pc
  integer, allocatable, dimension(:) :: iau_pc
  real(dp), dimension(:,:), allocatable  ::    resid
  real(dp) :: res,resnrm,bnorm,apnorm,gdot,beta,delta
  integer :: idolu,ind,ntot,L,J,I,IR,istat,JIT,LIT,kcol,jold,istart,iend,inode,ii
  integer, dimension(nnodes) :: iw
  real(dp) :: rms,tempc,RELRES
  real(dp),pointer::shift_real_ptr
  integer,pointer::shift_int_ptr

  type(parms_Operator),intent(inout),pointer ::op!ceb

  integer:: nx,np
  nx=1
  np=1

  allocate(phi(nsub,nnodes),STAT = istat)
  allocate(dq(nsub,nnodes),STAT = istat)
  allocate(D(nsub,nsub,nnodes),STAT = istat)
  allocate(H(k,k),STAT = istat)
  allocate(G(k),STAT = istat)
  allocate(C(k),STAT = istat)
  allocate(S(k),STAT = istat)
  allocate(resid(nsub,nnodes),STAT = istat)

  allocate(X(nsub*nnodes),STAT = istat)
  allocate(AP(nsub*nnodes),STAT = istat)
  allocate(F(nsub*nnodes,k+1),STAT = istat)

  ! Initialize u and store in dq
  call blank_2D_mat(dq,nsub,nnodes)

  ! First copy dq(nsub,nnodes) into x(nsub*nnodes)
  !this value is actually q not dq so initial u should actually be set to res
  !call DQX(nnodes,q,X,nsub)  ! Set initial guess X=u0 <- q
  call DQX(nnodes,dq,X,nsub)  ! Set initial guess X=u0 <- q

  IERROR=0
  ntot=nnodes*nsub

  CALL SCOPY(ntot,X,1,F(1:,1),1)  ! Copy u0 to F(.,1) ! at least initialize F

  ! Evaluate (b-A*u) and store in AP
  !generate ALU factorization
  call ILUPRE(nnodes,nnz,nsub,1,A,ALU,ia,ja,iau,phi,mypre,nnz_pc,ia_pc,ja_pc,iau_pc,fill_level,droptol)
  !call cond_inf_LU(ia,ja,iau,A,ALU,B,nnodes,nsub,nnz,mypre) !ceb test

  ! Compute r0 = b - Abar*u0
  !call findrms(ia,ja,A,q,B,phi,bnorm,nnodes,nsub,nnz)
  call findrms(ia,ja,A,dq,B,phi,bnorm,nnodes,nsub,nnz)
  !call findrms_prll(pcomm,ia,ja,A,q,B,phi,bnorm,nnodes,nsub,nnz)
  call DQX(nnodes,phi,AP,nsub) !move 2d res (phi) -> 1d ap

  ! bnorm = ||r0||
  IF (BNORM .le. tol) GO TO 2000 
  RESNRM=BNORM
  write(6,*)' ICYCLE=',icycle,' BNORM=',bnorm
  !---------------------------------
  
  ! Now solve Abar*u + r0 = 0
  DO  L=1,ICYCLE  ! Loop over cycles (icycle 1 means no restart)
     !{
     RES=RESNRM
     delta=RESNRM 

     g(1) = res !;//g=rms*e[1]
     delta=res

     if (g(1).eq.0.0) then
        write(6,*)"Exiting GMRES solve: initial solution is correct"
     else

        DO  J=1,K  ! Loop over dimension of Krylov subspace
           !{
           G(J)=RES      ! This is the right-hand side of the least squares problem
           ! It will start with beta*e1 = ||r0|| but new entries will
           ! get added and it will continually change as Givens rotations are applied
           ! Normalize v
           CALL SSCAL(ntot,1./delta,AP,1)! here AP = v[k] = v[k]/||res||
           ! Copy ap=v into F(*,NP+J) Note that NP = 1 so we will copy into 2nd column and higher
           CALL SCOPY(ntot,AP,1,F(1:ntot,NP+J),1)     !F(*,NP+J)=     v[k]/||res||

!do ii=1,ntot
!write(6,*)"1: F(",ii,")=",F(ii,NP+J)
!end do

           ! Copy v into X
           CALL SCOPY(ntot,AP,1,X,1) !copies ap=v to x        x=  v[k]/||res||

           ! begin Arnoldi process----------------------------
           ! u[k] = A*y[k] = A*(M^-1)v[k],   u=ap  AP=(A^-1)v= (M^-1)v
           ! y[k]=(M^-1)(v[k]/||res||),   dq=y[k]
           !call precondition(nnodes,nsub,nnz,ia,ja,iau,A,dq,X,nnz_pc,ia_pc,ja_pc,iau_pc,ALU,nsubiters,mypre,op)
           call XDQ(nnodes,phi,X,nsub)!ceb copy 1d x to 2d phi
           call precondition(nnodes,nsub,nnz,ia,ja,iau,A,dq,phi,nnz_pc,ia_pc,ja_pc,iau_pc,ALU,nsubiters,mypre,op)!ceb

           call matvect(ia,ja,A,dq,phi,nnodes,nsub,nnz,1)! u[k] = A*y[k],   phi = u[k] 
           !call matvect_prll(pcomm,ia,ja,A,dq,phi,nnodes,nsub,nnz,1)! u[k] = A*y[k],   phi = u[k] 

           call DQX(nnodes,phi,AP,nsub)!copy 2d phi to 1d x ! AP = u[k]
           
           ! Orthoganolize against previous vectors using modified Gram-Schmidt
           DO  I=1,J
              CALL SCOPY(ntot,F(1:,NP+I),1,X,1) ! copies f to x,    x=v[j]  
!do ii=1,ntot
!write(6,*)"1: x(",ii,")=",x(ii)
!end do
              H(I,J)=SDOT(ntot,X,1,AP,1)         ! h[j,k] = v[j]^T dot u[k]

              CALL SAXPY(ntot,-H(I,J),X,1,AP,1)  ! AP = u[k] = u[k]-h[j,k]*v[j]

           end DO
           ! end Arnoldi proceess------------------------------
           
           ! Apply Givens rotations to H for all the previous cos and sin
           DO  I=1,J-1
              delta=H(I,J)
              tempc= C(I)*delta + S(I)*H(I+1,J)     ! save before altering h[j+1,k]
              H(I+1,J)= -S(I)*delta + C(I)*H(I+1,J) !h[j+1,k] = -s[j]*delta + c[j]*h[j+1,k]
              H(I,J)=tempc                           !h[j,k] = c[j]*delta + s[j]*h[j+1,k]
           end DO
           
           ! do i=1,j-1
           !             beta = h(i,j)!// beta = h[i,j]
           !             h(i,j) = c(i)*beta + s(i)*h(i+1,j)!// h[i,j] = c[i]*delta + s[i]*h[i,j+1]
           !             h(i+1,j) = -s(i)*beta + c(i)*h(i+1,j)!// h[j+1,k] = -s[j]*delta + c[j]*h[j+1,k]
           !          end do
           
           ! Get current cos and sin and apply Givens to this row including RHS
           delta=SNRM2(ntot,AP,1) ! apnorm=||u[k]||
           tempc=delta
           
           CALL SROTG(H(J,J),TEMPC,C(J),S(J))
           !// Apply current Givens rotations to new column h[j,k] ==================================== //all local
           RES=-S(J)*G(J) 
           G(J)=C(J)*G(J) 
           
           RESNRM=ABS(RES)
           RELRES=RESNRM/BNORM 
           write(6,*)'krylov iter=',j,' resnrm=',resnrm !," tol=",tol," k=",k
           !write(6,*)'krylov iter=',j,' relres=',relres !," tol=",tol," k=",k
           JIT = J

           !IF (RELRES .LE. TOL .OR. J .EQ. K) GO TO 550
           IF (RESNRM .LE. TOL .OR. J .EQ. K) GO TO 550
           !}
        end DO ! End of loop over Krylov subspace
      
  
550     CONTINUE
        IF (ABS(H(J,J)) .EQ. 0.) THEN 
write(6,*) "Error: H diag==0 exiting"
           IERROR=-1 
           GO TO 2000
        END IF
        !========================================================





write(6,*) "exiting krylov loop "

        ! Solve Hy=g by backsubstitution
        ! After this, G holds y (recall that x = x0 + Z and Z = V*y)
        G(J)=G(J)/H(J,J)
!write(6,*) "g(",j,")=",g(j)
        DO  I=J-1,1,-1
           !G(I)=(G(I)- SDOT(J-I, G(I+1),1, H(I,I+1),K) )/H(I,I)
           do IR=I+1,J
              g(I) = g(I) - g(IR)*h(I,IR)!dot product
           end do
           g(I) = g(I)/H(I,I)
!write(6,*) "g(",i,")=",g(i)
        end DO
        


!do ii=1,ntot
!write(6,*)"2: x(",ii,")=",x(ii)
!end do

        ! We now know y so compute z=Vy and add to x

        ! u <- u + Vy (this gets put back into F(,.1)

        ! Recall that the rest of the V's are stored in F(.,2)
        CALL SCOPY(ntot,F(1:,NX),1,X,1)  ! Copy u0 into X (remember u0 is stored in F(.,1)

!do ii=1,ntot
!write(6,*)"3: x(",ii,")=",x(ii)
!end do


        !ceb we want z0 to be zero here?
        DO  I=1,J
           CALL SCOPY(ntot,F(1:,NP+I),1,AP,1)!ap=v[]
           CALL SAXPY(ntot,G(I),AP,1,X,1)!x=x+g*ap  u=u+g*v
        end DO

        CALL SCOPY(ntot,X,1,F(1:,NX),1) ! Copy u into F(,.1) !ceb
        !call XDQ(nnodes,dq,X,nsub) ! Copy X(.) into dq(.,.) (this puts u into dq)
        call XDQ(nnodes,phi,X,nsub)!ceb copy 1d x to 2d phi
 
        ! At this point, X, dq and F(,.1)contain u
        ! We need to compute x = (M^-1)*u
        ! Recall that we solved A*(M^-1)*u + r0 = 0 for u
        ! We now know u so we need to compute x from Mx = u
        ! We do this by applying the preconditoner
        
        !call QTOQ(nnodes,phi,dq,nsub) ! copy dq into phi 
        ! After qtoq, x, dq, and phi all contain u
        ! (this puts u into phi to use as the right-hand side for the preconditioner)
        ! x[k]=(M^-1)u, 
        call precondition(nnodes,nsub,nnz,ia,ja,iau,A,dq,phi,nnz_pc,ia_pc,ja_pc,iau_pc,ALU,nsubiters,mypre,op)
        
        ! After ilupre, dq contains the actual x whereas the phi and X arrays contain the current value of u. 
        ! Also note that F(.,1) also contains the current value of u
        do i=1,nnodes
           do j=1,nsub
write(6,*) "dq(",j,i,")",dq(j,i)
              q(j,i)=q(j,i)+dq(j,i)
           end do
        end do
        ! Now dq contains x so we can compute Ax + r0 to see if it converged
        ! To do this don't apply preconditioner. Result comes back in AP
        !call findrms_prll(pcomm,ia,ja, A, q, B, phi, resnrm,nnodes,nsub,nnz)
        !call findrms(ia,ja, A, q, B, phi, resnrm,nnodes,nsub,nnz)
        call findrms(ia,ja, A, dq, B, phi, resnrm,nnodes,nsub,nnz)
        write(6,*)"Residual at end of icycle",L,RESNRM
        RELRES=RESNRM/BNORM 

!write(6,*)" gmres chkpt 7 relres=",relres

        LIT = L
        !IF (RELRES .LE. TOL) GO TO 2000 
        IF (RESNRM .LE. TOL) GO TO 2000 
     end if
!write(6,*)" gmres chkpt 8"
  end DO! End of restart loop
  IERROR=1

2000 CONTINUE ! If we reached tolerance kick out to here
  !
  ! At this point, F(,.1) has u in it. It does not have x in it.
  ! It doesn't do any harm to copy into x but we don't need it
  ! so comment it out
  !     CALL SCOPY(ntot,F(1,NX),1,X,1) ! Copy F to X
!write(6,*)" gmres chkpt 9"

  ! Evaluate Abar*u + b (A*x + b) and store in AP
  ! Compute r0 = b + Abar*u0
  ! Result is in AP

  !write(6,105)lit,jit,resnrm
  !write(6,105)lit,jit,relres
105 format(1h 'gmres iters=',i6,' srch directions=',i6,' rms=',e14.7)
  !
  ! Copy x back into dq
  !     call xdq(nnodes,dq,X,nsub)

!write(6,*)" gmres chkpt 10"
  if(ALLOCATED(ALU))deallocate(ALU)
!write(6,*)" gmres chkpt 10.5"
  deallocate(D)
!write(6,*)" gmres chkpt 11"
  deallocate(phi)
!write(6,*)" gmres chkpt 13"
  deallocate(dq)
  deallocate(f)
  deallocate(x)
  deallocate(ap)
!write(6,*)" gmres chkpt 14"
  deallocate(H)
  deallocate(G)
  deallocate(C)
  deallocate(S)
write(6,*)" gmres chkpt final"

  return
end subroutine GMRES








!========================================================
!subroutine ILUK_factorization(p,nnodes,nsub,nnz,ia,iau,ja,A, ianew,iaunew,janew,list,oldtonew,jnewtojold,bwidth, nnzFill,iaFill,iauFill,jaFill,Afill)
subroutine ILUK_factorization(p,nnodes,nsub,nnz,A, ianew,iaunew,janew,jnewtojold,bwidth, nnzFill,iaFill,iauFill,jaFill,Afill)
  !parameters
  integer,intent(in) :: p,nnodes,nsub,nnz,bwidth
  integer,intent(inout) :: nnzfill
  !integer,intent(in) :: ianew(nnodes+1),iaunew(nnodes),janew(nnz),jnewtojold(nnz)
  !real(dp),intent(in) :: A(nsub,nsub,nnz)  !integer,intent(inout) :: iaFill(nnodes+1),iauFill(nnodes)
  !integer,intent(inout), dimension(:), allocatable :: jaFill
  !real(dp),intent(inout), dimension(:,:,:), allocatable ::  Afill
  !integer,intent(in) :: ianew(nnodes+1),iaunew(nnodes),janew(nnz),jnewtojold(nnz)

  !integer,intent(in) :: ianew(nnodes+1),iaunew(nnodes),janew(nnz),jnewtojold(nnz)
  real(dp),intent(in),pointer,dimension(:,:,:) :: A

  integer,intent(inout),pointer,dimension(:) :: iaFill,iauFill
  integer,intent(inout),pointer,dimension(:) :: jaFill
  real(dp),intent(inout),pointer,dimension(:,:,:) ::  Afill
  integer,intent(in),pointer,dimension(:) :: ianew,iaunew,janew,jnewtojold

  !local vars
  integer :: i,j,k,jj,kk,kstart,kend,jstart,jend,jcol,kindex,kcolumnStart,kcolumnEnd,jold,kcol
  integer :: jloopStart,jloopEnd,jindex,kstartk,kendk,jindexk
  integer :: ifill,ioldnode,nonZeroOnRow,icount
  integer, dimension(:), allocatable :: jscratch 
  real(dp) :: compareTop,compareLeft
  real(dp), dimension(:,:), allocatable ::  bigB


  !do i = 1,nnodes
  !  jstart = ianew(i)
  !  jend   = ianew(i+1) - 1
  !  do j = jstart,jend
  !    jcol = janew(j)
  !    write(46,'(2(i5,1x))')i,jcol
  !  end do
  !end do

! Now lets compute the fill pattern

! Now compute the same thing as in bigA (call it bigB) but use banded structure

      write(6,*)"Allocating bigB"
      allocate(bigB(nnodes,2*bwidth+1))
      write(6,*)"Finished allocating bigB"
      do i = 1,nnodes
        do k = 1,2*bwidth+1
          bigB(i,k) = 1.e10
        end do
      end do
      write(6,*)"Finished initializing bigB to infinity"

!     do i = 1,nnodes
!       kstart = max(1,i-bwidth)      ! First non-zero column in this row of bigB
!       kend   = min(nnodes,i+bwidth) ! Last non-zero column in this row of bigB
!       do k = kstart,kend
!         kindex = 1 + k - kstart     ! Index into bigB that is for column k
!         bigB(i,kindex) = 1.e10
!       end do
!     end do
!
! Put zeros in the locations where there are elements present in the original matrix
      do i = 1,nnodes
        kstart = max(1,i-bwidth)      ! First non-zero column in this row of bigB
        kend   = min(nnodes,i+bwidth) ! Last non-zero column in this row of bigB)
!
! Starting and ending columns for the matrix in compressed row storage
        jstart = ianew(i)
        jend   = ianew(i+1) - 1
        do j = jstart,jend
          jcol = janew(j)
          kindex = 1 + jcol - kstart  ! Index into bigB that is for column k
          bigB(i,kindex) = 0.
        end do
      end do
      write(6,*)"Finished setting non-zero elements to zero"

      write(6,*)" Now starting bigB "
      do i = 2,nnodes
!        if((i/100)*100.eq.i)write(6,*)"Working row ",i
        kstart = max(1,i-bwidth)       ! First non-zero column in this row of bigB
        kend   = min(nnodes,i+bwidth)  ! Last non-zero column in this row of bigB)

        kcolumnStart = kstart
        kcolumnEnd   = i-1
        do k = kcolumnStart,kcolumnEnd
          kindex = 1 + k - kstart
  
          if(bigB(i,kindex).le.p)then
            jloopStart = k+1   ! Loop from first element in row to the last element in the row
            jloopEnd   = kend
            do j = jloopStart,jloopEnd
              jindex = 1 + j - kstart ! Index into bigB array for element (i,j)
              compareLeft = bigB(i,kindex)  ! level for element on the same row(i) but different column(k) level(i,k)
              compareTop = 1.e10            ! level for element on the same column(j) but different row(k) level(k,j)
              kstartk = max(1,k-bwidth)    ! Starting column for row k
              kendk = min(nnodes,k+bwidth) ! Ending column for row k
              if(j.le.kendk)then ! If column j is less than or equal to the number of columns on row k
                jindexk = 1 + j - kstartk ! Index for column j into bigB on row k
                 compareTop = bigB(k,jindexk)
              end if
              bigB(i,jindex) = min(bigB(i,jindex), compareLeft + compareTop + 1.)
            end do
          end if
        end do
      end do
!
! Now output for plotting
      write(6,*)"Finished computing bigB"

      do i = 1,nnodes
        kstart = max(1,i-bwidth)       ! First non-zero column in this row of bigB
        kend   = min(nnodes,i+bwidth)  ! Last non-zero column in this row of bigB)
        do k = kstart,kend
          kindex = 1 + k - kstart      ! Index into bigB that is for column k
          if(bigB(i,kindex).le.p)write(48,'(2(i5,1x))')i,k
        end do
      end do
!     stop

! If we want to allow fill, get new arrays for compressed row storage

! Fill the new matrix using bigB
!
!---------------------------- fill using bigB ----------------------
!
! If we want to allow fill, get new arrays for compressed row storage

  !9983 continue 
!      goto 9984 
      ifill = 1
        allocate(jscratch(nnodes))
        !allocate(iaFill(nnodes+1))       ! The new ia
        !allocate(iauFill(nnodes))        ! The new iau

! Figure out how many non-zeros there are
        nnzFill = 0
        iaFill(1) = 1
        do i = 1,nnodes
          nonZeroOnrow = 0
          kstart = max(1,i-bwidth)       ! First non-zero column in this row of bigB
          kend   = min(nnodes,i+bwidth)  ! Last non-zero column in this row of bigB)
          do k = kstart,kend
            kindex = 1 + k - kstart
            if(bigB(i,kindex).le.p)then
              nnzFill = nnzFill + 1
              nonZeroOnRow = nonZeroOnRow + 1
            end if
          iaFill(i+1) = iaFill(i) + nonZeroOnRow
          end do
        end do
        allocate(jaFill(nnzFill)) ! Allocate enough space for adding a fill for every node
        write(6,*)"nnzFill = ",nnzFill

! Now fill jaFill
        icount = 0
        do i = 1,nnodes
          kstart = max(1,i-bwidth)       ! First non-zero column in this row of bigB
          kend   = min(nnodes,i+bwidth)  ! Last non-zero column in this row of bigB)
          do k = kstart,kend
            kindex = 1 + k - kstart
            if(bigB(i,kindex).le.p)then
              icount = icount + 1
              jaFill(icount) = k
              if(k.eq.i)iauFill(i) = icount
            end if 
          end do
        end do

! Now copy the A matrix into the new matrix (Afill)
       allocate(Afill(nsub,nsub,nnzFill))

       do i = 1,nnzFill
         do j = 1,nsub
           do k = 1,nsub
             Afill(j,k,i) = 0.
           end do
         end do
       end do

!
! Loop over rows and fill the new matrix values with the original ones
!
       do i = 1,nnodes
         !ioldnode = list(i)           ! This is the corresponding row in the original matrix (even before CM) (list is the same as newtoold)
         !jstart = ia(ioldnode)        ! Start of row in original matrix
         !jend   = ia(ioldnode+1) - 1  ! End of row in original matrix
         jstart = ianew(i)        ! Start of row in original matrix
         jend   = ianew(i+1) - 1  ! End of row in original matrix

         kstart = iaFill(i)           ! Start of row in new matrix
         kend   = iaFill(i+1) - 1     ! End of row in new matrix
!
! Now loop over the columns in the old row
! For each element of A, find where in the new row it belongs
!
         do j = jstart,jend ! Original row
           !jcol = ja(j)  ! Original column
           jcol = janew(j)  ! Original column

           !jcolnew = oldtonew(jcol) ! New location for old column

!ceb need to convert j to jold for access to A in pre cuthill-mckee ordering
           jold=jnewtojold(j) !ceb

           do k = kstart,kend ! New row
             kcol = jaFill(k) ! New column
! If new column matches the old column, put this value of A into Afill
             !if(kcol.eq.jcolnew)then
             if(kcol.eq.jcol)then
                do jj=1,nsub
                   do kk=1,nsub
                      Afill(jj,kk,k) = A(jj,kk,jold)
                   end do
                end do

               goto 5464 ! If we've found it, kick out
             end if
           end do ! k loop
 5464 continue

         end do ! jloop
       end do ! loop over rows

  9984 continue

! Write out the matrix
!

!      write(99,*)"Dimension"
!      write(99,*)nnodes
!      write(299)nnodes
!      write(99,*)" nnz"
!      write(99,*)nnzFill
!      write(299)nnzFill
!      write(99,*)"IA dimension=",nnodes+1
!      do 5000 i = 1,nnodes+1
!        write(99,*)iaFill(i)
!        write(299)iaFill(i)
! 5000 continue
!      write(99,*)"IAU dimension=",nnodes
!      do 5010 i = 1,nnodes
!      write(99,*)iauFill(i)
!      write(299)iauFill(i)
! 5010 continue
!      write(99,*)"JA nnz=",nnzFill
!      do 5080 i = 1,nnzFill
!         write(99,*)jaFill(i)
!         write(299)jaFill(i)
! 5080 continue
!      write(99,*)"A dimension=",3*nnzFill," (3*nnz)"
!      do 6000 i = 1,nnodes
!        jstart = iaFill(i)
!        jend   = iaFill(i+1) - 1
!        do 6010 j = jstart,jend
!          write(99,'(6(e14.7,1x))')Afill(1,1,j),Afill(1,2,j),Afill(1,3,j)
!          write(99,'(6(e14.7,1x))')Afill(2,1,j),Afill(2,2,j),Afill(2,3,j)
!          write(99,'(6(e14.7,1x))')Afill(3,1,j),Afill(3,2,j),Afill(3,3,j)
!          write(299)Afill(1,1,j),Afill(1,2,j),Afill(1,3,j)
!          write(299)Afill(2,1,j),Afill(2,2,j),Afill(2,3,j)
!          write(299)Afill(3,1,j),Afill(3,2,j),Afill(3,3,j)
! 6010   continue
! 6000 continue


!
! Write out the right-hand side
!
!      write(99,*)"Right hand side dimension=",nnodes
!      do 5060 i = 1,nnodes
!        index = list(i)
!        write(99,'(6(e14.7,1x))')(qnode(j,index),j=1,3)
!        write(299)(qnode(j,index),j=1,3)
! 5060 continue
!      write(99,*)"oldtonew =",nnodes
!      do i = 1,nnodes
!        write(99,*)oldtonew(i)
!        write(299)oldtonew(i)
!      end do
      
      return
end subroutine ILUK_factorization











!============================ REORDER ==============================
!
! Re-order A using Cuthill-McKee
!
!===================================================================
subroutine reOrder(nnodes,ndim,nnz,ia,ja,iau,ianew,iaunew,janew,list,oldtonew,jnewtojold,ireorder)
  integer, intent(in) :: nnodes,ndim,nnz,ireorder
  integer,intent(in):: ia(nnodes+1),iau(nnodes),ja(nnz)
  integer,intent(inout):: ianew(nnodes+1),iaunew(nnodes),janew(nnz),jnewtojold(nnz)
  integer,intent(inout):: oldtonew(nnodes)
  integer,intent(inout):: list(nnodes)

  integer :: tag(nnodes)
  integer i,j,k,istart,iend,icount,nonzeros,inode,oldnode,oldcol,jcount,newcol
  integer p,kstart,kend,jstart,jend
!write(6,*)"reOrder: checkpt 1"
  !
  do  i = 1,nnodes
     list(i) = i
     iaunew(i) = iau(i)
     ianew(i)  = ia(i)
  end do
!write(6,*)"reOrder: checkpt 2"

  !list(nnodes+1) = nnodes + 1
  ianew(nnodes+1) = ia(nnodes+1)
!write(6,*)"reOrder: checkpt 3"

  do  i = 1,nnz
     janew(i) = ja(i)
  end do
!write(6,*)"reOrder: checkpt 4"

  do  i = 1,nnodes
     oldtonew(i) = i
  end do
!write(6,*)"reOrder: checkpt 5"

  do  i = 1,nnz
     jnewtojold(i) = i
  end do
!write(6,*)"reOrder: checkpt 6"

  !ireorder = 0 
  if(ireorder.ne.1)return
  !
!write(6,*)"reOrder: checkpt 7"
  do  i = 1,nnodes
     tag(i) = 0
  end do

!write(6,*)"reOrder: checkpt 8"

  !
  ! Now go grab all the entries that connect to this node in the number list
  !
  !     write(6,*)"Entering reOrder"
  icount = 1
  istart = 1
  list(1) = 1
  tag(1)  = 1
6000 continue
  jstart = ia(list(istart))
  jend   = ia(list(istart)+ 1) - 1
  do  j = jstart,jend
     inode = ja(j)
     if(tag(inode).eq.0)then
        icount = icount + 1
        list(icount) = inode
        tag(inode) = 1
        if(icount.eq.nnodes)goto 7000
     end if
  end do
  istart = istart + 1
  goto 6000
7000 continue
  !list(nnodes+1) = nnodes + 1

  !     write(6,*)"Exiting reOrder"
  !
  ! At this point, list(i) is a list of the old node
  ! numbers stored in the new order
  ! For example, if list() = {1,3,4,6,5,2}
  ! then list(1) = old node 1
  !      list(2) = old node 3
  !      list(3) = old node 4
  !      list(4) = old node 6
  !      list(5) = old node 5
  !      list(6) = old node 2
  !
  ! Now we can construct new ia, iau, and ja
  !
!write(6,*)"reOrder: checkpt 8"

  ianew(1) = 1
  do  i = 1,nnodes
     oldnode = list(i)
     nonzeros = ia(oldnode + 1) - ia(oldnode)
     ianew(i+1) = ianew(i) + nonzeros
  end do
  !
  ! Fill oldtonew
  do  i = 1,nnodes
     oldtonew(list(i)) = i
  end do
  !
  ! Now get ja
  jcount = 0
  do  i = 1,nnodes
     oldnode = list(i)
     jstart = ia(oldnode)
     jend   = ia(oldnode + 1) - 1
     do  j = jstart,jend
        oldcol = ja(j)
        newcol = oldtonew(oldcol)
        jcount = jcount + 1
        janew(jcount) = newcol
     end do
     istart = ianew(i)
     iend   = ianew(i+1) - 1
     call SORTER(istart,iend,janew,iaunew,i)
  end do
  !
  ! Now get iau
  !
  !     do 8040 i = 1,nnodes
  !       jstart = ianew(i)
  !       jend   = ianew(i+1) - 1
  !       do 8050 j = jstart,jend
  !         if(janew(j).eq.i)iaunew(i) = j
  !8050   continue
  !8040 continue
  !iaunew(nnodes + 1) = -1

!write(6,*)"reOrder: checkpt 9"


  !map jnew to jold
  do i=1,nnodes
     p=oldtonew(i) 
     jstart=ia(i)
     jend=ia(i+1)-1 
     do j=jstart,jend
        oldnode=ja(j) 
        !//now loop over janew to find corresponding node
        kstart=ianew(p)
        kend=ianew(p+1)-1 !//printf("kstart=%d kend=%d\n",kstart,kend);fflush(stdout);//ceb
        do k=kstart,kend
           if (oldnode.eq.list(janew(k))) then
              jnewtojold(k)=j 
              goto 500 !//exit loop
           end if
        end do
500 continue
     end do
  end do



  !
  !     write(6,*)"New ia and iau"
  !     do 9000 i = 1,nnodes+1
  !       write(6,*)ianew(i),iaunew(i)
  !9000 continue
  !
  !     write(6,*)"New ja"
  !     do 9001 i = 1,nnodes
  !       jstart = ianew(i)
  !       jend = ianew(i+1) - 1
  !       write(6,*)(janew(j),j=jstart,jend)
  !9001 continue
  !     stop
  !
  return
end subroutine reOrder


!===================================================================
!
! Sort each of our bins
!
!===================================================================
subroutine SORTER(istart,iend,ja,iau,inode)
  integer, intent(in)::istart,iend,inode
  integer ja(1),iau(1),i,j,min,minsave,jsave
  do  i = istart,iend
     min = ja(i)
     minsave = ja(i)
     jsave = i
     do  j = i+1,iend
        if(ja(j).lt.min)then
           min = ja(j)
           jsave = j
        end if
     end do
     ja(i) = min
     ja(jsave) = minsave
     if(ja(i).eq.inode)iau(inode) = i
  end do
  return
end subroutine SORTER





!============================ BANDWIDTH ================================
!
! Finds the maximum bandwidth
!
!===================================================================
subroutine BANDWIDTH(ia,iau,ja,nnodes,nnz,bw)
  integer,intent(in) :: nnodes,nnz
  integer,intent(in) :: ia(nnodes + 1),iau(nnodes),ja(nnz)
  integer,intent(inout) ::  bw
  integer i,jstart,jend,lowerBand,upperBand

  bw = 0
  do  i = 1,nnodes
     jstart = ia(i)
     jend   = ia(i+1) - 1
     lowerBand  = i - ja(jstart)
     upperBand = ja(jend) - i
     bw = max(lowerBand,upperBand,bw)
  end do
  write(6,*)"Bandwidth = ",bw
  return
end subroutine BANDWIDTH






!c ===================================================================== 
!subroutine readlist(INFILE,OUTFIL,PLTFIL,DRVFIL,DPLTFL,PTBFIL,LSTFIL,LISTF,ITER,NLF)
subroutine readparams(sp)
  
  !character*64 :: matrix_file 
  integer :: mformat, method, precon, nsubiters, krylovdirs,restarts
  integer :: fill_BU, fill_BL, fill_EU, fill_LF, fill_S, fill_SL, fill_SU
  integer ::  block_size, nlev, mreord, rscale_B, cscale_B, rscale_S, cscale_S, pq_S
  real(dp)::convgtol, droptol
  
  !include 'REG6ST.FOR'!defines parameters
  type(solver_params_type),intent(inout) :: sp  
  CHARACTER*64 paramfile  
  CHARACTER*64 LSTFIL
  INTEGER NLF,NITER,ITER,LISTF,ierror
  character(len=130) :: matrix_file

  NAMELIST/INSET/ matrix_file, mformat, method, precon, nsubiters, krylovdirs, &
       restarts, convgtol, droptol, block_size, nlev, mreord,&
       fill_BU, fill_BL, fill_LF, fill_EU, fill_S, fill_SL, fill_SU, &
       rscale_B, cscale_B, rscale_S, cscale_S, pq_S
  !      data PRE1,PRE2,POST/'./dvlists/dvlist','dvlist','.in'/
  data paramfile/'input.in'/
  !call filename1(PRE1,iter,POST,CONV)
  nlf=1
  LSTFIL=paramfile
  OPEN(UNIT=NLF,FILE=LSTFIL,delim='apostrophe',STATUS='OLD',IOSTAT=ierror)
  !open(UNIT=444,FILE=input_file,form='formatted',STATUS='old',IOSTAT=ierror)
  READ(NLF,NML=INSET)
  LISTF=1
  print*,''
  print*,'FOUND input param file ',LSTFIL
  close (unit=NLF)
  
  print*,'input matrix file ',matrix_file


  !this is a messy way to do this but I can't use the sp member vars directly in the namelist
  sp%matrix_file=matrix_file

  print*,'input sp%matrix file ',sp%matrix_file

  sp%mformat= mformat
  sp%method=method
  sp%precon=precon
  sp%nsubiters=nsubiters
  sp%krylovdirs=krylovdirs
  sp%restarts=restarts
  sp%convgtol=convgtol
  sp%droptol=droptol
  sp%block_size=block_size
  sp%nlev=nlev
  sp%mreord=mreord
  sp%fill_BU=fill_BU
  sp%fill_BL=fill_BL
  sp%fill_LF=fill_LF
  sp%fill_EU=fill_EU
  sp%fill_S=fill_S
  sp%fill_SL=fill_SL
  sp%fill_SU=fill_SU
  sp%rscale_B=rscale_B
  sp%cscale_B=cscale_B
  sp%rscale_S=rscale_S
  sp%cscale_S=cscale_S
  sp%pq_S=pq_S
  return
3 print*,''
  print*,'CANNOT FIND input param file ',LSTFIL
  return 
end subroutine readparams
!c ===================================================================== 
























end module linalg







