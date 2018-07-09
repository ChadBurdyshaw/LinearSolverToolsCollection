module pilu_lib
  use kinddefs!, only: dp,ja_type,submat_type
  use blas
  use armslib
!  use matsol_lib
  implicit none
contains

!#define DBL_EPSILON 2.2204460492503131e-16 // double epsilon

!int pilu(p4ptr amat, csptr B, csptr C, double *droptol, int *lfil, csptr schur) 
subroutine pilu(amat, B, C, droptol, lfil, schur) 

  !{
  !/*---------------------------------------------------------------------- 
  !| PARTIAL ILUT -
  !| Converted to C so that dynamic memory allocation may be implememted
  !| in order to have no dropping in block LU factors.
  !|----------------------------------------------------------------------
  !| Partial block ILU factorization with dual truncation. 
  !|               
  !| note: amat is already in form B,F,E,C from reordering process                                     
  !| |  B   F  |        |    L      0  |   |  U   L^{-1} F |
  !| |         |   =    |              | * |               |
  !| |  E   C  |        | E U^{-1}  I  |   |  0       S    |                   
  !|                                                                      
  !| where B is a sub-matrix of dimension B->n.
  !| 
  !|----------------------------------------------------------------------
  !|
  !| on entry:
  !|========== 
  !| ( amat ) = Permuted matrix stored in a PerMat4 struct on entry -- 
  !|            Individual matrices stored in SparRow structs.
  !|            On entry matrices have C (0) indexing.
  !|            on return contains also L and U factors.
  !|            Individual matrices stored in SparRow structs.
  !|            On return matrices have C (0) indexing.
  !|
  !| lfil[1]  =  number nonzeros in L-part
  !| lfil[2]  =  number nonzeros in U-part
  !| lfil[3]  =  number nonzeros in L^{-1} F
  !| lfil[4]  =  not used
  !| lfil[5]  =  number nonzeros in Schur complement
  !|
  !| droptol[1] = threshold for dropping small terms in L during
  !|              factorization.
  !| droptol[2] = threshold for dropping small terms in U.
  !| droptol[3] = threshold for dropping small terms in L^{-1} F during
  !|              factorization.
  !| droptol[4] = threshold for dropping small terms in E U^{-1} during
  !|              factorization.
  !| droptol[5] = threshold for dropping small terms in Schur complement
  !|              after factorization is completed.
  !|
  !| On return:
  !|===========
  !|
  !| (schur)  = contains the Schur complement matrix (S in above diagram)
  !|            stored in SparRow struct with C (0) indexing.
  !|
  !|
  !|       integer value returned:
  !|
  !|             0   --> successful return.
  !|             1   --> Error.  Input matrix may be wrong.  (The 
  !|                         elimination process has generated a
  !|                         row in L or U whose length is > n.)
  !|             2   --> Memory allocation error.
  !|             5   --> Illegal value for lfil or last.
  !|             6   --> zero row in B block encountered.
  !|             7   --> zero row in [E C] encountered.
  !|             8   --> zero row in new Schur complement
  !|----------------------------------------------------------------------- 
  !| work arrays:
  !|=============
  !| jw, jwrev = integer work arrays of length B->n.
  !| w         = real work array of length B->n. 
  !| jw2, jwrev2 = integer work arrays of length C->n.
  !| w2          = real work array of length C->n. 
  !|----------------------------------------------------------------------- 
  !|     All processing is done using C indexing.
  !|--------------------------------------------------------------------*/

  !input parameters
  type(p4_type),intent(inout),pointer :: amat!should already be allocated
  type(cs_type),pointer :: B
  type(cs_type),pointer :: C
  type(cs_type),pointer :: schur
  real(dp),intent(in) :: droptol(:)
  integer,intent(in) :: lfil(:)
  !internal vars
  integer iisub,jjsub,kksub
  integer:: i, ii, j, jj, jcol, jpos, jrow, k, len, len2, lenu, lenl
  integer :: rmax, lsize, rsize,lrowz,rrowz,itmp,idiag
  integer :: test,istat,nsubrow,nsubcol
  integer :: fil_L,fil_U,fil_EU,fil_LF,fil_S,avgB_nnz,avgC_nnz,maxB_nnz,maxC_nnz
  real(dp) :: tmp, tnorm, tabs, tmax, drop0,drop1,drop2,drop3,drop4
  integer, target,dimension(:),allocatable :: jw, jwrev, jw2, jwrev2,lflen
  type(ja_type),target,dimension(:),allocatable :: lfja
  type(submat_type),target,dimension(:),allocatable :: lfma
  real(dp),target,dimension(:,:,:),allocatable :: w,w2
  type(ja_type), pointer :: lrowj, rrowj !pointers to ja type
  type(submat_type), pointer :: lrowm, rrowm !pointers to submat type
  real(dp), dimension(:),allocatable :: diag
  real(dp), dimension(:,:),allocatable:: Uij,Lij,s,t,fact

  real(dp),pointer,dimension(:,:,:) :: shift_w
  integer,pointer,dimension(:) :: shift_jw
  !/*-----------------------------------------------------------------------*/


    !ceb compute avg nnzrow for amat
    avgB_nnz=0
    maxB_nnz=0
    do i=1,B%n
       avgB_nnz=avgB_nnz+B%nnzrow(i)
       if(B%nnzrow(i)>maxB_nnz) maxB_nnz=B%nnzrow(i)
    end do
    avgB_nnz=avgB_nnz/B%n
    !avgB_nnz=maxB_nnz

    avgC_nnz=0
    maxC_nnz=0
    do i=1,C%n
       avgC_nnz=avgC_nnz+C%nnzrow(i)
       if(C%nnzrow(i)>maxC_nnz) maxC_nnz=C%nnzrow(i)
if(C%nnzrow(i).eq.0) write(6,*)"C%nnzrow(",i,")=",C%nnzrow(i)
    end do
    avgC_nnz=avgC_nnz/C%n
    !avgC_nnz=maxC_nnz

    !ceb define fill (number of nonzeros allowed) = [1+lfil()]*avg_nnz
    !     so that lfil=0 results in fil=avg_nnz per row
    !  this allows us to avoid having to change lfil as the size of the input matrix changes (and as we descend levels)
    !ceb though currently defined as int lfil could be real so that fil_() does not have to be multiple of avg_nnz
  fil_L  = int((1+lfil(1))*avgB_nnz)  
  fil_U  = int((1+lfil(2))*avgB_nnz) 
  fil_EU = int((1+lfil(3))*avgB_nnz)
  fil_LF = int((1+lfil(4))*avgB_nnz)
  fil_S  = int((1+lfil(5))*avgC_nnz) 

    write(6,*) "pilu: fil_L=",fil_L," fil_U=",fil_U," avgB_nnz=",avgB_nnz
    write(6,*) "pilu: fil_EU=",fil_EU," fil_LF=",fil_LF
    write(6,*) "pilu: fil_S=",fil_S,"avgC_nnz=",avgC_nnz

!write(6,*) "fil[0,1,2,4]=",fil_0,fil_1,fil_2,fil_4 

  drop0=droptol(1)
  drop1=droptol(2)
  drop2=droptol(3)
  drop3=droptol(4)
  drop4=droptol(5)
!write(6,*)"droptol[1,2,3,4,5]=",drop0,drop1,drop2,drop3,drop4

  nsubrow=amat%nsubrow
  nsubcol=amat%nsubcol

  !temporary vars
  allocate(diag(nsubrow),STAT=istat)
  allocate(s(nsubrow,nsubcol),STAT=istat)
  allocate(t(nsubrow,nsubcol),STAT=istat)
  allocate(fact(nsubrow,nsubcol),STAT=istat)
  allocate(Lij(nsubrow,nsubcol),STAT=istat)
  allocate(Uij(nsubrow,nsubcol),STAT=istat)

  lsize = amat%nB
  rsize = C%n
  rmax = rsize
  if (lsize>rsize) rmax=lsize

write(6,*)"calling pilu"


  test=0!ceb


  allocate(jw(rmax),STAT=istat) !jw = (int *) Malloc(rmax*sizeof(int), "pilu:1" );
  allocate(w(nsubrow,nsubcol,rmax),STAT=istat) !w = (double *) Malloc(rmax*sizeof(double), "pilu:2" );
  allocate(jwrev(rmax),STAT=istat)!jwrev = (int *) Malloc(rmax*sizeof(int), "pilu:3" );
  allocate(jw2(rmax),STAT=istat)  !jw2 = (int *) Malloc(rmax*sizeof(int), "pilu:4" );
  allocate(w2(nsubrow,nsubcol,rmax),STAT=istat)  !w2 = (double *) Malloc(rmax*sizeof(double), "pilu:5" );
  allocate(jwrev2(rmax),STAT=istat)!jwrev2 = (int *) Malloc(rmax*sizeof(int), "pilu:6" );
  if (fil_L < 0 .or. fil_U<0 .or. amat%L%n<=0) goto 9995

  allocate(lfma(lsize),STAT=istat)!lfma = (double **) Malloc(lsize*sizeof(double *), "pilu:7" ); 
  allocate(lfja(lsize),STAT=istat)!lfja = (int **) Malloc(lsize*sizeof(int *), "pilu:8" ); 
  allocate(lflen(rmax),STAT=istat)!lflen = (int *) Malloc(lsize*sizeof(int), "pilu:9" ); 

!write(6,*)"pilu checkpt after temp array allocation"


!ceb--
!  do i=1,B%n
!     do j=1,B%nnzrow(i)
!  !      write(94,*) i, B%pj(i)%cols(j)
!        write(94,*) i, B%pj(i)%cols(j),((B%pa(i)%submat(iisub,jjsub,j),iisub=1,amat%nsubrow),jjsub=1,amat%nsubcol)
!!        write(94,*) i, B%pd(i), B%pj(i)%cols(j)
!     end do
!  end do
!!  !write(6,*)""
!  do i=1,C%n
!     do j=1,C%nnzrow(i)
!        !write(96,*) i,C%pj(i)%cols(j)
!        write(96,*) i, C%pj(i)%cols(j),((C%pa(i)%submat(iisub,jjsub,j),iisub=1,amat%nsubrow),jjsub=1,amat%nsubcol)
!!        write(96,*) i,C%pd(i),C%pj(i)%cols(j)
!     end do
!  end do
!
!  do i=1,amat%F%n
!     do j=1,amat%F%nnzrow(i)
!        !write(98,*) i, amat%F%pj(i)%cols(j)
!        !write(98,*) i, amat%F%pj(i)%cols(j),((amat%F%pa(i)%submat(iisub,jjsub,j),iisub=1,amat%nsubrow),jjsub=1,amat%nsubcol)
!        write(98,*) i, amat%F%pj(i)%cols(j),((amat%F%pa(i)%submat(iisub,jjsub,j),iisub=1,amat%nsubrow),jjsub=1,amat%nsubcol)
!     end do
!  end do
!
!  do i=1,amat%E%n
!     do j=1,amat%E%nnzrow(i)
!        !write(99,*) i, amat%E%pj(i)%cols(j)
!        !write(99,*) i, amat%E%pj(i)%cols(j),((amat%E%pa(i)%submat(iisub,jjsub,j),iisub=1,amat%nsubrow),jjsub=1,amat%nsubcol)
!        write(99,*) i, amat%E%pj(i)%cols(j),((amat%E%pa(i)%submat(iisub,jjsub,j),iisub=1,amat%nsubrow),jjsub=1,amat%nsubcol)
!     end do
!  end do
!!--ceb


  do j=1,rmax
     jwrev(j) = -1
     jwrev2(j) = -1
  end do


!write(6,*)"pilu checkpt prior to main loop"

  !/*---------------------------------------------------------------------
  !|    beginning of first main loop - L, U, L^{-1}F calculations
  !|--------------------------------------------------------------------*/
  do ii=1,lsize
!write(6,*) "pilu: computing row ",ii
     lrowj => B%pj(ii)!lrowj is pointer
     lrowm => B%pa(ii)!lrowm is pointer
     lrowz =  B%nnzrow(ii)

     rrowj => amat%F%pj(ii)!rrowj is pointer
     rrowm => amat%F%pa(ii)!rrowm is pointer
     rrowz =  amat%F%nnzrow(ii)
     !/*---------------------------------------------------------------------
     !|   check for zero row in B block
     !|--------------------------------------------------------------------*/
     tnorm = 0.0 
     do k=1,lrowz
        tnorm = tnorm + L1Norm2(lrowm%submat,k,nsubrow) !tnorm = tnorm + abs(rowm(k)) 
     end do
     tnorm = tnorm / lrowz 
!if(tnorm>10.) write(6,*) "pilu: row ",ii,"tnorm=",tnorm
     if(abs(tnorm-DBL_EPSILON) .le. DBL_EPSILON*tnorm) goto  9996

     !/*---------------------------------------------------------------------
     !|     unpack B-block in arrays w, jw, jwrev
     !|     WE ASSUME THERE IS A DIAGONAL ELEMENT
     !|--------------------------------------------------------------------*/
     lenu = 1
     lenl = 0!lenl does not include diagonal
     jw(ii) = ii
     jwrev(ii) = ii
     call blank_3D_submat(w,ii,nsubrow,nsubcol) !w(ii) = 0.0!zero submatrix
     do j=1, lrowz
        jcol = lrowj%cols(j)
        if (jcol .lt. ii) then !left of diag
           lenl=lenl+1
           jpos=lenl
           jw(jpos) = jcol 
           jwrev(jcol) = jpos!ceb
        else if (jcol .eq. ii) then !on diag
           jpos=ii
        else                        !right of diag
           jpos = ii+lenu
           jw(jpos) = jcol
           jwrev(jcol) = jpos
           lenu=lenu+1
        end if
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              w(iisub,jjsub,jpos)=lrowm%submat(iisub,jjsub,j)
           end do
        end do  
     end do
     !/*---------------------------------------------------------------------
     !|     unpack F-block in arrays w2, jw2, jwrev2 
     !|     (all entries are in U portion)
     !|--------------------------------------------------------------------*/
     len2 = 0
!     if (test.eq.1 .and. C%n.ne.0) then
        do j=1,rrowz
           len2=len2+1
           jpos=len2
           jcol = rrowj%cols(j)
           jw2(jpos) = jcol
           jwrev2(jcol) = jpos
           !w2(len2) = rrowm(j)!copy submatrix
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 w2(iisub,jjsub,jpos)=rrowm%submat(iisub,jjsub,j)
              end do
           end do
        end do
!     end if!test
     
     !/*---------------------------------------------------------------------
     !|     Eliminate previous rows -  
     !|--------------------------------------------------------------------*/
     len = 0
     jj=0!ceb
     !loop over block column entries in row ii of L
     !====================================================================
     do !jj=1,lenl
        jj=jj+1
        if (jj.gt.lenl) exit
        !/*---------------------------------------------------------------------
        !|    in order to do the elimination in the correct order we must select
        !|    the smallest column index among jw(k), k=jj+1, ..., lenl.
        !|--------------------------------------------------------------------*/
        jrow = jw(jj)
        k = jj
        !/*---------------------------------------------------------------------
        !|     determine smallest column index
        !|--------------------------------------------------------------------*/
        do j=jj+1,lenl
           if (jw(j) < jrow) then
              jrow = jw(j)
              k = j
           end if
        end do
        if (k .ne. jj) then 
           !write(6,*)"1: pilu row ",ii," exchanging w col ",k," with ",jj   
           !/*   exchange in jw   */
           j = jw(jj)
           jw(jj) = jw(k)
           jw(k) = j
           !/*   exchange in jwrev   */
           jwrev(jrow) = jj
           jwrev(j) = k
           !/*   exchange in w   */
           call swap_submat(w,jj,k,nsubrow)
        end if
        !/*---------------------------------------------------------------------
        !|     zero out element in row.
        !|--------------------------------------------------------------------*/
        jwrev(jrow) = -1
        !/*---------------------------------------------------------------------
        !|     get the multiplier for row to be eliminated (jrow).
        !|--------------------------------------------------------------------*/
        lrowm => amat%U%pa(jrow)
        lrowj => amat%U%pj(jrow)
        lrowz = amat%U%nnzrow(jrow)
        !fact = w(jj) * U(0)!submatrix mult
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              fact(iisub,jjsub)=0.0
              do kksub=1,nsubrow !assumes square submatrix
                 fact(iisub,jjsub) = fact(iisub,jjsub) + w(iisub,kksub,jj)*lrowm%submat(kksub,jjsub,1) !equivalent L(ii,jj)*U(diagj) (pivot)
              end do
           end do
        end do
        
!if(L1Norm(fact,nsubrow)>10)write(6,*)"pilu: ii:",ii," pivot norm(",jrow,")=",L1Norm(fact,nsubrow)

        if (L1Norm(fact,nsubrow) <= drop0 ) then !/*   DROPPING IN L   */
!write(6,*)"pilu: dropping norm ",L1Norm(fact,nsubrow)," in L  row,column(",ii,",",jrow,")"
           cycle
        end if
        
        !/*---------------------------------------------------------------------
        !|     combine current row and row jrow
        !|--------------------------------------------------------------------*/
        do k=2,lrowz
           !s = fact * rowm(k)
           jcol = lrowj%cols(k)
           jpos = jwrev(jcol)  
           if (jpos == -1) then 
              !/*---------------------------------------------------------------------
              !|     this is a fill-in element
              !|--------------------------------------------------------------------*/
              if (jcol >= ii) then ! right of diag
                 !/*---------------------------------------------------------------------
                 !|     dealing with U
                 !|--------------------------------------------------------------------*/
                 if(lenu < fil_U) then
                 !idiag=B%pd(ii)
                 !if(idiag<=0)idiag=1
                 !if(lenu < fil_U+(B%nnzrow(ii)-idiag+1) ) then
                    jpos= ii+ lenu
                    jw(jpos) = jcol
                    jwrev(jcol) = jpos
                    lenu=lenu+1
                    call blank_3D_submat(w,jpos,nsubrow,nsubcol)
                 end if
              else ! left of diag
                 !/*---------------------------------------------------------------------
                 !|     dealing  with L
                 !|--------------------------------------------------------------------*/ 
                 !if(lenl < fil_L) then
                 !idiag=B%pd(ii)
                 !if(idiag<=0)idiag=1
                 !if(lenl < fil_L+(B%nnzrow(ii)-idiag+1) ) then
                 !if(lenl < fil_L+(idiag-1)) then
                 if(lenl < fil_L) then
                    lenl=lenl+1
                    jpos=lenl
                    jw(jpos) = jcol
                    jwrev(jcol) = jpos
                    call blank_3D_submat(w,jpos,nsubrow,nsubcol)
                 end if
              end if! if fill in left or right of diag
           !else
              !/*---------------------------------------------------------------------
              !|     this is not a fill-in element 
              !|--------------------------------------------------------------------*/
           end if!if fill in element
  
           !perform row operation w(ii)-=pivot*U(k)
           if(jpos.ne.-1) then
              do iisub=1,nsubrow
                 do jjsub=1,nsubcol
                    do kksub=1,nsubrow
                       w(iisub,jjsub,jpos) = w(iisub,jjsub,jpos) -  fact(iisub,kksub)*lrowm%submat(kksub,jjsub,k)
                    end do
                 end do
              end do
           end if

           !if ((lenu .gt. lsize) .or. (lenl .gt. lsize)) goto 9991 !cycle?
           if ((lenu .gt. rmax) .or. (lenl .gt. rmax)) goto 9991 !cycle?
        end do
        
        !/*---------------------------------------------------------------------
        !|     store this pivot element
        !|--------------------------------------------------------------------*/
        len=len+1
        jpos=len
        jw(jpos) = jrow
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              w(iisub,jjsub,jpos)= fact(iisub,jjsub)
           end do
        end do        

!ceb this will have to be modified to include L portion in block diagonal
!if (test.eq.1 .and. C%n.ne.0) then  

        rrowj => lfja(jrow)
        rrowm => lfma(jrow)
        rrowz = lflen(jrow)      
        !/*---------------------------------------------------------------------
        !|     dealing  with  L^{-1} F
        !|--------------------------------------------------------------------*/
        do k=1,rrowz 
           !s = fact * rrowm[k];
           jcol = rrowj%cols(k)
           jpos = jwrev2(jcol)
           !/*---------------------------------------------------------------------
           !|     this is a fill-in element
           !|--------------------------------------------------------------------*/
           if (jpos == -1) then
              if(len2<fil_LF) then
              !if(len2<fil_2+amat%F%nnzrow(ii)) then
                 len2=len2+1
                 jpos=len2
                 jw2(jpos) = jcol
                 jwrev2(jcol) = jpos
                 call blank_3D_submat(w2,jpos,nsubrow,nsubcol)
              end if
              !/*---------------------------------------------------------------------
              !|     this is not a fill-in element
              !|--------------------------------------------------------------------*/
           !else
              !   jpos = jwrev2(jcol)!ceb as above
           end if
           
           if(jpos.ne.-1) then 
              do iisub=1,nsubrow
                 do jjsub=1,nsubcol
                    do kksub=1,nsubrow
                       w2(iisub,jjsub,jpos) = w2(iisub,jjsub,jpos) - fact(iisub,kksub)*rrowm%submat(kksub,jjsub,k)
                    end do
                 end do
              end do
           end if           
        end do
        
!end if!test            
     end do!end loop over block cols

!write(6,*)"calling pilu checkpt 1 end col loop"

     !/*---------------------------------------------------------------------
     !|     reset nonzero indicators
     !|--------------------------------------------------------------------*/
     do j=1,len2    !/*  L^{-1} F block  */
        jwrev2(jw2(j)) = -1
     end do
     do j=1,lenl    !/*  L block  */
        jwrev(jw(j)) = -1
     end do
     do j=1,lenu    !/*  U block  */
        jwrev(jw(ii+j-ofst)) = -1
     end do

     !/*---------------------------------------------------------------------
     !|     done reducing this row, now store L
     !|--------------------------------------------------------------------*/
     !lenl=len
     len=lenl
     if (lenl> fil_L) len= fil_L !lenl = len > fil0 ? fil0 : len;
     amat%L%nnzrow(ii) = len

     !/* quick sort */
!write(6,*)"pilu qsplit checkpt 1"
     if (lenl > len) call qsplit(w,jw,lenl,len,nsubrow)!quick sort alg (sorts based on L1norm of submats)
     if (len > 0) then
!write(6,*)"allocate L(",ii,",",lenl,")"
        allocate(amat%L%pj(ii)%cols(len),STAT=istat)
        allocate(amat%L%pa(ii)%submat(nsubrow,nsubcol,len),STAT=istat)
        do k=1,len
           amat%L%pj(ii)%cols(k)=jw(k)!copy jw into new U_ja elements
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 amat%L%pa(ii)%submat(iisub,jjsub,k)=w(iisub,jjsub,k)!copy w into new U elements, U[k]=w[ii+(k-1)]
              end do
           end do
        end do
     end if


     !/*---------------------------------------------------------------------
     !|     store the diagonal element of U
     !|     dropping in U if size is less than drop1 * diagonal entry
     !|--------------------------------------------------------------------*/


     !/*---------------------------------------------------------------------
     !|     apply dropping strategy to U  (entries after diagonal)
     !|--------------------------------------------------------------------*/
!ceb: This following section only works when jpos=ii+len not jpos=len as in the arms source code
!  this is because the iith diagonal starts in w(ii), keeping L to the left and U to the right
     tnorm=L1Norm2(w,ii,nsubrow)!ceb
!if(tnorm>10.0) write(6,*)"pilu checkpoint 30 LU tnorm(",ii,")=",tnorm

     len = 0
     do j=2,lenu !if j col element is dropped len keeps track of what col entry we are filling
        if ( L1Norm2(w,ii+j-ofst,nsubrow) > drop1*tnorm ) then !use matrix norm
           len=len+1
           jpos=ii+len
           jw(jpos) = jw(ii+j-ofst)
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 w(iisub,jjsub,jpos)= w(iisub,jjsub,ii+j-ofst)
             end do
           end do
        else
           !write(6,*) "pilu checkpoint 31 :Dropping small column element  norm(w(", ii+j-ofst, "))=", L1Norm2(w,ii+j-ofst,nsubrow),&
           !     " < tol*tnorm(",ii,")=",drop1*tnorm
        end if
     end do

     lenu=len+1!assign new length after dropping
     len=lenu!ceb
     if (len > fil_U) len=fil_U !lenu = len+1 > fil_1 ? fil_1 : len+1;
     amat%U%nnzrow(ii) = len
!write(6,*)"pilu: total B_LU(",ii,") nnzrow=",amat%U%nnzrow(ii)+amat%L%nnzrow(ii)
     jpos = lenu !jpos = lenu - 1
!write(6,*)"pilu qsplit checkpt 2"
     if (jpos < len) call qsplit(w,jw,len,jpos,nsubrow)!quick sort alg (sorts based on L1norm of submats)

     if (len>0) then
!write(6,*)"allocate U(",ii,",",lenu,")"
        allocate(amat%U%pa(ii)%submat(nsubrow,nsubcol,len),STAT=istat)
        allocate(amat%U%pj(ii)%cols(len),STAT=istat)
     end if



     !/*---------------------------------------------------------------------
     !|     copy the rest of U (all but diagonal)
     !|--------------------------------------------------------------------*/
     lenu = 1!diagonal already exists !ceb
     do k=2,len! (skips diagonal)
          lenu=lenu+1 !ceb
        amat%U%pj(ii)%cols(k)=jw(ii+k-ofst)!copy jw into new U_pj elements
        do iisub=1,nsubrow
           do jjsub=1,nsubrow
               amat%U%pa(ii)%submat(iisub,jjsub,k)=w(iisub,jjsub,ii+k-ofst)!copy w into new U elements, U[k]=w[ii+(k-1)]
           end do
        end do
     end do

     tnorm=L1Norm2(w,ii,nsubrow)
     if(abs(tnorm-DBL_EPSILON) <= DBL_EPSILON*tnorm) then
        do iisub=1,nsubrow
           w(iisub,iisub,ii) = (0.0001+drop0)*tnorm 
        end do
     end if



!ceb-- replaced diagonal section
     !===========================================================
     !//compute L,U for pivot (diagonal of combined LU matrix)
     ! Lij,Uij submatrices refer only to the diagonal block element
     tmp=0.0
     !diagp=iau(ii)!need diag of B(ii)
     do jjsub=1,nsubrow
        do kksub=1,nsubcol
           if (jjsub.eq.kksub) then!//set diagonal (sub diagonal)
              tmp=w(jjsub,jjsub,ii)!A[p,p](jj,jj)
              do iisub=1,jjsub-1
                 tmp = tmp - Lij(iisub,jjsub)*Uij(jjsub,iisub) ! A[p,p](jj,jj)=A[p,p](jj,jj)- L[pj](ii,jj)*U[pj](jj,ii)
              end do
              Lij(kksub,kksub)= 1.0/tmp !inverse diag
           else if(kksub<jjsub) then!//set lower
              tmp = w(jjsub,kksub,ii)!A[p,p](jjsub,kksub) L 
              do iisub=1,kksub-1
                 tmp = tmp - Lij(iisub,jjsub)*Uij(kksub,iisub)
              end do
              Lij(kksub,jjsub)=tmp
           else if(kksub>jjsub) then !//set upper
              tmp = w(jjsub,kksub,ii)!A[p,p](jjsub,kksub) U
              do iisub=1,jjsub-1
                 tmp = tmp - Lij(iisub,jjsub)*Uij(kksub,iisub)
              end do
              Uij(kksub,jjsub) = tmp*Lij(jjsub,jjsub)
           end if
        end do
     end do
     
     !//set alu[diagp] to identity
     do jjsub=1,nsubrow
        do kksub=1,nsubcol
           amat%U%pa(ii)%submat(jjsub,kksub,1)=0.0
        end do
        amat%U%pa(ii)%submat(jjsub,jjsub,1)=1.0
     end do
     
     !//write LU into U[diagp]
     do kksub=1,nsubcol !loop subcols
        do jjsub=1,nsubcol-1 !get D[1,..,neqn-1] ! loop subcols
           diag(jjsub)=amat%U%pa(ii)%submat(jjsub,kksub,1) !Identity matrix row jj 
           do iisub=1,jjsub-1
              diag(jjsub) = diag(jjsub) - Lij(iisub,jjsub)*diag(iisub)
           end do
           diag(jjsub) = diag(jjsub)*Lij(jjsub,jjsub)!ceb not sure if this is right
        end do
        
        do jjsub=nsubcol,1,-1 !get alu[diagp,jj,kk] ! loop subcols
           if(jjsub.eq.nsubcol) then     
              !D=D-L*D (left of subdiag)
              do iisub=1,jjsub-1
                 amat%U%pa(ii)%submat(jjsub,kksub,1) = amat%U%pa(ii)%submat(jjsub,kksub,1)  - Lij(iisub,jjsub)*diag(iisub)
              end do
              amat%U%pa(ii)%submat(jjsub,kksub,1) = amat%U%pa(ii)%submat(jjsub,kksub,1)*Lij(jjsub,jjsub)
              !D=
           else      
              !D=D-U*D (right of sub diag)
              amat%U%pa(ii)%submat(jjsub,kksub,1)=diag(jjsub)
              do iisub=jjsub+1,nsubrow
                 amat%U%pa(ii)%submat(jjsub,kksub,1) = amat%U%pa(ii)%submat(jjsub,kksub,1) - &
                      Uij(iisub,jjsub)*amat%U%pa(ii)%submat(iisub,kksub,1)
              end do
           end if
        end do
     end do
   
     !===========================================================
     amat%U%pj(ii)%cols(1) = ii!ceb moved from prior to diag formation
     !ceb may have to include computation of [E U^{-1}] and  [L^{-1} F] here since U(ii,1)
     !now contains portions of L and diag 
     !===========================================================


     !write(6,*)"pilu checkpoint 38"
     
     !/*---------------------------------------------------------------------
     !|     copy  L^{-1} F
     !|--------------------------------------------------------------------*/
     len = 0
     !write(6,*)"pilu checkpoint 38.1 len2=",len2
     do j=1,len2
        !if(ii.eq.1 ) write(6,*)"checkpt 38.9: j=",j
        if ( L1Norm2(w2,j,nsubrow) > drop2*tnorm ) then !use matrix norm
           !w(len) = w2(j);!copy submatrix
           !overwrites w,jw
           len=len+1
           jpos=len
           jw(jpos) = jw2(j)
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 w(iisub,jjsub,jpos)= w2(iisub,jjsub,j)
              end do
           end do
        else
           !write(6,*)" Dropped element ",j," in inv(L)F norm=",L1Norm2(w2,j,nsubrow)," tol=drop2*tnorm=",drop2*tnorm
        end if
     end do
     
     !ceb is lenu right here?
     lenu=len!+ofst!ceb
     !if(lenu>fil_LF+amat%F%nnzrow(ii)) lenu=fil_LF+amat%F%nnzrow(ii) !lenu = len > fil2 ? fil2 : len;
     !if(lenu>fil_LF+amat%F%nnzrow(ii)) lenu=fil_LF+amat%F%nnzrow(ii) !lenu = len > fil2 ? fil2 : len;
     if(lenu>fil_LF) lenu=fil_LF !lenu = len > fil2 ? fil2 : len;
     
     
     !if(lenu.eq.0) write(6,*)"pilu checkpoint 38.1 ii=",ii," lenu=",lenu
     !write(6,*)"pilu qsplit checkpt 3"
     if (len > lenu) call qsplit(w,jw,len,lenu,nsubrow)!quick sort alg (sorts based on L1norm of submats)
     !if (len > lenu) call qsplit(w,jw,len+ofst,lenu,nsubrow)!quick sort alg (sorts based on L1norm of submats)
     lflen(ii) = lenu
     if (lenu > 0) then
        !write(6,*)"allocate lf(",ii,",",lenu,")"
        allocate(lfja(ii)%cols(lenu),STAT=istat)
        allocate(lfma(ii)%submat(nsubrow,nsubcol,lenu),STAT=istat)
        !copy elements
        do k=1,lenu 
           lfja(ii)%cols(k)=jw(k)!copy jw into new U_pj elements
           !write(6,*)"checkpt 38.9: lenu=",lenu," lfja(ii)%cols(",k,")=",jw(k)
           do iisub=1,nsubrow
              do jjsub=1,nsubrow
                 lfma(ii)%submat(iisub,jjsub,k)=w(iisub,jjsub,k)!copy w into new U elements, U[k]=w[ii+(k-1)]
              end do
           end do
        end do
     end if
          
  end do !end loop over rows

  !write(6,*)"pilu checkpoint 39"







     
  !write(6,*)"pilu checkpoint 40 begin second main loop EU calc"
  
  !/*---------------------------------------------------------------------
  !|    beginning of second main loop   E U^{-1} and Schur complement
  !|--------------------------------------------------------------------*/
  do ii=1,rsize
     !write(6,*)"pilu EU computing row ",ii
     lrowj => amat%E%pj(ii)
     lrowm => amat%E%pa(ii)
     lrowz =  amat%E%nnzrow(ii)
     
     rrowj => C%pj(ii)
     rrowm => C%pa(ii)
     rrowz =  C%nnzrow(ii)
     
     !/*---------------------------------------------------------------------
     !|    determine if there is a zero row in [ E C ]
     !|--------------------------------------------------------------------
     !do k=1,lrowz
     !   if (L1Norm2(lrowm%submat,k,nsubrow) .ne. 0.0) then
     !      write(6,*) "Error: zero row in E"
     !      goto 42
     !   end if
     !end do
     do k=1,rrowz
        if (L1Norm2(rrowm%submat,k,nsubrow) .ne. 0.0) then
           goto 42
        end if
     end do
     goto 9996
42   continue
     !*/
     
     !/*---------------------------------------------------------------------
     !|     unpack E in arrays w, jw, jwrev
     !|--------------------------------------------------------------------*/
     lenl = 0
     do j=1,lrowz
        lenl=lenl+1
        jpos=lenl
        jcol = lrowj%cols(j)
        jw(jpos) = jcol
        jwrev(jcol) = jpos
        !w(lenl) = lrowm(j)!copy submatrix
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              w(iisub,jjsub,jpos)= lrowm%submat(iisub,jjsub,j)
           end do
        end do
     end do
     
     !/*---------------------------------------------------------------------
     !|     unpack C in arrays w2, jw2, jwrev2    
     !|--------------------------------------------------------------------*/
     lenu = 0
     do j=1,rrowz
        lenu=lenu+1
        jpos=lenu!+ofst !ceb
        jcol = rrowj%cols(j)
        jw2(jpos) = jcol
        jwrev2(jcol) = jpos
        !w2(lenu) = rrowm(j)!copy submatrix
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              w2(iisub,jjsub,jpos)= rrowm%submat(iisub,jjsub,j)
           end do
        end do
     end do
     
     !write(6,*)"pilu checkpoint 40.1: lenu(",ii,")=",lenu
     
     !/*---------------------------------------------------------------------
     !|     eliminate previous rows
     !|--------------------------------------------------------------------*/
     !write(6,*)"pilu checkpoint 41"
     len = 0
     jj=0
     do !jj=1,lenl
        jj=jj+1
        if (jj.gt.lenl) exit
        !/*---------------------------------------------------------------------
        !|    in order to do the elimination in the correct order we must select
        !|    the smallest column index among jw(k), k=jj+1, ..., lenl.
        !|--------------------------------------------------------------------*/
        jrow = jw(jj)
        k = jj
        !/*---------------------------------------------------------------------
        !|     determine smallest column index
        !|--------------------------------------------------------------------*/
        do j=jj+1,lenl
           if (jw(j) < jrow) then
              jrow = jw(j)
              k = j
           end if
        end do
        if (k .ne. jj) then    
           !/*   exchange in jw   */
           j = jw(jj)
           jw(jj) = jw(k)
           jw(k) = j
           !/*   exchange in jwrev   */
           jwrev(jrow) = jj
           jwrev(j) = k
           !write(6,*)"2: pilu row ",ii," exchanging w col ",k," with ",jj   
           !/*   exchange in w   */
           call swap_submat(w,jj,k,nsubrow)
        end if
        
        !/*---------------------------------------------------------------------
        !|     zero out element in row.
        !|--------------------------------------------------------------------*/
        jwrev(jrow) = -1
        
        
        !/*---------------------------------------------------------------------
        !|     get the multiplier for row to be eliminated (jrow).
        !|--------------------------------------------------------------------*/
        lrowm => amat%U%pa(jrow)
        lrowj => amat%U%pj(jrow)
        lrowz =  amat%U%nnzrow(jrow)
        
        !fact = w(jj) * lrowm(0)!submatrix multiply
        !ceb this factor actually includes some L contribution 
        !  from both the lower half of the block diagonal and the sub diagonal itself
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              fact(iisub,jjsub)=0.
              do kksub=1,nsubrow !assumes square submatrix
                 fact(iisub,jjsub) = fact(iisub,jjsub) + w(iisub,kksub,jj)*lrowm%submat(kksub,jjsub,1)
              end do
           end do
        end do
        
        !write(6,*)"EU: pivot(",ii,jj,")=",fact(1,1)
        
        if (L1Norm(fact,nsubrow) <= drop3 ) then !/*  DROPPING IN E U^{-1}   */
           !              write(6,*)"pilu: DROPPING norm ",L1Norm(fact,nsubrow)," IN E U^{-1}: row,column(",ii,",",jrow,")"
           cycle
        end if
        
        rrowj => lfja(jrow)
        rrowm => lfma(jrow)
        rrowz =  lflen(jrow)
        
        !/*---------------------------------------------------------------------
        !|     combine current row and row jrow   -   first  E U^{-1}
        !|--------------------------------------------------------------------*/
        do k=2,lrowz !(U off diags only)
           !s = fact * lrowm(k)!submatrix multiply
           jcol = lrowj%cols(k)
           jpos = jwrev(jcol)
           
           if (jpos .eq. -1) then
              !/*---------------------------------------------------------------------
              !|     fill-in element
              !|--------------------------------------------------------------------*/
              if (lenl < fil_EU) then
                 !if (lenl < fil_EU+amat%E%nnzrow(ii)) then
                 lenl=lenl+1
                 jpos=lenl
                 jw(jpos) = jcol
                 jwrev(jcol) = jpos
                 call blank_3D_submat(w,jpos,nsubrow,nsubcol)
              end if
              !else
              !/*---------------------------------------------------------------------
              !|     this is not a fill-in element 
              !|--------------------------------------------------------------------*/
           end if
           
           if(jpos.ne.-1) then
              do iisub=1,nsubrow
                 do jjsub=1,nsubcol
                    do kksub=1,nsubrow
                       w(iisub,jjsub,jpos) = w(iisub,jjsub,jpos) -  fact(iisub,kksub)*lrowm%submat(kksub,jjsub,k)
                    end do
                 end do
              end do
           end if
           
           !write(6,*) "pilu checkpt 11"
           !if (lenl .gt. lsize) goto 9991 !cycle?
           !if (lenl .gt. lsize) goto 9991 !cycle?
           if (lenl .gt. rmax) goto 9991 !cycle?
        end do
        
        
        !/*---------------------------------------------------------------------
        !|     incorporate into Schur complement   C - (E U^{-1}) (L^{-1} F)
        !|--------------------------------------------------------------------*/
        do k=1,rrowz
           !s = fact * rrowm(k)!submatrix multiply
           jcol = rrowj%cols(k)
           jpos = jwrev2(jcol)
           !/*---------------------------------------------------------------------
           !|     this is not a fill-in element 
           !|--------------------------------------------------------------------*/
           if (jpos == -1) then
              !if(lenu < fil_S+C%nnzrow(ii) ) then
              if(lenu < fil_S) then
                 !w2(jpos) = -s!subtract submatrix
                 lenu=lenu+1 
                 jpos=lenu
                 jw2(jpos) = jcol
                 jwrev2(jcol) = jpos
                 call blank_3D_submat(w2,jpos,nsubrow,nsubcol)
              end if
              !else
              !/*---------------------------------------------------------------------
              !|     this is not a fill-in element
              !|--------------------------------------------------------------------*/
           end if
           
           if(jpos.ne.-1) then
              do iisub=1,nsubrow
                 do jjsub=1,nsubcol
                    do kksub=1,nsubrow
                       w2(iisub,jjsub,jpos) = w2(iisub,jjsub,jpos) - fact(iisub,kksub)*rrowm%submat(kksub,jjsub,k)
                       !write(6,*)"EULF(",ii,jpos,")=",fact(iisub,kksub)*rrowm%submat(kksub,jjsub,k)
                       
                    end do
                 end do
              end do
           end if
           !write(6,*) "pilu checkpt 10"
           !           if (lenu .gt. rsize) goto 9991 !cycle?
        end do
        
        !/*---------------------------------------------------------------------
        !|     store this pivot element
        !|--------------------------------------------------------------------*/
        !w(len) = fact!copy submatrix
        len=len+1
        jpos=len!ceb
        jw(jpos) = jrow
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              w(iisub,jjsub,jpos)= fact(iisub,jjsub)
           end do
        end do
        !if(jw(jpos).eq.ii) write(6,*)"checkpt 43.8  jw(",ii,"),(",jpos,")=",jw(jpos)
     end do
     
     !write(6,*)"pilu checkpoint 43.9: lenu(",ii,")=",lenu
     
     !/*---------------------------------------------------------------------
     !|     reset nonzero indicators
     !|--------------------------------------------------------------------*/
     do j=1,lenu    !/*  Schur complement  */
        jwrev2(jw2(j)) = -1
     end do
     
     do j=1,lenl   !/*  E U^{-1} block  */
        jwrev(jw(j)) = -1
     end do
     
     !/*---------------------------------------------------------------------
     !|     done reducing this row, now throw away row of E U^{-1}
     !|     and apply a dropping strategy to the Schur complement.
     !|
     !|     store the diagonal element of Schur Complement
     !|
     !|     drop in Schur complement if size less than drop4*tnorm
     !|     where tnorm is the size of the maximum entry in the row
     !|--------------------------------------------------------------------*/
     
     !write(6,*)"pilu checkpoint 44: lenu(",ii,")=",lenu
     !write(6,*)"pilu checkpoint 44.5: C%nnzrow(",ii,")=",C%nnzrow(ii)
     !if(ii.le.lsize)write(6,*)"pilu checkpoint 44.6: amat%F%nnzrow(",ii,")=",amat%F%nnzrow(ii)
     !write(6,*)"pilu checkpoint 44.7: amat%E%nnzrow(",ii,")=",amat%E%nnzrow(ii)
     !(6,*)"pilu checkpoint 44.5: F%nnzrow(",ii,")=",F%nnzrow(ii)
     !write(6,*)"pilu checkpoint 44.5: E%nnzrow(",ii,")=",E%nnzrow(ii)
     
     tnorm = 0.0 
     tmax  = 0.0
     do j=1,lenu
        tabs = L1Norm2(w2,j,nsubrow)
        if (tmax < tabs) tmax = tabs
        tnorm = tnorm + tabs
     end do
     tnorm=tnorm/lenu
     
     
     !if(tnorm>10.0) write(6,*)"pilu checkpoint 45: Schur tnorm(",ii,")=",tnorm
     
     if(abs(tnorm-DBL_EPSILON) <= DBL_EPSILON*abs(tnorm))then
        
        len = 1!+ofst!ceb
        !w(0) = 1.0!set submatrix diagonal to identity 
        jw(1) = ii
        do jjsub=1,nsubrow
           w(jjsub,jjsub,1)=1.0
        end do
        
        
        !if(jw(1).eq.ii)  write(6,*)"checkpt 45.5 jw(",1,")=",jw(1)
     else 
        
        len = 0
        !/*     tabs = drop4*tmax*(tmax/tnorm); */
        tabs = drop4*tmax*tmax/(tnorm*lenu)
        
        do j=1,lenu
           if (L1Norm2(w2,j,nsubrow) > tabs) then!submatrix norm
              len=len+1
              jpos=len
              jw(jpos) = jw2(j)
              !w(len) = w2(j)!copy submatrix
              do iisub=1,nsubrow
                 do jjsub=1,nsubcol
                    w(iisub,jjsub,jpos)= w2(iisub,jjsub,j)
                 end do
              end do
           end if
        end do
        
     end if
     
     lenu=len
     if(lenu>fil_S) lenu=fil_S !lenu = len > fil4 ? fil4 : len;
     !if(lenu>fil_S+C%nnzrow(ii)) lenu=fil_S+C%nnzrow(ii) !lenu = len > fil4 ? fil4 : len;
     schur%nnzrow(ii) = lenu
     jpos = lenu
     !if(lenu.eq.0) write(6,*)"pilu checkpoint 46.5: schur lenu(",ii,")=",lenu
     
     !if (jpos < len+ofst) call qsplit(w,jw,len+ofst,jpos,nsubrow)
     if (jpos < len) then
        !write(6,*)"pilu qsplit checkpt 4"
        call qsplit(w,jw,len,jpos,nsubrow)
     end if
     if(lenu>0) then
        !write(6,*)"allocate schur(",ii,",",lenu,")"
        allocate(schur%pa(ii)%submat(nsubrow,nsubcol,lenu),STAT=istat)
        allocate(schur%pj(ii)%cols(lenu),STAT=istat)
     end if
     
     
     !/*---------------------------------------------------------------------
     !|     copy ---
     !|--------------------------------------------------------------------*/
     !copy elements
     do k=1,jpos 
        schur%pj(ii)%cols(k)=jw(k)!copy jw into new schur_pj elements
        !if(ii.eq.20) 
        !write(6,*) "pj(",ii,k,")=",schur%pj(ii)%cols(k)
        schur%pd(ii)=-1
        if(ii.eq.jw(k)) schur%pd(ii)=k
        do iisub=1,nsubrow
           do jjsub=1,nsubrow
              schur%pa(ii)%submat(iisub,jjsub,k)=w(iisub,jjsub,k)!copy w into new schur elements, schur[k]=w[ii+(k-1)]?
           end do
        end do
     end do
  end do !end loop over rows
  
  
  
  !/*---------------------------------------------------------------------
  !|     end main loop - now do cleanup
  !|--------------------------------------------------------------------*/
  !write(6,*)"pilu checkpoint 49.5 end of main loop"
  
  
  
  nullify(lrowj)
  nullify(lrowm)
  nullify(rrowj)
  nullify(rrowm)

  deallocate(diag,STAT = istat)
  deallocate(s,STAT = istat)
  deallocate(t,STAT = istat)
  deallocate(fact,STAT = istat)
  deallocate(Lij,STAT = istat)
  deallocate(Uij,STAT = istat)
  deallocate(jw,STAT = istat)
  deallocate(w,STAT = istat)
  deallocate(jwrev,STAT = istat)
  deallocate(jw2,STAT = istat)
  deallocate(w2,STAT = istat)
  deallocate(jwrev2,STAT = istat)
  !write(6,*)"pilu checkpoint 50"
  do i=1,lsize
     if (lflen(i) > 0) then
        deallocate(lfma(i)%submat,STAT = istat)
        deallocate(lfja(i)%cols,STAT = istat)
     end if
  end do
  !write(6,*)"pilu checkpoint 51"
  
  deallocate(lfma,STAT = istat)!or just nullify here?  ceb
  deallocate(lfja,STAT = istat)
  deallocate(lflen,STAT = istat)
  !write(6,*)"pilu checkpoint 52"
  


  !do i=1,lsize
  !      write(96,*) i,"L:",(amat%L%pj(i)%cols(j),j=1,amat%L%nnzrow(i)),"U:",(amat%U%pj(i)%cols(j),j=1,amat%U%nnzrow(i))
  !      !write(76,*) i,"U:",(amat%U%pj(i)%cols(j),j=1,amat%U%nnzrow(i))
  !end do

  !  do i=1,lsize
  !     do j=1,amat%L%nnzrow(i)
  !        !write(103,*) i, amat%L%pj(i)%cols(j)
  !        write(103,*) i, amat%L%pj(i)%cols(j),((amat%L%pa(i)%submat(iisub,jjsub,j),iisub=1,amat%L%nsubrow),jjsub=1,amat%L%nsubcol)
  !     end do
  !     do j=1,amat%U%nnzrow(i)
  !        !write(104,*) i, amat%U%pj(i)%cols(j)
  !        write(104,*) i, amat%U%pj(i)%cols(j),((amat%U%pa(i)%submat(iisub,jjsub,j),iisub=1,amat%U%nsubrow),jjsub=1,amat%U%nsubcol)
  !     end do
  !  end do
  
  
  !  do i=1,schur%n
  !     do j=1,schur%nnzrow(i)
  !        !write(105,*) i, schur%pj(i)%cols(j)
  !        write(105,*) i, schur%pj(i)%cols(j),((schur%pa(i)%submat(iisub,jjsub,j),iisub=1,schur%nsubrow),jjsub=1,schur%nsubcol)
  !     end do
  !  end do
  
  
  
  !/*---------------------------------------------------------------------
  !|     done  --  correct return
  !|--------------------------------------------------------------------*/
  return! 0;
  
9991 continue
  !/*  Incomprehensible error. Matrix must be wrong.  */
  write(6,*) "Incomprehensible error. Matrix must be wrong."
  return !1;
  !/* 9992:
  !   Memory allocation error.  
  !   return 2; */
9995 continue
  !/*  illegal value for lfil or last entered  */
  write(6,*) "illegal value for lfil or last entered"
  return! 5;
9996 continue
  !/*  zero row encountered  */
  write(6,*) "Error zero row encountered index= ",ii
  return! 6;
end subroutine pilu
!/*---------------------------------------------------------------------
!|     end of pilut
!|--------------------------------------------------------------------*/

end module pilu_lib
