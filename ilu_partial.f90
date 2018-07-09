module pilu
  use kinddefs!, only: dp,ja_type,submat_type
  use blas
  implicit none
contains
!#define DBL_EPSILON 2.2204460492503131e-16 // double epsilon
type(dp) PARAMETER DBL_EPSILON 2.2204460492503131e-16 // double epsilon

!int pilu(p4ptr amat, csptr B, csptr C, double *droptol, int *lfil, csptr schur) 
subroutine pilu(amat, B, C, droptol, lfil,  schur) 
  !{
  !/*---------------------------------------------------------------------- 
  !| PARTIAL ILUT -
  !| Converted to C so that dynamic memory allocation may be implememted
  !| in order to have no dropping in block LU factors.
  !|----------------------------------------------------------------------
  !| Partial block ILU factorization with dual truncation. 
  !|                                                                      
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
  !| lfil[0]  =  number nonzeros in L-part
  !| lfil[1]  =  number nonzeros in U-part
  !| lfil[2]  =  number nonzeros in L^{-1} F
  !| lfil[3]  =  not used
  !| lfil[4]  =  number nonzeros in Schur complement
  !|
  !| droptol[0] = threshold for dropping small terms in L during
  !|              factorization.
  !| droptol[1] = threshold for dropping small terms in U.
  !| droptol[2] = threshold for dropping small terms in L^{-1} F during
  !|              factorization.
  !| droptol[3] = threshold for dropping small terms in E U^{-1} during
  !|              factorization.
  !| droptol[4] = threshold for dropping small terms in Schur complement
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
  type(p4ptr),intent(inout) :: amat
  type(csptr),intent(inout) :: B
  type(csptr),intent(inout) :: C
  type(csptr),intent(inout) :: schur
  type(dp),intent(in) :: droptol(5)
  integer,intent(in) :: lfil(5)

  integer ofst,iisub,jjsub,kkusb
  integer:: i, ii, j, jj, jcol, jpos, jrow, k, len, len2, lenu, lenl
  integer : rmax, lsize, rsize,lrowz,rrowz
  integer :: istat,nsubrow,nsubcol,fil0,fil_1,fil2,fil4
  type(dp) :: tmp, tnorm, tabs, tmax, drop0,drop1,drop2,drop3,drop4
  integer, dimension(:),allocatable :: jw, jwrev
  integer, dimension(:),allocatable :: jw2, jwrev2,lflen
  type(ja_type),dimension(:),allocatable :: lfja
  type(submat_type),dimension(:),allocatable :: lfma
  type(submat_type):: w,w2
  type(ja_type), pointer :: lrowj, rrowj !pointers to ja type
  type(submat_type), pointer :: lrowm, rrowm !pointers to submat type
  type(dp) dimension(:),allocatable :: d
  type(dp) dimension(:,:),allocatable:: Uij,Lij,s,t,fact

  !/*-----------------------------------------------------------------------*/
  ofst=1!c index offset for fortran
  fil0=lfil(0)
  fil_1=lfil(1)
  fil2=lfil(2)
  fil4=lfil(4)
 
  drop0=droptol(0)
  drop1=droptol(1)
  drop2=droptol(2)
  drop3=droptol(3)
  drop4=droptol(4)

  nsubrow=amat%nsubrow
  nsubcol=amat%nsubcol

  !temporary vars
  allocate(d(nsubrow),STAT=istat)
  allocate(s(nsubrow,nsubcol),STAT=istat)
  allocate(t(nsubrow,nsubcol),STAT=istat)
  allocate(fact(nsubrow,nsubcol),STAT=istat)
  allocate(Lij(nsubrow,nsubcol),STAT=istat)
  allocate(Uij(nsubrow,nsubcol),STAT=istat)

  lsize = amat%nB
  rsize = C%n
  rmax = rsize
  if (lsize>rsize) rmax=lsize
  
  allocate(jw(rmax),STAT=istat) !jw = (int *) Malloc(rmax*sizeof(int), "pilu:1" );
  allocate(w(rmax),STAT=istat) !w = (double *) Malloc(rmax*sizeof(double), "pilu:2" );
  allocate(jwrev(rmax),STAT=istat)!jwrev = (int *) Malloc(rmax*sizeof(int), "pilu:3" );
  allocate(jw2(rmax),STAT=istat)  !jw2 = (int *) Malloc(rmax*sizeof(int), "pilu:4" );
  allocate(w2(rmax),STAT=istat)  !w2 = (double *) Malloc(rmax*sizeof(double), "pilu:5" );
  allocate(jwrev2(rmax),STAT=istat)!jwrev2 = (int *) Malloc(rmax*sizeof(int), "pilu:6" );
  if (fil0 < 0 .or. fil_1<0 .or. amat%L%n<=0) goto 9995
  allocate(lfma(lsize),STAT=istat)!lfma = (double **) Malloc(lsize*sizeof(double *), "pilu:7" ); 
  allocate(lfja(lsize),STAT=istat)!lfja = (int **) Malloc(lsize*sizeof(int *), "pilu:8" ); 
  allocate(lflen(rmax),STAT=istat)!lflen = (int *) Malloc(lsize*sizeof(int), "pilu:9" ); 
  !/*---------------------------------------------------------------------
  !|    beginning of first main loop - L, U, L^{-1}F calculations
  !|--------------------------------------------------------------------*/
  do j=0+ofst,rmax-1+ofst
     jwrev(j) = -1
  end do
  do j=0+ofst,rmax-1+ofst
     jwrev2(j) = -1
  end do
  do ii=0+ofst,lsize-1+ofst 
     lrowj = B%pj(ii)!lrowj is pointer
     lrowm = B%pa(ii)!lrowm is pointer
     lrowz = B%nnzrow(ii)
     rrowj = amat%F%pj(ii)!rrowj is pointer
     rrowm = amat%F%pa(ii)!rrowm is pointer
     rrowz = amat%F%nnzrow(ii)
     !/*---------------------------------------------------------------------
     !|   check for zero row in B block
     !|--------------------------------------------------------------------*/
     do k=0+ofst,lrowz-1+ofst
        !if(fabs(lrowm(k)-DBL_EPSILON) > DBL_EPSILON*fabs(lrowm(k))) goto 41;
        if(abs(L1Norm2(lrowm(k))-DBL_EPSILON) > DBL_EPSILON*L1Norm(lrowm(k))) goto 41
     end do
     goto 9996 
     !/*---------------------------------------------------------------------
     !|     unpack B-block in arrays w, jw, jwrev
     !|     WE ASSUME THERE IS A DIAGONAL ELEMENT
     !|--------------------------------------------------------------------*/
41   :
     lenu = 1
     lenl = 0!lenl does not include diagonal?
     !w(ii) = 0.0!zero submatrix
     call blank_3D_submat(w,ii,nsubrow,nsubcol,rmax)
     jw(ii) = ii
     jwrev(ii) = ii
     do j=0+ofst, lrowz-1+ofst
        jcol = lrowj(j)
        !t = lrowm(j)
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              t(iisub,jjsub)=lrowm(iisub,jjsub,j)
           end do
        end do

        if (jcol < ii) then
           jpos=lenl+ofst!ceb
           !jw(lenl) = jcol !ceb  what is lenl here? must be >0
           jw(jpos) = jcol !ceb  what is lenl here? must be >0
           !w(lenl) = t!copy submatrix
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 !w(iisub,jjsub,lenl)=t(iisub,jjsub)
                 w(iisub,jjsub,jpos)=t(iisub,jjsub)
              end do
           end do
           !jwrev(jcol) = lenl;
           jwrev(jcol) = jpos!ceb
           lenl=lenl+1
        else if (jcol == ii) then
           jpos=ii!ceb
           !w(ii) = t!copy submatrix
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 w(iisub,jjsub,jpos)=t(iisub,jjsub)
              end do
           end do
        else 
           jpos = ii+lenu;
           jw(jpos) = jcol;
           !w(jpos) = t!copy submatrix
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 w(iisub,jjsub,jpos)=t(iisub,jjsub)
              end do
           end do
           jwrev(jcol) = jpos
           lenu=lenu+1
        end if
     end do
     !/*---------------------------------------------------------------------
     !|     unpack F-block in arrays w2, jw2, jwrev2 
     !|     (all entries are in U portion)
     !|--------------------------------------------------------------------*/
     len2 = 0
     do j=0+ofst,rrowz-1+ofst
        jcol = rrowj(j)
        jw2(len2+ofst) = jcol
        !w2(len2) = rrowm(j)!copy submatrix
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              w(iisub,jjsub,len2+ofst)=rrowm(iisub,jjsub,j)
           end do
        end do
        jwrev2(jcol+ofst) = len2+ofst!ceb
        len2=len2+1
     end do
     !/*---------------------------------------------------------------------
     !|     Eliminate previous rows -  
     !|--------------------------------------------------------------------*/
     len = 0
     jj=0!ceb
     do !jj=0+ofst,lenl-1+ofst
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
        do j=jj+1,lenl-1+ofst
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
           !/*   exchange in w   */
           !s = w(jj)!copy submatrix
           !w(jj) = w(k)!copy submatrix
           !w(k) = s!copy submatrix
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 s(iisub,jjsub)=w(iisub,jjsub,jj)!s = w[jj];
                 w(iisub,jjsub,jj)=w(iisub,jjsub,k)!w[jj] = w[k];
                 w(iisub,jjsub,k)=s(iisub,jjsub)!w[k] = s;
              end do
           end do
        end if
        !/*---------------------------------------------------------------------
        !|     zero out element in row.
        !|--------------------------------------------------------------------*/
        jwrev(jrow) = -1
        !/*---------------------------------------------------------------------
        !|     get the multiplier for row to be eliminated (jrow).
        !|--------------------------------------------------------------------*/
        lrowm = amat%U%pa(jrow)
        !fact = w(jj) * lrowm(0)!submatrix mult
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              fact(iisub,jjsub)=0.
              do kksub=1,nsubrow !assumes square submatrix
                 fact(iisub,jjsub) = fact(iisub,jjsub) + w(iisub,kksub,jj)*lrowm%submat(kksub,jjsub,1)
              end do
           end do
        end do
        !if (fabs(fact) > drop0 ) then   !/*   DROPPING IN L   */
        if (L1Norm(fact,nsubrow) > drop0 ) then   !/*   DROPPING IN L   */
           lrowj = amat->U->pj(jrow)
           lrowz = amat->U->nnzrow(jrow)
           rrowj = lfja(jrow)
           rrowm = lfma(jrow)
           rrowz = lflen(jrow)
           !/*---------------------------------------------------------------------
           !|     combine current row and row jrow
           !|--------------------------------------------------------------------*/
           do k=1+ofst,lrowz-1+ofst
              !s = fact * lrowm(k)!submatrix mult
              do iisub=1,nsubrow
                 do jjsub=1,nsubcol
                    s(iisub,jjsub)=0.
                    do kksub=1,nsubrow!assumes square submatrix
                       s(iisub,jjsub)=s(iisub,jjsub)+fact(iisub,kksub)*lrowm%submat(kksub,jjsub,k)
                    end do
                 end do
              end do
              j = lrowj(k)
              jpos = jwrev(j)
              !/*---------------------------------------------------------------------
              !|     dealing with U
              !|--------------------------------------------------------------------*/
              if (j >= ii) then
                 !/*---------------------------------------------------------------------
                 !|     this is a fill-in element
                 !|--------------------------------------------------------------------*/
                 if (jpos == -1) then
                    if (lenu > lsize) then !printf("U  row = %d\n",ii
                       goto 9991
                    endif
                    i = ii + lenu
                    jw(i) = j
                    jwrev(j) = i
                    !w(i) = - s!copy submatrix
                    do iisub=1,nsubrow
                       do jjsub=1,nsubcol
                          w(iisub,jjsub,i)=-s(iisub,jjsub)
                       end do
                    end do
                    lenu=lenu+1
                    !/*---------------------------------------------------------------------
                    !|     this is not a fill-in element 
                    !|--------------------------------------------------------------------*/
                 else
                    !w(jpos) -= s;!copy submatrix
                    do iisub=1,nsubrow
                       do jjsub=1,nsubcol
                          w(iisub,jjsub,jpos)= w(iisub,jjsub,jpos)-s(iisub,jjsub)
                       end do
                    end do
                 end if
                 !/*---------------------------------------------------------------------
                 !|     dealing  with L
                 !|--------------------------------------------------------------------*/
              else 
                 !/*---------------------------------------------------------------------
                 !|     this is a fill-in element
                 !|--------------------------------------------------------------------*/
                 if (jpos == -1) then
                    if (lenl > lsize) then!{printf("L  row = %d\n",ii);
                       goto 9991
                    end if
                    jw(lenl) = j;
                    jwrev(j) = lenl;
                    !w(lenl) = - s;!copy submatrix
                    do iisub=1,nsubrow
                       do jjsub=1,nsubcol
                          w(iisub,jjsub,lenl)= -s(iisub,jjsub)
                       end do
                    end do
                    lenl=lenl+1;
                    !/*---------------------------------------------------------------------
                    !|     this is not a fill-in element 
                    !|--------------------------------------------------------------------*/
                 else
                    !w(jpos) -= s!subtract submatrix
                    do iisub=1,nsubrow
                       do jjsub=1,nsubcol
                          w(iisub,jjsub,jpos)= w(iisub,jjsub,jpos)-s(iisub,jjsub)
                       end do
                    end do
                 end if
              end if
           end do
           !/*---------------------------------------------------------------------
           !|     dealing  with  L^{-1} F
           !|--------------------------------------------------------------------*/
           do k=0+ofst,rrowz-1+ofst 
              !s = fact * rrowm(k)!submatrix multiply
              do iisub=1,nsubrow
                 do jjsub=1,nsubcol
                    s(iisub,jjsub)=0.
                    do kksub=1,nsubrow !assumes square submatrix
                       s(iisub,jjsub) = s(iisub,jjsub) + fact(iisub,kksub)*rrowm%submat(kksub,jjsub,k)
                    end do
                 end do
              end do
              j = rrowj(k)
              jpos = jwrev2(j)
              !/*---------------------------------------------------------------------
              !|     this is a fill-in element
              !|--------------------------------------------------------------------*/
              if (jpos == -1) then
                 jw2(len2) = j
                 jwrev2(j) = len2
                 !w2(len2) = - s!subtract submatrix
                 do iisub=1,nsubrow
                    do jjsub=1,nsubcol
                       w2(iisub,jjsub,len2)= -s(iisub,jjsub)
                    end do
                 end do
                 len2=len2+1
                 !/*---------------------------------------------------------------------
                 !|     this is not a fill-in element
                 !|--------------------------------------------------------------------*/
              else
                 !w2(jpos) -= s!subtract submatrix
                 do iisub=1,nsubrow
                    do jjsub=1,nsubcol
                       w2(iisub,jjsub,jpos)= w2(iisub,jjsub,jpos)-s(iisub,jjsub)
                    end do
                 end do
              end if
           end do
           !/*---------------------------------------------------------------------
           !|     store this pivot element
           !|--------------------------------------------------------------------*/
           !w(len+ofst) = fact!copy submatrix
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 w(iisub,jjsub,len+ofst)= fact(iisub,jjsub)
              end do
           end do
           jw(len+ofst)  = jrow!ceb
           len=len+1
        end if
     end do!end loop over block cols: 1..lenl
     !/*---------------------------------------------------------------------
     !|     reset nonzero indicators
     !|--------------------------------------------------------------------*/
     do j=0+ofst,len2+ofst    !/*  L^{-1} F block  */
        jwrev2(jw2(j)) = -1
     end do
     do j=0+ofst,lenl-1+ofst    !/*  L block  */
        jwrev(jw(j)) = -1
     end do
     do j=0+ofst,lenu-1+ofst    !/*  U block  */
        jwrev(jw(ii+j)) = -1
     end do
     !/*---------------------------------------------------------------------
     !|     done reducing this row, now store L
     !|--------------------------------------------------------------------*/
     lenl=len
     if (len > fil0) lenl=fil0
     !lenl = len > fil0 ? fil0 : len;
     amat%L%nnzrow(ii) = lenl

     !if (lenl < len) call qsplitC(w, jw, len, lenl)
     itmp=1
     if (len > lenl) call qsplitC(w,jw,len+ofst,lenl+ofst,nsubrow,itmp)!quick sort alg (sorts based on L1norm of submats)
     !ceb check above

     if (len > 0) then
        !amat->L->pj(ii) = (int *) Malloc(lenl*sizeof(int), "pilu:10" ); 
        allocate(amat%L%pj(ii)%cols(lenl),STAT=istat)
        !amat->L->pa(ii) = (double *) Malloc(lenl*sizeof(double), "pilu:11" );
        allocate(amat%L%pa(ii)%submat(nsubrow,nsubcol,lenl),STAT=istat)
        !copy entries
        !memcpy(amat->L->pj(ii), jw, lenl*sizeof(int));
        !memcpy(amat->L->pa(ii), w, lenl*sizeof(double));
        do k=1,lenl
           L%pj(ii)%cols(k)=jw(k)!copy jw into new U_ja elements
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 L%pa(ii)%submat(iisub,jjsub,k)=w(iisub,jjsub,k))!copy w into new U elements, U[k]=w[ii+(k-1)]
              end do
           end do
        end do
     end if

     !/*---------------------------------------------------------------------
     !|     store the diagonal element of U
     !|     dropping in U if size is less than drop1 * diagonal entry
     !|--------------------------------------------------------------------*/

     !t = w(ii)!copy submatrix
     do iisub=1,nsubrow
        do jjsub=1,nsubcol
           t(iisub,jjsub)=w(iisub,jjsub,ii)
        end do
     end do
     !tnorm = abs(t)!get matrix norm
     tnorm= L1Norm(t,nsubrow)
     len = 0
     do j=1+ofst,lenu-1+ofst
        !if ( abs(w(ii+j)) > drop1*tnorm ) then !use matrix norm
        if ( L1Norm(w(ii+j),nsubrow) > drop1*tnorm ) then !use matrix norm
           !w(len) = w(ii+j)!copy submatrix
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 w(iisub,jjsub,len+ofst)= w(iisub,jjsub,ii+j)
              end do
           end do
           jw(len+ofst) = jw(ii+j)
           len=len+1
        end if
     end do

     lenu=len+1
     if (len+1 > fil_1) lenu=fil_1
     !lenu = len+1 > fil_1 ? fil_1 : len+1;

     amat%U%nnzrow(ii) = lenu
     jpos = (lenu-1)+ofst !jpos = lenu - 1
     !if (jpos < len) call qsplitC(w, jw, len, jpos)
     itmp=ii+1
     if (len > (lenu-1)) call qsplitC(w,jw,len+ofst,jpos,neqn,itmp)!quick sort alg (sorts based on L1norm of submats)
     !ceb check above

     !amat->U->pa(ii) = (double *) Malloc(lenu*sizeof(double), "pilu:12" );
     allocate(amat%U%pa(ii)%submat(nsubrow,nsubcol,lenu),STAT=istat)
     !amat->U->pj(ii) = (int *) Malloc(lenu*sizeof(int), "pilu:13" );
     allocate(amat%U%pj(ii)%cols(lenu),STAT=istat)


     !if(abs(L1Norm(t,nsubrow)-DBL_EPSILON) <= DBL_EPSILON*L1Norm(t,nsubcol))then
     !   !t=(0.0001+drop1);
     !   do iisub,nsubrow
     !      !do jjsub,nsubcol
     !      t(iisub,iisub)=(0.0001+drop1)
     !      !end do
     !   end do
     !end if
     !amat%U%pa(ii)%submat(,,0) = 1.0 / t;!this is a matrix invers
     !amat%U%pa(ii)%submat(iisub,jjsub,0) = inv_t()!this is a matrix invers

     !amat%U%pj(ii)%cols(0+ofst) = ii !ceb moved to follow diag formation

     !/*---------------------------------------------------------------------
     !|     copy the rest of U
     !|--------------------------------------------------------------------*/
     !memcpy(&amat->U->pj(ii)(1), jw, jpos*sizeof(int));
     !memcpy(&amat->U->pa(ii)(1), w, jpos*sizeof(double));
     !/* copy the rest of U */   !ceb do we need to deal with diagonal here? 
     do k=2,jpos! (skips diagonal)
        U%pj(ii)%cols(k)=jw(k)!copy jw into new U_pj elements
        do iisub=1,nsubrow
           do jjsub=1,nsubrow
              U%pa(ii)%submat(iisub,jjsub,k)=w(iisub,jjsub,k)!copy w into new U elements, U[k]=w[ii+(k-1)]
           end do
        end do
     end do


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
              Lij(kksub,kksub)= 1.0/tmp 
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
           U%pa(ii)%submat(jjsub,kksub,1)=0.0
        end do
        U%pa(ii)%submat(jjsub,jjsub,1)=1.0
     end do
     
     !//write LU into U[diagp]
     do kksub=1,nsubcol !loop subcols
        !{  
        do jjsub=1,nsubcol-1 !get D[1,..,neqn-1] ! loop subcols
           d(jjsub)=U%pa(ii)%submat(jjsub,kksub,1) !Identity matrix row jj 
           do iisub=1,jjsub-1
              d(jjsub) = d(jjsub) - Lij(iisub,jjsub)*d(iisub)
           end do
           d(jjsub) = d(jjsub)*Lij(jjsub,jjsub)!ceb not sure if this is right
        end do
        
        do jjsub=nsubcol,1,-1 !get alu[diagp,jj,kk] ! loop subcols
           if(jjsub.eq.nsubcol) then     
              !D=D-L*D (left of subdiag)
              do iisub=1,jjsub-1
                 U%pa(ii)%submat(jjsub,kksub,1) = U%pa(ii)%submat(jjsub,kksub,1)  - Lij(iisub,jjsub)*d(iisub)
              end do
              U%pa(ii)%submat(jjsub,kksub,1) = U%pa(ii)%submat(jjsub,kksub,1)*Lij(jjsub,jjsub)
              !D=
           else      
              !D=D-U*D (right of sub diag)
              U%pa(ii)%submat(jjsub,kksub,1)=d(jjsub)
              do iisub=jjsub+1,nsubrow
        U%pa(ii)%submat(jjsub,kksub,1) = U%pa(ii)%submat(jjsub,kksub,1) - Uij(iisub,jjsub)*U%pa(ii)%submat(iisub,kksub,1)
              end do
           end if
        end do
        !}
     end do
   
     amat%U%pj(ii)%cols(0+ofst) = ii!ceb moved from prior to diag formation
  
     !U%pj(ii)%cols(0+ofst) = jw(ii) !rowj(0) = jw(ii);
     U%pj(ii)%cols(0+ofst) = ii !rowj(0) = jw(ii);
     
     !===========================================================
!--ceb





     !/*---------------------------------------------------------------------
     !|     copy  L^{-1} F
     !|--------------------------------------------------------------------*/
     len = 0;
     do j=0+ofst,len2-1+ofst
        if ( L1Norm(w2,j,nsubrow) > drop2*tnorm ) then !use matrix norm
           !w(len) = w2(j);!copy submatrix
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 w(iisub,jjsub,len)= w2(iisub,jjsub,j)
              end do
           end do
           jw(len) = jw2(j)
           len=len+1
        end if
     end do

     len=fil2
     if(len>fil2) lenu=len
     !lenu = len > fil2 ? fil2 : len;

     !if (lenu < len) call qsplitC(w, jw, len, lenu);
     itmp=ii+1
     if ((lenu-1) < len) call qsplitC(w,jw,len+ofst,lenu,nsubrow,itmp)!quick sort alg (sorts based on L1norm of submats)

     lflen(ii) = lenu
     
     if (lenu > 0) then
        !lfja(ii)  = (int *) Malloc(lenu*sizeof(int), "pilu:14" );
        allocate(lfja(ii)%cols(lenu),STAT=istat)
        !lfma(ii)  = (double *) Malloc(lenu*sizeof(double), "pilu:15" ); 
        allocate(lfma(ii)%submat(nsubrow,nsubcol,lenu),STAT=istat)
        !copy elements
        !memcpy(lfma(ii), w, lenu*sizeof(double));
        !memcpy(lfja(ii), jw, lenu*sizeof(int)); 
     end if
  end do

  !/*---------------------------------------------------------------------
  !|    beginning of second main loop   E U^{-1} and Schur complement
  !|--------------------------------------------------------------------*/
  do ii=0+ofst,rsize-1+ofst
     lrowj = amat%E%pj(ii)
     lrowm = amat%E%pa(ii)
     lrowz = amat%E%nnzrow(ii)
     rrowj = C%pj(ii)
     rrowm = C%pa(ii)
     rrowz = C%nnzrow(ii)
     
     !/*---------------------------------------------------------------------
     !|    determine if there is a zero row in [ E C ]
     !|--------------------------------------------------------------------
     !     for (k=0; k<lrowz; k++)
     !       if (lrowm(k) != 0.0) goto 42;
     !     for (k=0; k<rrowz; k++)
     !       if (rrowm(k) != 0.0) goto 42;
     !     goto 9997;
     !     42:
     !*/
     
     !/*---------------------------------------------------------------------
     !|     unpack E in arrays w, jw, jwrev
     !|--------------------------------------------------------------------*/
     lenl = 0
     do j=0+ofst,lrowz-1+ofst
        jcol = lrowj(j)
        jw(lenl) = jcol
        !w(lenl) = lrowm(j)!copy submatrix
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              w(iisub,jjsub,lenl)= lrowm(iisub,jjsub,j)
           end do
        end do
        jwrev(jcol) = lenl
        lenl=lenl+1
     end do
     !/*---------------------------------------------------------------------
     !|     unpack C in arrays w2, jw2, jwrev2    
     !|--------------------------------------------------------------------*/
     lenu = 0
     do j=0+ofst,rrowz-1+ofst
        jcol = rrowj(j)
        jw2(lenu) = jcol
        !w2(lenu) = rrowm(j)!copy submatrix
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              w2(iisub,jjsub,lenu)= rrowm(iisub,jjsub,j)
           end do
        end do
        jwrev2(jcol) = lenu
        lenu=lenu+1
     end do
     !/*---------------------------------------------------------------------
     !|     eliminate previous rows
     !|--------------------------------------------------------------------*/
     len = 0
     do jj=0+ofst,lenl-1+ofst
        !/*---------------------------------------------------------------------
        !|    in order to do the elimination in the correct order we must select
        !|    the smallest column index among jw(k), k=jj+1, ..., lenl.
        !|--------------------------------------------------------------------*/
        jrow = jw(jj)
        k = jj
        !/*---------------------------------------------------------------------
        !|     determine smallest column index
        !|--------------------------------------------------------------------*/
        do j=jj+1,lenl-1+ofst
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
           !/*   exchange in w   */
           !s = w(jj)!copy submatrix
           !w(jj) = w(k)!copy submatrix
           !w(k) = s!copy submatrix
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 s(iisub,jjsub)=w(iisub,jjsub,jj)!s = w[jj];
                 w(iisub,jjsub,jj)=w(iisub,jjsub,k)!w[jj] = w[k];
                 w(iisub,jjsub,k)=s(iisub,jjsub)!w[k] = s;
              end do
           end do
        end if
        !/*---------------------------------------------------------------------
        !|     zero out element in row.
        !|--------------------------------------------------------------------*/
        jwrev(jrow) = -1
        !/*---------------------------------------------------------------------
        !|     get the multiplier for row to be eliminated (jrow).
        !|--------------------------------------------------------------------*/
        lrowm = amat%U%pa(jrow)
        !fact = w(jj) * lrowm(0)!submatrix multiply
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              fact(iisub,jjsub)=0.
              do kksub=1,nsubrow !assumes square submatrix
                 fact(iisub,jjsub) = fact(iisub,jjsub) + w(iisub,kksub,jj)*lrowm%submat(kksub,jjsub,1)
              end do
           end do
        end do
        if ( abs(fact) > drop3 ) then!use matrix norm      !/*  DROPPING IN E U^{-1}   */
           lrowj = amat%U%pj(jrow)
           lrowz = amat%U%nnzrow(jrow)
           rrowj = lfja(jrow)
           rrowm = lfma(jrow)
           rrowz = lflen(jrow)
           !/*---------------------------------------------------------------------
           !|     combine current row and row jrow   -   first  E U^{-1}
           !|--------------------------------------------------------------------*/
           do k=1+ofst,lrowz-1+ofst 
              !s = fact * lrowm(k)!submatrix multiply
              do iisub=1,nsubrow
                 do jjsub=1,nsubcol
                    s(iisub,jjsub)=0.
                    do kksub=1,nsubrow !assumes square submatrix
                       s(iisub,jjsub) = s(iisub,jjsub) + fact(iisub,kksub)*lrowm%submat(kksub,jjsub,k)
                    end do
                 end do
              end do
              j = lrowj(k)
              jpos = jwrev(j)
              !/*---------------------------------------------------------------------
              !|     fill-in element
              !|--------------------------------------------------------------------*/
              if (jpos .eq. -1) then
                 if (lenl > lsize) then! {printf(" E U^{-1}  row = %d\n",ii);
                    goto 9991
                 end if
                 jw(lenl) = j
                 jwrev(j) = lenl
                 !w(lenl) = - s!subtract submatrix
                 do iisub=1,nsubrow
                    do jjsub=1,nsubcol
                       w(iisub,jjsub,lenl)= s(iisub,jjsub)
                    end do
                 end do
                 lenl=lenl+1
                 !/*---------------------------------------------------------------------
                 !|     this is not a fill-in element 
                 !|--------------------------------------------------------------------*/
              else
                 !w(jpos) -= s!subtract submatrix
                 do iisub=1,nsubrow
                    do jjsub=1,nsubcol
                       w(iisub,jjsub,jpos)= w(iisub,jjsub,jpos)-s(iisub,jjsub)
                    end do
                 end do
              end if
           end do
           !/*---------------------------------------------------------------------
           !|     incorporate into Schur complement   C - (E U^{-1}) (L^{-1} F)
           !|--------------------------------------------------------------------*/
           do k=0+ofst,rrowz-1+ofst
              !s = fact * rrowm(k)!submatrix multiply
              do iisub=1,nsubrow
                 do jjsub=1,nsubcol
                    s(iisub,jjsub)=0.
                    do kksub=1,nsubrow !assumes square submatrix
                       s(iisub,jjsub) = s(iisub,jjsub) + fact(iisub,kksub)*rrowm%submat(kksub,jjsub,k)
                    end do
                 end do
              end do
              j = rrowj(k)
              jpos = jwrev2(j)
              !/*---------------------------------------------------------------------
              !|     this is not a fill-in element 
              !|--------------------------------------------------------------------*/
              if (jpos == -1) then
                 jw2(lenu) = j
                 jwrev2(j) = lenu
                 !w2(lenu) = - s!copy submatrix
                 do iisub=1,nsubrow
                    do jjsub=1,nsubcol
                       w2(iisub,jjsub,lenu)= -s(iisub,jjsub)
                    end do
                 end do
                 lenu=lenu+1
                 !/*---------------------------------------------------------------------
                 !|     this is not a fill-in element
                 !|--------------------------------------------------------------------*/
              else
                 !w2(jpos) -= s!subtract submatrix
                 do iisub=1,nsubrow
                    do jjsub=1,nsubcol
                       w2(iisub,jjsub,jpos)= w2(iisub,jjsub,jpos)-s(iisub,jjsub)
                    end do
                 end do
              end if
           end do
           !/*---------------------------------------------------------------------
           !|     store this pivot element
           !|--------------------------------------------------------------------*/
           !w(len) = fact!copy submatrix
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 w(iisub,jjsub,len)= fact(iisub,jjsub)
              end do
           end do
           jw(len) = jrow
           len=len+1
        end if
     end do
     !/*---------------------------------------------------------------------
     !|     reset nonzero indicators
     !|--------------------------------------------------------------------*/
     do j=0+ofst,lenu-1+ofst    !/*  Schur complement  */
        jwrev2(jw2(j)) = -1
     end do
     do j=0+ofst,lenl-1+ofst    !/*  E U^{-1} block  */
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
     tnorm = 0.0 
     tmax  = 0.0
     do j=0+ofst,lenu-1+ofst 
        !tabs = fabs(w2(j)) !matrix norm
        tabs = L1Norm(w2(j),nsubrow)
        if (tmax < tabs) tmax = tabs
        tnorm = tnorm + tabs
     end do
     !/* if (fabs(w2(j)) > tnorm) tnorm =  fabs(w2(j)); */
     if(abs(tnorm-DBL_EPSILON) <= DBL_EPSILON*abs(tnorm))then
        len = 1
        !w(0) = 1.0!set submatrix diagonal to identity 
        do jjsub=1,nsubrow
           !do kksub=1,nsubcol
           !   w(jjsub,kksub,1)=0.0
           !end do
           w(jjsub,jjsub,1)=1.0
        end do

        jw(0+ofst) = ii
     else {
        len = 0
        !/*     tabs = drop4*tmax*(tmax/tnorm); */
        tabs = drop4*tmax*tmax/( tnorm * (double) lenu)
        do j=0+ofst,lenu-1+ofst
           if (fabs(w2(j)) > tabs) then!submatrix norm
              !w(len) = w2(j)!copy submatrix
              do iisub=1,nsubrow
                 do jjsub=1,nsubcol
                    w(iisub,jjsub,len)= w2(iisub,jjsub,j)
                 end do
              end do
              jw(len) = jw2(j)
              len=len+1
           end if
        end do
     end if

     lenu=fil4
     if(len>fil4) lenu=len
     !lenu = len > fil4 ? fil4 : len;
     schur%nnzrow(ii) = lenu
     jpos = lenu

     !if (jpos < len) call qsplitC(w, jw, len, jpos)

     itmp=ii+1
     if ((lenu-1) < len) call qsplitC(w,jw,len+ofst,jpos,neqn,itmp)!quick sort alg (sorts based on L1norm of submats)
     !ceb check above

     !schur->pa(ii) = (double *) Malloc(lenu*sizeof(double), "pilu:16" );
     allocate(schur%pa(ii)%submat(nsubrow,nsubcol,lenu),STAT=istat)
     !schur->pj(ii) = (int *) Malloc(lenu*sizeof(int), "pilu:17" );
     allocate(schur%pj(ii)%cols(lenu),STAT=istat)
     !/*---------------------------------------------------------------------
     !|     copy ---
     !|--------------------------------------------------------------------*/
     !copy elements
     !memcpy(&schur->pj(ii)(0), jw, jpos*sizeof(int));
     !memcpy(&schur->pa(ii)(0), w, jpos*sizeof(double));
  end do
  !/*---------------------------------------------------------------------
  !|     end main loop - now do cleanup
  !|--------------------------------------------------------------------*/
  deallocate(d,STAT = istat)
  deallocate(s,STAT = istat)
  deallocate(t,STAT = istat)
  deallocate(Lij,STAT = istat)
  deallocate(Uij,STAT = istat)
  deallocate(jw,STAT = istat)
  deallocate(w,STAT = istat)
  deallocate(jwrev,STAT = istat)
  deallocate(jw2,STAT = istat)
  deallocate(w2,STAT = istat)
  deallocate(jwrev2,STAT = istat)
  do i=0+ofst,lsize-1+ofst
     if (lflen(i) > 0) then
        deallocate(lfma(i)%submat,STAT = istat)
        deallocate(lfja(i)%cols,STAT = istat)
     end if
  end do
  deallocate(lfma,STAT = istat)
  deallocate(lfja,STAT = istat)
  deallocate(lflen,STAT = istat)
  !/*---------------------------------------------------------------------
  !|     done  --  correct return
  !|--------------------------------------------------------------------*/
  return! 0;

9991 :
  !/*  Incomprehensible error. Matrix must be wrong.  */
  write(6,*) "Incomprehensible error. Matrix must be wrong."
  return 1;
  !/* 9992:
  !   Memory allocation error.  
  !   return 2; */
9995 :
  !/*  illegal value for lfil or last entered  */
  write(6,*) "illegal value for lfil or last entered"
  return! 5;
9996 :
  !/*  zero row encountered  */
  write(6,*) "zero row encountered index= ",ii+1
  return! 6;
end subroutine pilu
!/*---------------------------------------------------------------------
!|     end of pilut
!|--------------------------------------------------------------------*/

end module pilu
