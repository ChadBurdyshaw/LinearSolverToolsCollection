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
  type(p4_type),pointer :: amat!should already be allocated
  type(cs_type),pointer :: B
  type(cs_type),pointer :: C
  type(cs_type),pointer :: schur
  real(dp),intent(in) :: droptol(7)
  integer,intent(in) :: lfil(5)
  !internal vars
  integer iisub,jjsub,kksub
  integer:: i, ii, j, jj, jcol, jpos, jrow, k, len, len2, lenu, lenl
  integer :: rmax, lsize, rsize,lrowz,rrowz,itmp
  integer :: test,istat,nsubrow,nsubcol,fil_0,fil_1,fil_2,fil_4
  real(dp) :: tmp, tnorm, tabs, tmax, drop0,drop1,drop2,drop3,drop4
  integer, target,dimension(:),allocatable :: jw, jwrev, jw2, jwrev2,lflen
  type(ja_type),target,dimension(:),allocatable :: lfja
  type(submat_type),target,dimension(:),allocatable :: lfma
  real(dp),target,dimension(:,:,:),allocatable :: w,w2
  type(ja_type), pointer :: lrowj, rrowj !pointers to ja type
  type(submat_type), pointer :: lrowm, rrowm !pointers to submat type
  real(dp), dimension(:),allocatable :: d
  real(dp), dimension(:,:),allocatable:: Uij,Lij,s,t,fact

  real(dp),pointer,dimension(:,:,:) :: shift_w
  integer,pointer,dimension(:) :: shift_jw
  !/*-----------------------------------------------------------------------*/
  fil_0 =lfil(1)
  fil_1 =lfil(2)
  fil_2 =lfil(3)
  fil_4 =lfil(4)
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
  !write(6,*)"0 pilu: rmax=",rmax


  test=0!ceb


  allocate(jw(rmax),STAT=istat) !jw = (int *) Malloc(rmax*sizeof(int), "pilu:1" );
  allocate(w(nsubrow,nsubcol,rmax),STAT=istat) !w = (double *) Malloc(rmax*sizeof(double), "pilu:2" );
  allocate(jwrev(rmax),STAT=istat)!jwrev = (int *) Malloc(rmax*sizeof(int), "pilu:3" );
  allocate(jw2(rmax),STAT=istat)  !jw2 = (int *) Malloc(rmax*sizeof(int), "pilu:4" );
  allocate(w2(nsubrow,nsubcol,rmax),STAT=istat)  !w2 = (double *) Malloc(rmax*sizeof(double), "pilu:5" );
  allocate(jwrev2(rmax),STAT=istat)!jwrev2 = (int *) Malloc(rmax*sizeof(int), "pilu:6" );
  if (fil_0 < 0 .or. fil_1<0 .or. amat%L%n<=0) goto 9995

  allocate(lfma(lsize),STAT=istat)!lfma = (double **) Malloc(lsize*sizeof(double *), "pilu:7" ); 
  allocate(lfja(lsize),STAT=istat)!lfja = (int **) Malloc(lsize*sizeof(int *), "pilu:8" ); 
  allocate(lflen(rmax),STAT=istat)!lflen = (int *) Malloc(lsize*sizeof(int), "pilu:9" ); 


!ceb--
  do i=1,B%n
     do j=1,B%nnzrow(i)
        write(95,*) i,B%pj(i)%cols(j)
        write(97,*) ((B%pa(i)%submat(iisub,jjsub,j),iisub=1,B%nsubrow),jjsub=1,B%nsubcol)
     end do
  end do
!  !write(6,*)""
!  do i=1,C%n
!     do j=1,C%nnzrow(i)
!        !write(6,*) "pilu: C%pj(",i,")%cols(",j,")=",C%pj(i)%cols(j)
!        write(67,*) i,C%pj(i)%cols(j)
!     end do
!  end do
!--ceb















if (test.eq.1) then


    if (fil_0<0 .or. fil_1<0 .or. rmax<=0) goto  9996

    do j=1,rmax 
       jwrev(j) = -1 
       !iprev(j) = j 
    end do

    !/*---------------------------------------------------------------------
    !|    beginning of main loop - L, U calculations
    !|--------------------------------------------------------------------*/
    do ii=1,rmax 
       lrowj => B%pj(ii) 
       lrowm => B%pa(ii) 
       lrowz =  B%nnzrow(ii) 
       !/*---------------------------------------------------------------------
       !|   check for zero row 
       !|--------------------------------------------------------------------*/
       tnorm = 0.0 
       do k=1,lrowz
          tnorm = tnorm + L1Norm2(lrowm%submat,k,nsubrow) !tnorm = tnorm + abs(rowm(k)) 
       end do
       tnorm = tnorm / lrowz 
!if(tnorm>10.)
! write(6,*) "1: row ",ii,"tnorm=",tnorm
       if(abs(tnorm-DBL_EPSILON) .le. DBL_EPSILON*tnorm) goto  9996

       !/*---------------------------------------------------------------------
       !|     unpack amat in arrays w, jw, jwrev
       !|     WE ASSUME THERE IS A DIAGONAL ELEMENT
       !|--------------------------------------------------------------------*/
       lenu = 1 
       lenl = 0 
       jw(ii) = ii 
       jwrev(ii) = ii 
       call blank_3D_submat(w,ii,nsubrow,nsubcol) !,rmax)!w(ii) = 0.0 
       do j=1,lrowz 
          !jcol = iprev(rowj%cols(j)) 
          jcol = lrowj%cols(j) 
          if (jcol < ii) then 
             lenl=lenl+1 
             jpos=lenl
             jw(jpos) = jcol 
             jwrev(jcol) = jpos 
          elseif (jcol .eq. ii) then
             jpos=ii
          else 
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
       !|     eliminate previous rows
       !|--------------------------------------------------------------------*/
       len = 0 
       jj=0
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
!write(6,*)"exchanging col ",k," with ",jj   
             !/*   exchange in jw   */
             j = jw(jj) 
             jw(jj) = jw(k) 
             jw(k) = j 
             !/*   exchange in jwrev   */
             jwrev(jrow) = jj 
             jwrev(j) = k 
             !/*   exchange in w   */
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
          lrowm => amat%U%pa(jrow) 
          lrowj => amat%U%pj(jrow) 
          lrowz = amat%U%nnzrow(jrow) 
          !fact = w(jj) * rowm(0) = row A(ii)*U_d(ii)
          do iisub=1,nsubrow
             do jjsub=1,nsubcol
                fact(iisub,jjsub)=0.0
                do kksub=1,nsubrow !assumes square submatrix
                   fact(iisub,jjsub) = fact(iisub,jjsub) + w(iisub,kksub,jj)*lrowm%submat(kksub,jjsub,1)
                end do
             end do
          end do

!if(L1Norm(fact,nsubrow)>10)write(6,*)"2: ii:",ii," pivot norm(",jrow,")=",L1Norm(fact,nsubrow),"drop1=",drop1

          if (L1Norm(fact,nsubrow) <= drop1 ) then  !/*  DROPPING IN L  */
             !write(6,*)"dropping in L  row,column(",ii,",",jrow,")",L1Norm(fact,nsubrow),"<=",drop1
             cycle
          end if

          !/*---------------------------------------------------------------------
          !|     combine current row and row jrow
          !|--------------------------------------------------------------------*/
          do k=2,lrowz 
             !s = fact * rowm(k)
             !jcol = iprev(rowj%cols(k))    !/*  new column number  */
             jcol = lrowj%cols(k)    !/*  new column number  */
             jpos = jwrev(jcol)        
             if(jpos .eq. -1) then 
                !/*---------------------------------------------------------------------
                !|     this is a fill-in element
                !|--------------------------------------------------------------------*/
                if (jcol .ge. ii) then
                   !/*---------------------------------------------------------------------
                   !|     dealing with U
                   !|--------------------------------------------------------------------*/
                   if(lenu < fil_1) then
                      jpos = ii + lenu
                      jw(jpos) = jcol
                      jwrev(jcol) = jpos 
                      lenu=lenu+1 
                      call blank_3D_submat(w,jpos,nsubrow,nsubcol)
                   end if                 
                else !left of diag
                   !/*---------------------------------------------------------------------
                   !|     dealing  with L
                   !|--------------------------------------------------------------------*/
                   if(lenl < fil_0) then
                      lenl=lenl+1 
                      jpos=lenl
                      jw(jpos) = jcol 
                      jwrev(jcol) = jpos 
                      call blank_3D_submat(w,jpos,nsubrow,nsubcol)
                   end if
                end if
             !else
                !/*---------------------------------------------------------------------
                !|     this is not a fill-in element 
                !|--------------------------------------------------------------------*/
             end if

             !perform row operation w(ii)-=pivot*U(k)
             if(jpos.ne.-1) then
                do iisub=1,nsubrow
                   do jjsub=1,nsubcol
                      do kksub=1,nsubrow
                         w(iisub,jjsub,jpos) = w(iisub,jjsub,jpos) - fact(iisub,kksub)*lrowm%submat(kksub,jjsub,k)
                      end do
                   end do
                end do
             end if
             
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

       end do!end loop over block cols


       !/*---------------------------------------------------------------------
       !|     reset nonzero indicators
       !|--------------------------------------------------------------------*/
       do j=1,lenu    !/*  U block  */
          jwrev(jw(ii+j-ofst)) = -1 
       end do


       !/*---------------------------------------------------------------------
       !|     done reducing this row, now store L
       !|--------------------------------------------------------------------*/       lenl = len
       len = lenl
       if(lenl>fil_0) len=fil_0 !len = lenl > fil5 ? fil5 : lenl 
       amat%L%nnzrow(ii) = len


     !/* quick sort */
       !if (lenl > len) call qsplit(w,jw,lenl+ofst,len+ofst,nsubrow)!quick sort alg (sorts based on L1norm of submats)
        if (lenl > len) call qsplit(w,jw,lenl,len,nsubrow)!quick sort alg (sorts based on L1norm of submats)
    !ceb from ilut.f90:   if (len > lenl) call qsplit(w,jw,len+ofst,lenl+ofst,neqn)!quick sort alg 
       

!write(6,*)" row ",ii," length of L = ",len 
       if (len> 0) then
!write(6,*)"allocate L(",ii,",",len,")"
          allocate(amat%L%pj(ii)%cols(len),STAT=istat)
          allocate(amat%L%pa(ii)%submat(nsubrow,nsubcol,len),STAT=istat)
          do k=1,len
             !amat%L%pj(ii)%cols(k)=iperm(jw(k))
             amat%L%pj(ii)%cols(k)=jw(k)
             do iisub=1,nsubrow
                do jjsub=1,nsubcol
                   amat%L%pa(ii)%submat(iisub,jjsub,k)=w(iisub,jjsub,k)!copy w into new U elements, U[k]=w[ii+(k-1)]
                end do
             end do
          end do
       end if

       !/*---------------------------------------------------------------------
       !|     apply dropping strategy to U  (entries after diagonal)
       !|--------------------------------------------------------------------*/


       tnorm=L1Norm2(w,ii,nsubrow)!ceb
!if(tnorm>10.) 
write(6,*)"3: tnorm(",ii,")=",tnorm

       len = 0 
       do j=2,lenu
          !if ( fabs(w(ii+j)) > drop6*tnorm ) then
          if ( L1Norm2(w,ii+j-ofst,nsubrow) > drop1*tnorm ) then !use matrix norm
             len=len+1
             jw(ii+len) = jw(ii+j-ofst) 
             !w(ii+len) = w(ii+j) 
             do iisub=1,nsubrow
                do jjsub=1,nsubcol
                   w(iisub,jjsub,ii+len)= w(iisub,jjsub,ii+j-ofst)
                end do
             end do
          end if
       end do

       lenu = len+1
       len=lenu
       if(len>fil_1) len=fil_1   !len = lenu > fil6 ? fil6 : lenu 
       amat%U%nnzrow(ii) = len 
       !if (lenu > len+1) call qsplitC(&w(ii+1), &jw(ii+1), lenu-1, len) 
       !ceb use shift_array?
       !if (len+1 < lenu) call qsplitC(w,jw,lenu,len+ofst,nsubrow,ii+1)!quick sort alg (sorts based on L1norm of submats)
       if (len+1 < lenu) call qsplitC(w,jw,lenu-1,len,nsubrow,ii+1)!quick sort alg (sorts based on L1norm of submats)
       if (len>0) then
!write(6,*)"allocate U(",ii,",",len,")"
          allocate(amat%U%pa(ii)%submat(nsubrow,nsubcol,len),STAT=istat)
          allocate(amat%U%pj(ii)%cols(len),STAT=istat)
       end if













       !/*---------------------------------------------------------------------
       !|     now store U in original coordinates
       !|--------------------------------------------------------------------*/
!write(6,*)"5 ilutp: ii=",ii

     !/*---------------------------------------------------------------------
     !|     copy the rest of U (all but diagonal)
     !|--------------------------------------------------------------------*/
       !memcpy(&amat%U%pa(ii)submat(1), &w(ii+1), (len-1)*sizeof(double)) 
       !/* copy the rest of U */ 
       lenu = 1!diagonal already exists
       do k=2,len !do k=ii+1, ii+len-1
          lenu=lenu+1
!write(6,*)"ii=",ii," k=",k," lenu=",lenu
          !amat%U%pj(ii)%cols(lenu) = iperm(jw(k-ofst))
          !amat%U%pj(ii)%cols(k) = iperm(jw(ii+k-ofst))!ceb
          amat%U%pj(ii)%cols(k) = jw(ii+k-ofst)
          do iisub=1,nsubrow
             do jjsub=1,nsubcol
                !amat%U%pa(ii)%submat(iisub,jjsub,lenu+1)=w(iisub,jjsub,k)!copy w into new U elements, U[k]=w[ii+(k-1)]
                amat%U%pa(ii)%submat(iisub,jjsub,k)=w(iisub,jjsub,ii+k-ofst)!copy w into new U elements, U[k]=w[ii+(k-1)]
             end do
          end do
       end do

       tnorm=L1Norm2(w,ii,nsubrow)
       if(abs(tnorm-DBL_EPSILON) <= DBL_EPSILON*tnorm) then
          do iisub=1,nsubrow
             w(iisub,iisub,ii) = (0.0001+drop1)*tnorm 
          end do
       end if


       !===========================================================
       !//compute L,U for pivot (diagonal of combined LU matrix)
       ! Lij,Uij submatrices refer only to the diagonal block element
       tmp=0.0
       !dp=iau(ii)
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
             amat%U%pa(ii)%submat(jjsub,kksub,1)=0.0
          end do
          amat%U%pa(ii)%submat(jjsub,jjsub,1)=1.0
       end do 

       !//write LU into U[diagp]
       do kksub=1,nsubcol !loop subcols
          do jjsub=1,nsubcol-1 !get D[1,..,neqn-1] ! loop subcols
             d(jjsub)=amat%U%pa(ii)%submat(jjsub,kksub,1) !Identity matrix row jj 
             do iisub=1,jjsub-1
                d(jjsub) = d(jjsub) - Lij(iisub,jjsub)*d(iisub)
             end do
             d(jjsub) = d(jjsub)*Lij(jjsub,jjsub)!ceb not sure if this is right
          end do
          do jjsub=nsubcol,1,-1 !get alu[diagp,jj,kk] ! loop subcols
             if(jjsub.eq.nsubcol) then     
                !D=D-L*D (left of subdiag)
                do iisub=1,jjsub-1
                   amat%U%pa(ii)%submat(jjsub,kksub,1) = amat%U%pa(ii)%submat(jjsub,kksub,1)  - Lij(iisub,jjsub)*d(iisub)
                end do
                amat%U%pa(ii)%submat(jjsub,kksub,1) = amat%U%pa(ii)%submat(jjsub,kksub,1)*Lij(jjsub,jjsub)
                !D=
             else      
                !D=D-U*D (right of sub diag)
                amat%U%pa(ii)%submat(jjsub,kksub,1) = d(jjsub)
                do iisub=jjsub+1,nsubrow
                   amat%U%pa(ii)%submat(jjsub,kksub,1) = amat%U%pa(ii)%submat(jjsub,kksub,1) - &
                        Uij(iisub,jjsub)*amat%U%pa(ii)%submat(iisub,kksub,1)
                end do
             end if
          end do
       end do
       
       !U%ja(ii)%cols(1) = jw(ii) !rowj(0) = jw(ii);
       !===========================================================
       amat%U%pj(ii)%cols(1) = ii

    end do



  do i=1,amat%n
     do j=1,amat%L%nnzrow(i)
        write(96,*) i,amat%L%pj(i)%cols(j)
     end do
     do k=1,amat%U%nnzrow(i)
        write(96,*) i,amat%U%pj(i)%cols(k)
     end do
  end do















else!===========================================================
































  do j=1,rmax
     jwrev(j) = -1
     jwrev2(j) = -1
  end do

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
if(tnorm>10.) write(6,*) "pilu: row ",ii,"tnorm=",tnorm
     if(abs(tnorm-DBL_EPSILON) .le. DBL_EPSILON*tnorm) goto  9996

     !/*---------------------------------------------------------------------
     !|     unpack B-block in arrays w, jw, jwrev
     !|     WE ASSUME THERE IS A DIAGONAL ELEMENT
     !|--------------------------------------------------------------------*/
     lenu = 1
     lenl = 0!lenl does not include diagonal?
     jw(ii) = ii
     jwrev(ii) = ii
     call blank_3D_submat(w,ii,nsubrow,nsubcol) !w(ii) = 0.0!zero submatrix
     do j=1, lrowz
        jcol = lrowj%cols(j)
        if (jcol .lt. ii) then !left of diag
           lenl=lenl+1
           jpos=lenl!ceb
           jw(jpos) = jcol 
           jwrev(jcol) = jpos!ceb
        else if (jcol .eq. ii) then !on diag
           jpos=ii!ceb 
        else                        !right of diag
           jpos = ii+lenu!ceb
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
     if (test.eq.1 .and. C%n.ne.0) then
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
     end if!test
     
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
!write(6,*)"exchanging col ",k," with ",jj   
           !/*   exchange in jw   */
           j = jw(jj)
           jw(jj) = jw(k)
           jw(k) = j
           !/*   exchange in jwrev   */
           jwrev(jrow) = jj
           jwrev(j) = k
           !/*   exchange in w   */
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
        lrowm => amat%U%pa(jrow)
        lrowj => amat%U%pj(jrow)
        lrowz = amat%U%nnzrow(jrow)
        !fact = w(jj) * U(0)!submatrix mult
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              fact(iisub,jjsub)=0.0
              do kksub=1,nsubrow !assumes square submatrix
                 fact(iisub,jjsub) = fact(iisub,jjsub) + w(iisub,kksub,jj)*lrowm%submat(kksub,jjsub,1) !equivalent L(ii,jj)*U(diagnj) (pivot)
              end do
           end do
        end do
        
!if(L1Norm(fact,nsubrow)>10)
!write(6,*)"2: ii:",ii," pivot norm(",jrow,")=",L1Norm(fact,nsubrow)

        if (L1Norm(fact,nsubrow) <= drop0 ) then !/*   DROPPING IN L   */
           !write(6,*)"dropping in L  row,column(",ii,",",jrow,")"
           cycle
        end if
        
        !/*---------------------------------------------------------------------
        !|     combine current row and row jrow
        !|--------------------------------------------------------------------*/
        do k=2,lrowz
           !s = fact * rowm(k)
           jcol = lrowj%cols(k)
           jpos = jwrev(jcol)  
           if (jpos .eq. -1) then 
              !/*---------------------------------------------------------------------
              !|     this is a fill-in element
              !|--------------------------------------------------------------------*/
              if (jcol .ge. ii) then ! right of diag
                 !/*---------------------------------------------------------------------
                 !|     dealing with U
                 !|--------------------------------------------------------------------*/
                 if(lenu < fil_1) then
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
                 if(lenl < fil_0) then
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
if (test.eq.1 .and. C%n.ne.0) then  

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
           if (jpos .eq. -1) then
              len2=len2+1
              jpos=len2
              jw2(jpos) = jcol
              jwrev2(jcol) = jpos
              call blank_3D_submat(w2,jpos,nsubrow,nsubcol)
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
        
end if!test            







     end do!end loop over block cols


     !/*---------------------------------------------------------------------
     !|     reset nonzero indicators
     !|--------------------------------------------------------------------*/
     do j=1,len2    !/*  L^{-1} F block  */
        jwrev2(jw2(j)) = -1
     end do
if (test.eq.1) then
     do j=1,lenl    !/*  L block  */
        jwrev(jw(j)) = -1
     end do
end if
     do j=1,lenu    !/*  U block  */
        jwrev(jw(ii+j-ofst)) = -1
     end do

     !/*---------------------------------------------------------------------
     !|     done reducing this row, now store L
     !|--------------------------------------------------------------------*/
     !lenl=len
     len=lenl
     !if (lenl>fil_0) lenl= fil_0!lenl = len > fil0 ? fil0 : len;
     if (lenl>fil_0) len= fil_0!lenl = len > fil0 ? fil0 : len;
     !len = lenl!ceb
     amat%L%nnzrow(ii) = len
     !amat%L%nnzrow(ii) = lenl

     !if (lenl < len) call qsplitC(w, jw, len, lenl)
     !/* quick sort */
     !if (lenl < len) call qsplit(w,jw,len,lenl,nsubrow)!quick sort alg (sorts based on L1norm of submats)
     if (lenl > len) call qsplit(w,jw,lenl,len,nsubrow)!quick sort alg (sorts based on L1norm of submats)
     !if (lenl > len) call qsplit(w,jw,lenl+ofst,len+ofst,nsubrow)!quick sort alg (sorts based on L1norm of submats)
     !if (lenl > 0) then

     if (len > 0) then
!write(6,*)"allocate L(",ii,",",lenl,")"
        !allocate(amat%L%pj(ii)%cols(lenl),STAT=istat)
        allocate(amat%L%pj(ii)%cols(len),STAT=istat)
        !allocate(amat%L%pa(ii)%submat(nsubrow,nsubcol,lenl),STAT=istat)
        allocate(amat%L%pa(ii)%submat(nsubrow,nsubcol,len),STAT=istat)
        !do k=1,lenl
        do k=1,len
           amat%L%pj(ii)%cols(k)=jw(k)!copy jw into new U_ja elements
           !amat%L%pj(ii)%cols(k)=jwrev(jw(k))!copy jw into new U_ja elements
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
           !jpos=len
           !jw(jpos) = jw(ii+j-ofst)
           jw(ii+len) = jw(ii+j-ofst) 
           !ceb ilut:  jw(ii+len) = jw(ii+(k-ofst))
          !w(len) = w(ii+j)!copy submatrix
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 !w(iisub,jjsub,jpos)= w(iisub,jjsub,ii+j-ofst)
                 w(iisub,jjsub,ii+len)= w(iisub,jjsub,ii+j-ofst)
                 !ceb ilut: w(iisub,jjsub,(ii+len)) = w(iisub,jjsub,ii+k-ofst)
             end do
           end do
        !else
           !write(6,*) "pilu checkpoint 31 :Dropping small column element  norm(w(", ii+j-ofst, "))=", L1Norm2(w,ii+j-ofst,nsubrow),&
           !     " < tol*tnorm(",ii,")=",drop1*tnorm
        end if
     end do

     lenu=len+1!assign new length after dropping
     len=lenu!ceb
     !if (lenu > fil_1) lenu=fil_1 !lenu = len+1 > fil_1 ? fil_1 : len+1;
     if (len > fil_1) len=fil_1 !lenu = len+1 > fil_1 ? fil_1 : len+1;
     amat%U%nnzrow(ii) = len

     jpos = lenu !jpos = lenu - 1
     !if (jpos < len+ofst) call qsplit(w,jw,len+ofst,jpos,nsubrow)!quick sort alg (sorts based on L1norm of submats)
     if (len+1 < lenu) call qsplitC(w,jw,lenu-1,len,nsubrow,ii+1)!quick sort alg (sorts based on L1norm of submats)
     if (lenu>0) then
     !if (len>0) then
!        write(6,*)"allocate U(",ii,",",lenu,")"
        !allocate(amat%U%pa(ii)%submat(nsubrow,nsubcol,lenu),STAT=istat)
        allocate(amat%U%pa(ii)%submat(nsubrow,nsubcol,len),STAT=istat)
        !allocate(amat%U%pj(ii)%cols(lenu),STAT=istat)
        allocate(amat%U%pj(ii)%cols(len),STAT=istat)
     end if








     !/*---------------------------------------------------------------------
     !|     copy the rest of U (all but diagonal)
     !|--------------------------------------------------------------------*/
     !!do k=1,jpos-ofst! (skips diagonal)
       lenu = 1!diagonal already exists !ceb
     !do k=2,jpos! (skips diagonal)
     do k=2,len! (skips diagonal)
          lenu=lenu+1 !ceb
        !amat%U%pj(ii)%cols(k)=jw(k-ofst)!copy jw into new U_pj elements
        !   amat%U%pj(ii)%cols(k) = iperm(jw(ii+k-ofst))
        amat%U%pj(ii)%cols(k)=jw(ii+k-ofst)!copy jw into new U_pj elements
        !amat%U%pj(ii)%cols(k)=jwrev(jw(ii+k-ofst))!copy jw into new U_pj elements
        do iisub=1,nsubrow
           do jjsub=1,nsubrow
               !amat%U%pa(ii)%submat(iisub,jjsub,k)=w(iisub,jjsub,k-ofst)!copy w into new U elements, U[k]=w[ii+(k-1)]
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
           d(jjsub)=amat%U%pa(ii)%submat(jjsub,kksub,1) !Identity matrix row jj 
           do iisub=1,jjsub-1
              d(jjsub) = d(jjsub) - Lij(iisub,jjsub)*d(iisub)
           end do
           d(jjsub) = d(jjsub)*Lij(jjsub,jjsub)!ceb not sure if this is right
        end do
        
        do jjsub=nsubcol,1,-1 !get alu[diagp,jj,kk] ! loop subcols
           if(jjsub.eq.nsubcol) then     
              !D=D-L*D (left of subdiag)
              do iisub=1,jjsub-1
                 amat%U%pa(ii)%submat(jjsub,kksub,1) = amat%U%pa(ii)%submat(jjsub,kksub,1)  - Lij(iisub,jjsub)*d(iisub)
              end do
              amat%U%pa(ii)%submat(jjsub,kksub,1) = amat%U%pa(ii)%submat(jjsub,kksub,1)*Lij(jjsub,jjsub)
              !D=
           else      
              !D=D-U*D (right of sub diag)
              amat%U%pa(ii)%submat(jjsub,kksub,1)=d(jjsub)
              do iisub=jjsub+1,nsubrow
                 amat%U%pa(ii)%submat(jjsub,kksub,1) = amat%U%pa(ii)%submat(jjsub,kksub,1) - &
                      Uij(iisub,jjsub)*amat%U%pa(ii)%submat(iisub,kksub,1)
              end do
           end if
        end do
     end do
   
     !U%ja(ii)%cols(1) = jw(ii) !rowj(0) = jw(ii);
     !===========================================================
     amat%U%pj(ii)%cols(1) = ii!ceb moved from prior to diag formation

     !ceb may have to include computation of [E U^{-1}] and  [L^{-1} F] here since U(ii,1)
     !now contains portions of L and diag 


     !===========================================================








     !write(6,*)"pilu checkpoint 38"
     
     !ceb need to include diagonal contribution of L^-1 to F
     if (test.eq.1 .and. C%n.ne.0) then
        !/*---------------------------------------------------------------------
        !|     copy  L^{-1} F
        !|--------------------------------------------------------------------*/
        len = 0
        write(6,*)"pilu checkpoint 38.1 len2=",len2
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
        if(lenu>fil_2) lenu=fil_2 !lenu = len > fil2 ? fil2 : len;
        !if(lenu.eq.0) write(6,*)"pilu checkpoint 38.1 ii=",ii," lenu=",lenu
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
     end if!test
     

  end do !end loop over rows

!write(6,*)"pilu checkpoint 39"










  !if (test.eq.0 .and. C%n.ne.0) then
  if (test.eq.1 .and. C%n.gt.0 ) then
     
     write(6,*)"pilu checkpoint 40"
     
     !/*---------------------------------------------------------------------
     !|    beginning of second main loop   E U^{-1} and Schur complement
     !|--------------------------------------------------------------------*/
     do ii=1,rsize
        lrowj => amat%E%pj(ii)
        lrowm => amat%E%pa(ii)
        lrowz =  amat%E%nnzrow(ii)
        
        rrowj => C%pj(ii)
        rrowm => C%pa(ii)
        rrowz =  C%nnzrow(ii)
        
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
              !/*   exchange in w   */
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
           
           if (L1Norm(fact,nsubrow) <= drop3 ) then !/*  DROPPING IN E U^{-1}   */
              !write(6,*)"DROPPING IN E U^{-1}: row,column(",ii,",",jrow,")"
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
                 if (lenl .lt. fil_0) then
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
              if (jpos .eq. -1) then
                 !if(lenu .lt. fil_1) then
                 if(lenu .lt. fil_1) then
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
        !write(6,*)"pilu checkpoint 44.6: amat%F%nnzrow(",ii,")=",amat%F%nnzrow(ii)
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
        
        !ceb pretty sure this schur diagonal is not correct (missing LU combined diag contribution)
        
        if(tnorm>10.0) write(6,*)"pilu checkpoint 45: Schur tnorm(",ii,")=",tnorm
        
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
        
        !write(6,*)"pilu checkpoint 46 len(",ii,")=",len
        lenu=len!+ofst!ceb
        if(lenu>fil_4) lenu=fil_4 !lenu = len > fil4 ? fil4 : len;
        schur%nnzrow(ii) = lenu
        jpos = lenu
        !write(6,*)"pilu checkpoint 46.5: lenu(",ii,")=",lenu
        
        !if (jpos < len+ofst) call qsplit(w,jw,len+ofst,jpos,nsubrow)
        if (jpos < len) then
           !write(6,*) "pilu checkpt 47: calling qsplit"
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
           do iisub=1,nsubrow
              do jjsub=1,nsubrow
                 schur%pa(ii)%submat(iisub,jjsub,k)=w(iisub,jjsub,k)!copy w into new schur elements, schur[k]=w[ii+(k-1)]?
              end do
           end do
        end do
     end do !end loop over rows
     




  do i=1,lsize
     do j=1,amat%L%nnzrow(i)
        write(96,*) i,amat%L%pj(i)%cols(j)
     end do
     do j=1,amat%U%nnzrow(i)
        write(96,*) i,amat%U%pj(i)%cols(j)
     end do
  end do




     
  end if!test 
  








  !/*---------------------------------------------------------------------
  !|     end main loop - now do cleanup
  !|--------------------------------------------------------------------*/
  nullify(lrowj)
  nullify(lrowm)
  nullify(rrowj)
  nullify(rrowm)
write(6,*)"pilu checkpoint 49.5"

  deallocate(d,STAT = istat)
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
end if


  !do i=1,lsize
  !      write(96,*) i,"L:",(amat%L%pj(i)%cols(j),j=1,amat%L%nnzrow(i)),"U:",(amat%U%pj(i)%cols(j),j=1,amat%U%nnzrow(i))
  !      !write(76,*) i,"U:",(amat%U%pj(i)%cols(j),j=1,amat%U%nnzrow(i))
  !end do



  !do i=1,schur%n
  !   do j=1,schur%nnzrow(i)
  !      write(78,*) i,schur%pj(i)%cols(j) 
  !   end do
  !end do

  !do j=1,schur%n
  !      write(79,*) j,(schur%pj(j)%cols(k) ,k=1,schur%nnzrow(j))
  !end do



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
  write(6,*) "zero row encountered index= ",ii+1
  return! 6;
end subroutine pilu
!/*---------------------------------------------------------------------
!|     end of pilut
!|--------------------------------------------------------------------*/

end module pilu_lib
