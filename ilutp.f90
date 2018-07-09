module ilutp
  use kinddefs
  !use linalg
  use blas
  use armslib
  implicit none
contains
  
  !#define DBL_EPSILON 2.2204460492503131e-16 // double epsilon
  subroutine ilutpC(amat, droptol, lfil, permtol, mband, ilusch)
    !{
    !/*---------------------------------------------------------------------- 
    !| ILUTP -- ILUT with column pivoting -- adapted from ILUTP [Sparskit] 
    !| Converted to C so that dynamic memory allocation may be implememted.
    !| All indexing is in C format.
    !|----------------------------------------------------------------------
    !| ILUT factorization with dual truncation. 
    !|----------------------------------------------------------------------
    !|
    !| on entry:
    !|========== 
    !| ( amat ) =  Matrix stored in SparRow struct.
    !|
    !| lfil[6]  =  number nonzeros in L-part
    !| lfil[7]  =  number nonzeros in U-part     ( lfil >= 0 )
    !|
    !| droptol[5] = threshold for dropping small terms in L during
    !|              factorization.
    !| droptol[6] = threshold for dropping small terms in U.
    !|
    !| permtol  =  tolerance ratio used to  determine whether or not to permute
    !|             two columns.  At step i columns i and j are permuted when
    !|
    !|                     abs(a(i,j))*permtol > abs(a(i,i))
    !|
    !|           [0 --> never permute; good values 0.1 to 0.01]
    !|
    !| mband    =  permuting is done within a band extending to mband
    !|	      diagonals only. 
    !|             mband = 0 --> no pivoting. 
    !|	      mband = n --> pivot is searched in whole column 
    !|
    !|
    !| On return:
    !|===========
    !|
    !| (ilusch) =  Contains L and U factors in an LUfact struct.
    !|             Individual matrices stored in SparRow structs.
    !|             On return matrices have C (0) indexing.
    !|
    !| iperm    =  reverse permutation array.
    !|
    !|       integer value returned:
    !|
    !|             0   --> successful return.
    !|             1   --> Error.  Input matrix may be wrong.  (The 
    !|                         elimination process has generated a
    !|                         row in L or U whose length is > n.)
    !|             2   --> Memory allocation error.
    !|             5   --> Illegal value for lfil.
    !|             6   --> zero row encountered.
    !|----------------------------------------------------------------------- 
    !| work arrays:
    !|=============
    !| jw, jwrev = integer work arrays of length n
    !| w         = real work array of length n. 
    !|----------------------------------------------------------------------- 
    !|    
    !|--------------------------------------------------------------------*/
    type(cs_type),pointer:: amat !should already be allocated
    type(ilut_type),pointer::ilusch 
    real(dp),intent(in) :: droptol(7)
    integer,dimension(:)::lfil
    real(dp):: permtol
    integer:: mband

    !local vars
    integer :: i, ii, j, jj, jcol, jpos, jrow, k, len, lenu, lenl
    integer :: iisub,jjsub,kksub,itmp,istat,test,idiag,avg_nnz,max_nnz
    integer::rmax,lrowz,imax,icut,dec,decnt,fil_SL,fil_SU,nsubrow,nsubcol
    real(dp) :: tnorm,wnorm,xmax,xmax0,tmp,w_norm 
    integer,dimension(:),allocatable::jwrev,iprev 

    integer,target,dimension(:),allocatable :: jw 
    real(dp),target,dimension(:,:,:),allocatable :: w
 
    type(ja_type),pointer:: lrowj !=NULL,  
    integer,pointer,dimension(:) :: iperm
    type(submat_type),pointer::lrowm !=NULL

    ! temporary vars for copying submatrix data
    real(dp), dimension(:,:),allocatable:: Uij,Lij,s,t,fact
    real(dp), dimension(:),allocatable:: d

    integer,pointer,dimension(:):: shift_jw
    real(dp),pointer,dimension(:,:,:):: shift_w

    real(dp)::drop5,drop6
    nsubrow=amat%nsubrow
    nsubcol=amat%nsubcol
    iperm=>ilusch%perm2 

write(6,*)"calling ilutp"

  !write(6,*)"ilutp: droptol[1,2,3,4,5,6,7]=",droptol(1),droptol(2),droptol(3),droptol(4),droptol(5),droptol(6),droptol(7)

    decnt=0 

    !ceb compute avg nnzrow for amat
    avg_nnz=0
    max_nnz=0
    do i=1,amat%n
       avg_nnz=avg_nnz+amat%nnzrow(i)
       if(amat%nnzrow(i)>max_nnz) max_nnz=amat%nnzrow(i)
    end do
    avg_nnz=avg_nnz/amat%n
    !avg_nnz=max_nnz

    !ceb define fill (number of nonzeros allowed) = [1+lfil()]*avg_nnz
    !     so that lfil=0 results in fil=avg_nnz per row
    !  this allows us to avoid having to change lfil as the size of the input matrix changes (and as we descend levels)
    !ceb though currently defined as int lfil could be real so that fil_() does not have to be multiple of avg_nnz
    fil_SL= ((1+lfil(6))*avg_nnz)  
    fil_SU= ((1+lfil(7))*avg_nnz)
    drop5=droptol(6)
    drop6=droptol(7)

    !write(6,*) "ilutpC: fil_SL=",fil_SL," fil_SU=",fil_SU," avg_nnz=",avg_nnz

    test=0!ceb partial pivoting currently turned off

    rmax = amat%n 
    ilusch%n = rmax
!write(6,*)"0 ilutp: rmax=",rmax

    if (rmax .eq. 0) then
       return!(0)  
    else!    if (rmax > 0) then 
       allocate(Uij(nsubrow,nsubcol),STAT=istat)
       allocate(Lij(nsubrow,nsubcol),STAT=istat)
       allocate(s(nsubrow,nsubcol),STAT=istat)
       allocate(t(nsubrow,nsubcol),STAT=istat)
       allocate(fact(nsubrow,nsubcol),STAT=istat)
       allocate(d(nsubrow),STAT=istat)
       allocate(jw(rmax),STAT=istat) 
       allocate(w(nsubrow,nsubcol,rmax),STAT=istat) 
       allocate(jwrev(rmax),STAT=istat) 
       allocate(iprev(rmax),STAT=istat) 
    end if


!write(6,*)"ilutp: checkpt 1 after temp array allocation"



!  do i=1,amat%n
!     do j=1,amat%nnzrow(i)
!  !!!      write(85,*) i,amat%pj(i)%cols(j)
!        write(184,*) i,amat%pj(i)%cols(j),amat%pd(i),((amat%pa(i)%submat(iisub,jjsub,j),iisub=1,amat%nsubrow),jjsub=1,amat%nsubcol)
!  !      write(84,*) i,amat%pd(i),amat%pj(i)%cols(j),amat%nnzrow(i)
!     end do
!  end do



    if (fil_SL<0 .or. fil_SU<0 .or. rmax<=0) goto  998 

    do j=1,rmax 
       jwrev(j) = -1 
       iprev(j) = j 
       !ilusch%L%nnzrow(j)=0
       !ilusch%U%nnzrow(j)=0
    end do

    !/*---------------------------------------------------------------------
    !|    beginning of main loop - L, U calculations
    !|--------------------------------------------------------------------*/
    do ii=1,rmax 
write(6,*) "ilutp: computing row ",ii
       lrowj => amat%pj(ii) 
       lrowm => amat%pa(ii) 
       lrowz =  amat%nnzrow(ii) 
       !/*---------------------------------------------------------------------
       !|   check for zero row 
       !|--------------------------------------------------------------------*/
       tnorm = 0.0 
       do k=1,lrowz
          tnorm = tnorm + L1Norm2(lrowm%submat,k,nsubrow) !tnorm = tnorm + abs(rowm(k)) 
       end do
       tnorm = tnorm / lrowz 
!if(tnorm>10.) write(6,*) "1: row ",ii,"tnorm=",tnorm
       if(abs(tnorm-DBL_EPSILON) .le. DBL_EPSILON*tnorm) goto  9991 

       !/*---------------------------------------------------------------------
       !|     unpack amat in arrays w, jw, jwrev
       !|     WE ASSUME THERE IS A DIAGONAL ELEMENT
       !|--------------------------------------------------------------------*/
       lenu = 1 
       lenl = 0 
       jw(ii) = ii 
       jwrev(ii) = ii 
       call blank_3D_submat(w,ii,nsubrow,nsubcol) !w(ii) = 0.0 
       do j=1,lrowz 
          jcol = iprev(lrowj%cols(j)) 
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

!write(6,*)"ilutp: row ",ii," loading jw(",jpos,")=",jw(jpos)


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

!write(6,*)"ilutp: row ",ii," looping column jw(",jj,")=",jrow

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
!write(6,*)"row ",ii," exchanging col ",k," with ",jj   
             !/*   exchange in jw   */
             j = jw(jj) 
             jw(jj) = jw(k) 
             jw(k) = j 

!write(6,*)"ilutp: jw(",jj,")=",jw(jj)," jw(",k,")=",jw(k)

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
          lrowm => ilusch%U%pa(jrow) 
          lrowj => ilusch%U%pj(jrow) 
          lrowz =  ilusch%U%nnzrow(jrow) 
          !fact = w(jj) * rowm(0) = row A(ii)*U_d(ii)
          do iisub=1,nsubrow
             do jjsub=1,nsubcol
                fact(iisub,jjsub)=0.0
                do kksub=1,nsubrow !assumes square submatrix
                   fact(iisub,jjsub) = fact(iisub,jjsub) + w(iisub,kksub,jj)*lrowm%submat(kksub,jjsub,1)
                end do
             end do
          end do

!if(L1Norm(fact,nsubrow)>10)write(6,*)"2:ilutp ii:",ii," pivot norm(",jrow,")=",L1Norm(fact,nsubrow)

          if (L1Norm(fact,nsubrow) <= drop5 ) then  !/*  DROPPING IN L  */
!write(6,*)"ilutp dropping norm ",L1Norm(fact,nsubrow)," in L row,column(",ii,",",jrow,")",L1Norm(fact,nsubrow),"<=",drop5
             cycle
          end if

          idiag=amat%pd(ii)
          if(idiag<=0)idiag=1
          !/*---------------------------------------------------------------------
          !|     combine current row and row jrow
          !|--------------------------------------------------------------------*/
          do k=2,lrowz 
             jcol = iprev(lrowj%cols(k))    !/*  new column number  */
             jpos = jwrev(jcol)        
!write(6,*)"ilutp: combining row ",ii," and column ",jcol

             if(jpos == -1) then 
                !/*---------------------------------------------------------------------
                !|     this is a fill-in element
                !|--------------------------------------------------------------------*/
                if (jcol >= ii) then
                   !/*---------------------------------------------------------------------
                   !|     dealing with U
                   !|--------------------------------------------------------------------*/
!write(6,*)"ilutp: amat%nnzrow(",ii,")-amat%pd(",ii,")+1=",amat%nnzrow(ii)-amat%pd(ii)+1

                   !if(lenu < fil_SU) then
                      jpos = ii + lenu
                      lenu=lenu+1 

                      jw(jpos) = jcol
                      jwrev(jcol) = jpos 
                      call blank_3D_submat(w,jpos,nsubrow,nsubcol)
                   !end if
                else !left of diag
                   !/*---------------------------------------------------------------------
                   !|     dealing  with L
                   !|--------------------------------------------------------------------*/
!write(6,*)"ilutp: amat%pd(",ii,")-1=",amat%pd(ii)-1
                   !if(lenl < fil_SL) then
                      lenl=lenl+1 
                      jpos=lenl

                      jw(jpos) = jcol 
                      jwrev(jcol) = jpos 
                      call blank_3D_submat(w,jpos,nsubrow,nsubcol)
                   !end if
                end if
             !else
                !/*---------------------------------------------------------------------
                !|     this is not a fill-in element 
                !|--------------------------------------------------------------------*/
             end if

             !perform row operation w(ii)-=pivot*U(k)
             if(jpos.ne.-1) then

!write(6,*)"ilutp: col index ",k," new or updated jw(",jpos,")=",jw(jpos)


                do iisub=1,nsubrow
                   do jjsub=1,nsubcol
                      do kksub=1,nsubrow
                         w(iisub,jjsub,jpos) = w(iisub,jjsub,jpos) - fact(iisub,kksub)*lrowm%submat(kksub,jjsub,k)
                      end do
                   end do
                end do

             end if
             
             if ((lenu .gt. rmax) .or. (lenl .gt. rmax)) goto 994 !cycle?
             !if ((lenu .gt. fil_SU+(amat%nnzrow(ii)-idiag+1)) .or. (lenl .gt. fil_SL)) goto 994 !cycle?
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
       !|--------------------------------------------------------------------*/
       len = lenl
       !write(6,*)"2: ilutp: amat%pd(",ii,")-1=",amat%pd(ii)-1
       !idiag=amat%pd(ii)
       !if(idiag<=0)idiag=1
       if(lenl>fil_SL) len=fil_SL !len = lenl > fil5 ? fil5 : lenl 
       !if (lenl > (fil_SL+idiag-1)) len = fil_SL+idiag-1
       !if (lenl > (fil_SL+amat%nnzrow(ii)-idiag+1)) len = fil_SL+amat%nnzrow(ii)-idiag+1

       ilusch%L%nnzrow(ii) = len 
       !write(6,*)"L%nnzrow(",ii,")=",ilusch%L%nnzrow(ii)
       !/* quick sort */
!write(6,*)"ilutp: qsplit checkpt 1 lenl=",lenl," len=",len
       if (lenl > len) call qsplit(w,jw,lenl,len,nsubrow)!quick sort alg (sorts based on L1norm of submats)
       if (len> 0) then
!write(6,*)"allocate L(",ii,",",len,")"
          allocate(ilusch%L%pj(ii)%cols(len),STAT=istat)
          allocate(ilusch%L%pa(ii)%submat(nsubrow,nsubcol,len),STAT=istat)
          do k=1,len
             ilusch%L%pj(ii)%cols(k)=iperm(jw(k))
             do iisub=1,nsubrow
                do jjsub=1,nsubcol
                   ilusch%L%pa(ii)%submat(iisub,jjsub,k)=w(iisub,jjsub,k)!copy w into new U elements, U[k]=w[ii+(k-1)]
                end do
             end do
          end do
       end if

       !/*---------------------------------------------------------------------
       !|     apply dropping strategy to U  (entries after diagonal)
       !|--------------------------------------------------------------------*/

       tnorm=L1Norm2(w,ii,nsubrow)!ceb
!if(tnorm>10.)write(6,*)"3: tnorm(",ii,")=",tnorm

       len = 0 
       do j=2,lenu
          jcol=ii+j-ofst
          if ( L1Norm2(w,jcol,nsubrow) > drop6*tnorm ) then !use matrix norm
             len=len+1
             jw(ii+len) = jw(jcol) 

!write(6,*)"ilutpC: dropping in U  jw(",ii+len,")=",jw(ii+len)

             do iisub=1,nsubrow
                do jjsub=1,nsubcol
                   w(iisub,jjsub,ii+len)= w(iisub,jjsub,jcol) !w(ii+len) = w(ii+j) 
                end do
             end do
          else
           !write(6,*) "ilutp :U Dropping small column element norm(w(", jcol, "))=", L1Norm2(w,jcol,nsubrow),&
           !     " < tol*tnorm(",ii,")=",drop6*tnorm
          end if
       end do

       lenu = len+1
       len=lenu
       if(lenu>fil_SU) len=fil_SU   !len = lenu > fil6 ? fil6 : lenu 

!write(6,*)"2 ilutp: amat%nnzrow(",ii,")-amat%pd(",ii,")+1=",amat%nnzrow(ii)-amat%pd(ii)+1

       !idiag=amat%pd(ii)
       !if(idiag<=0)idiag=1
       !if(lenu>(fil_SU+amat%nnzrow(ii)-idiag+1)) len=fil_SU+(amat%nnzrow(ii)-idiag+1)   !len = lenu > fil6 ? fil6 : lenu 
       ilusch%U%nnzrow(ii) = len 

!write(6,*)"U%nnzrow(",ii,")=",ilusch%U%nnzrow(ii),"  L%nnzrow(",ii,")=",ilusch%L%nnzrow(ii),"total LU nnzrow=",&
!              ilusch%U%nnzrow(ii)+ilusch%L%nnzrow(ii)
!write(6,*)"ilutpC: total LU(",ii,") nnzrow=",ilusch%U%nnzrow(ii)+ilusch%L%nnzrow(ii)

!do i=1,ii+len-1
!write(6,*)"ilutpC: pre shift jw(",i,")=",jw(i)
!end do


!write(6,*)"ilutp: qsplit checkpt 3 len=",len," lenu=",lenu
       shift_w=>w(:,:,ii+1:)
       shift_jw=>jw(ii+1:)
       !if (len+1 < lenu) call qsplit(shift_w,shift_jw,lenu-1,len,nsubrow)!quick sort alg (sorts based on L1norm of submats)
       !if (len < lenu) call qsplit(shift_w,shift_jw,lenu,len,nsubrow)!quick sort alg (sorts based on L1norm of submats)
       if (len < lenu) call qsplit(shift_w,shift_jw,lenu-1,len,nsubrow)!quick sort alg (sorts based on L1norm of submats)
       nullify(shift_w)
       nullify(shift_jw)
!write(6,*)"ilutp: qsplit checkpt 4 len=",len


!do i=1,ii+len-1
!write(6,*)"ilutpC: new shifted jw(",i,")=",jw(i)
!end do







       if (len>0) then
!write(6,*)"allocate U(",ii,",",len,")"
          allocate(ilusch%U%pa(ii)%submat(nsubrow,nsubcol,len),STAT=istat)
          allocate(ilusch%U%pj(ii)%cols(len),STAT=istat)
       end if





if(test.eq.0) then
       !/*---------------------------------------------------------------------
       !|     determine next pivot
       !|--------------------------------------------------------------------*/
       !/* HERE  - all lines with dec included for counting pivots */
       imax = ii 
       xmax = L1Norm2(w,imax,nsubrow)
!write(6,*) "ilutp: check pivot for row ",ii," diag norm:",xmax
       xmax0 = xmax 
       !/* icut = ii - 1 + mband - ii % mband   */ 
       icut = ii-1 + mband  
       dec = 0 
       do k=ii+1, ii+len-ofst
          wnorm = L1Norm2(w,k,nsubrow) 
!write(6,*)"ilutp: checking col: jw(",k,"):",jw(k)," elem norm:",wnorm
          if ( (wnorm > xmax) .and. (wnorm*permtol > xmax0) .and. (jw(k) .le. icut) ) then
             imax = k 
!write(6,*)"ilutp:  found new pivot: row ",ii," col jw(",imax,")=",jw(imax)," wnorm=",wnorm,"> xmax=",xmax
             xmax = wnorm 
             dec = 1 
          end if
       end do
       if (dec .eq. 1) then
          decnt=decnt+1 
          !/*---------------------------------------------------------------------
          !|     exchange w's
          !|--------------------------------------------------------------------*/
          !ceb this causes errors in solution of big test matrix (maybe from block structure error?)
          call swap_submat(w,ii,imax,nsubrow)

          !/*---------------------------------------------------------------------
          !|     update iperm and reverse iperm (iprev)
          !|--------------------------------------------------------------------*/
          j = jw(imax) 
          i = iperm(ii)
          iperm(ii) = iperm(j)
          iperm(j) = i
          iprev(iperm(ii)) = ii
          iprev(iperm(j)) = j
       end if
end if



       !/*---------------------------------------------------------------------
       !|     now store U in original coordinates
       !|--------------------------------------------------------------------*/
!write(6,*)"ilutp: checkpt, store off diag U"

       wnorm=L1Norm2(w,ii,nsubrow)
       if(abs(wnorm-DBL_EPSILON) <= DBL_EPSILON*wnorm) then
          do iisub=1,nsubrow
             w(iisub,iisub,ii) = (0.0001+drop6)*tnorm 
          end do
       end if

     !/*---------------------------------------------------------------------
     !|     copy the rest of U (all but diagonal)
     !|--------------------------------------------------------------------*/
       !memcpy(&ilusch%U%pa(ii)submat(1), &w(ii+1), (len-1)*sizeof(double)) 
       !/* copy the rest of U */ 
       lenu = 1!diagonal already exists
       !do k=2,len 
       do k=ii+1, ii+len-ofst
          lenu=lenu+1

!write(6,*)"jw(",k,")=",jw(k),"  len=",len

          ilusch%U%pj(ii)%cols(lenu) = iperm(jw(k))
          do iisub=1,nsubrow
             do jjsub=1,nsubcol
                ilusch%U%pa(ii)%submat(iisub,jjsub,k-ii+1)=w(iisub,jjsub,k)!copy w into new U elements, U[k]=w[ii+(k-1)]
             end do
          end do
       end do

!write(6,*)"ilutp: checkpt, compute diag"



       !===========================================================
       !//compute L,U for pivot (diagonal of combined LU matrix)
       ! Lij,Uij submatrices refer only to the diagonal block element
       tmp=0.0
       !diagp=iau(ii)
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
             ilusch%U%pa(ii)%submat(jjsub,kksub,1)=0.0
          end do
          ilusch%U%pa(ii)%submat(jjsub,jjsub,1)=1.0
       end do 

       !//write LU into U[diagp]
       do kksub=1,nsubcol !loop subcols
          do jjsub=1,nsubcol-1 !get D[1,..,neqn-1] ! loop subcols
             d(jjsub)=ilusch%U%pa(ii)%submat(jjsub,kksub,1) !Identity matrix row jj 
             do iisub=1,jjsub-1
                d(jjsub) = d(jjsub) - Lij(iisub,jjsub)*d(iisub)
             end do
             d(jjsub) = d(jjsub)*Lij(jjsub,jjsub)!ceb not sure if this is right
          end do
          do jjsub=nsubcol,1,-1 !get alu[diagp,jj,kk] ! loop subcols
             if(jjsub.eq.nsubcol) then     
                !D=D-L*D (left of subdiag)
                do iisub=1,jjsub-1
                   ilusch%U%pa(ii)%submat(jjsub,kksub,1) = ilusch%U%pa(ii)%submat(jjsub,kksub,1)  - Lij(iisub,jjsub)*d(iisub)
                end do
                ilusch%U%pa(ii)%submat(jjsub,kksub,1) = ilusch%U%pa(ii)%submat(jjsub,kksub,1)*Lij(jjsub,jjsub)
                !D=
             else      
                !D=D-U*D (right of sub diag)
                ilusch%U%pa(ii)%submat(jjsub,kksub,1) = d(jjsub)
                do iisub=jjsub+1,nsubrow
                   ilusch%U%pa(ii)%submat(jjsub,kksub,1) = ilusch%U%pa(ii)%submat(jjsub,kksub,1) - &
                        Uij(iisub,jjsub)*ilusch%U%pa(ii)%submat(iisub,kksub,1)
                end do
             end if
          end do
       end do
       
       !U%ja(ii)%cols(1) = jw(ii) !rowj(0) = jw(ii);
       ilusch%U%pj(ii)%cols(1) = ii
       !===========================================================

    end do
!write(6,*)"ilutp: checkpt 100 after main loop"


    !permute U,L columns
    if(decnt>0) then !iprev has been permuted
       call cpermC(ilusch%U, iprev)
       call cpermC(ilusch%L, iprev)
    end if
    !/*---------------------------------------------------------------------
    !|     end main loop - now do clean up
    !|--------------------------------------------------------------------*/



!    do i=1,ilusch%n
!!write(6,*)"ilutp row ",i," L%rownnz(i)=",ilusch%L%nnzrow(i)," U%rownnz(i)=",ilusch%U%nnzrow(i)!
!
!       do j=1,ilusch%L%nnzrow(i)
!          write(88,*) i,ilusch%L%pj(i)%cols(j)
!       end do
!       do k=1,ilusch%U%nnzrow(i)
!          write(88,*) i,ilusch%U%pj(i)%cols(k)
!       end do
!   end do



!    do i=1,ilusch%n
!!write(6,*)"ilutc row ",i," L%rownnz(i)=",ilusch%L%nnzrow(i)," U%rownnz(i)=",ilusch%U%nnzrow(i)
!       do j=1,ilusch%L%nnzrow(i)
!          write(188,*) i,ilusch%L%pj(i)%cols(j), &
!               ((ilusch%L%pa(i)%submat(iisub,jjsub,j),iisub=1,ilusch%L%nsubrow),jjsub=1,ilusch%L%nsubcol)
!       end do
!       do j=1,ilusch%U%nnzrow(i)
!          write(288,*) i,ilusch%U%pj(i)%cols(j), & 
!               ((ilusch%U%pa(i)%submat(iisub,jjsub,j),iisub=1,ilusch%U%nsubrow),jjsub=1,ilusch%U%nsubcol)
!       end do
!    end do


!write(6,*)"14 ilutp:"

    if (rmax > 0) then
       deallocate(jw)
       deallocate(w)
       deallocate(jwrev)
       deallocate(iprev)
    end if
write(6,*) "ilutp: There were ",decnt," pivot changes"
    !/*---------------------------------------------------------------------
    !|     done  --  correct return
    !|--------------------------------------------------------------------*/
    return!(0)
994 continue
    !/*  Incomprehensible error. Matrix must be wrong.  */
    return!(1)
    !/*  997 continue
    !   Memory allocation error.  
    !   return!(2)   */
998 continue
    !/*  illegal value for lfil or last entered  */
    return!(5)
9991 continue
    !/*  zero row encountered  */
    return!(6)
    !}
  end subroutine ilutpC
  !/*---------------------------------------------------------------------
  !|     end of ilutpC
  !|--------------------------------------------------------------------*/














  !int ilutD(csptr amat, double *droptol, int *lfil, ilutptr ilusch)
  subroutine ilutD(amat, droptol, lfil, ilusch)
    !{
    !/*---------------------------------------------------------------------- 
    !| ILUT -
    !| Converted to C so that dynamic memory allocation may be implememted.
    !| All indexing is in C format.
    !|----------------------------------------------------------------------
    !| ILUT factorization with dual truncation. 
    !|
    !| March 1, 2000 - dropping in U: keep entries > tau * diagonal entry
    !|----------------------------------------------------------------------
    !|
    !| on entry:
    !|========== 
    !| ( amat ) =  Matrix stored in SparRow struct.
    !| (ilusch) =  Pointer to ILUTfac struct
    !|
    !| lfil[6]  =  number nonzeros in L-part 
    !| lfil[7]  =  number nonzeros in U-part     ( lfil >= 0 )
    !|
    !| droptol[5] = threshold for dropping small terms in L during 
    !|              factorization. 
    !| droptol[6] = threshold for dropping small terms in U.
    !|
    !| On return:
    !|===========
    !|
    !| (ilusch) =  Contains L and U factors in an LUfact struct.
    !|             Individual matrices stored in SparRow structs.
    !|             On return matrices have C (0) indexing.
    !|
    !|       integer value returned:
    !|
    !|             0   --> successful return.
    !|             1   --> Error.  Input matrix may be wrong.  (The 
    !|                         elimination process has generated a
    !|                         row in L or U whose length is > n.)
    !|             2   --> Memory allocation error.
    !|             5   --> Illegal value for lfil or last.
    !|             6   --> zero row encountered.
    !|----------------------------------------------------------------------- 
    !| work arrays:
    !|=============
    !| jw, jwrev = integer work arrays of length n
    !| w         = real work array of length n. 
    !|----------------------------------------------------------------------- 
    !|    
    !|--------------------------------------------------------------------*/
    type(cs_type),pointer ::amat !should already be allocated
    type(ilut_type),pointer::ilusch
    real(dp),dimension(:) :: droptol(7)
    integer,dimension(:) :: lfil 
    
    !local vars
    integer i, j, ii, jj, jcol, jpos, jrow,k
    integer :: iisub,jjsub,kksub,itmp,istat,idiag
    integer,dimension(:),allocatable:: jw, jwrev 
    integer ::len, lenu, lenl, rmax, fil_SL, fil_SU 
    real(dp) drop5,drop6,tnorm,rownorm,wnorm,tmp

    !temporary vars to copy submatrix data
    real(dp),dimension(:),allocatable::d
    real(dp),dimension(:,:),allocatable::Lij,Uij,t,s,fact
    real(dp),dimension(:,:,:),allocatable::w

    integer ::rowz,nsubrow,nsubcol
    type(ja_type),pointer ::rowj
    type(submat_type),pointer::rowm
    
    fil_SL=lfil(6)
    fil_SU=lfil(7)
    nsubrow=amat%nsubrow
    nsubcol=amat%nsubcol
    ilusch%n = amat%n
    rmax = amat%n 
    drop5=droptol(6)
    drop6=droptol(7)
    

  !do i=1,amat%n
  !   do j=1,amat%nnzrow(i)
  !!!      write(85,*) i,amat%pj(i)%cols(j)
  !      write(188,*) i,amat%pj(i)%cols(j),((amat%pa(i)%submat(iisub,jjsub,j),iisub=1,amat%nsubrow),jjsub=1,amat%nsubcol)
  !   end do
  !end do

    if (rmax .eq. 0) then
       return!(0)  
    else !if (rmax > 0) then
       allocate(d(nsubrow),STAT=istat)
       allocate(Lij(nsubrow,nsubcol),STAT=istat)
       allocate(Uij(nsubrow,nsubcol),STAT=istat)
       allocate(t(nsubrow,nsubcol),STAT=istat)
       allocate(s(nsubrow,nsubcol),STAT=istat)
       allocate(fact(nsubrow,nsubcol),STAT=istat)
       allocate(jw(rmax),STAT=istat)
       allocate(w(nsubrow,nsubcol,rmax),STAT=istat)
       allocate(jwrev(rmax),STAT=istat)
    end if

    if (fil_SL<0 .or. fil_SU<0 .or. rmax<=0) goto  9995 

    do j=1,rmax
       jwrev(j) = -1 
    end do

    !/*---------------------------------------------------------------------
    !|    beginning of first main loop - L, U, L^{-1}F calculations
    !|--------------------------------------------------------------------*/
    do ii=1,rmax
!write(6,*) "ilutd: computing row ",ii
       rowj => amat%pj(ii) 
       rowm => amat%pa(ii) 
       rowz =  amat%nnzrow(ii) 
       !/*---------------------------------------------------------------------
       !|   check for zero row 
       !|--------------------------------------------------------------------*/
       tnorm = 0.0 
       do k=1,rowz
          tnorm = tnorm + L1Norm2(rowm%submat,k,nsubrow) !tnorm = tnorm + abs(rowm(k)) 
!write(6,*) "ilutd check for zero row: row ",ii,"tnorm=",tnorm
       end do
       tnorm = tnorm / rowz 

!if(tnorm>10.) write(6,*) "1: row ",ii,"tnorm=",tnorm
!write(6,*) "1: row ",ii,"abs(tnorm-DBL_EPSILON): DBL_EPSILON*tnorm:",abs(tnorm-DBL_EPSILON),DBL_EPSILON*tnorm
       if(abs(tnorm-DBL_EPSILON) <= DBL_EPSILON*tnorm) goto  9996 

       !/*---------------------------------------------------------------------
       !|     unpack amat in arrays w, jw, jwrev
       !|     WE ASSUME THERE IS A DIAGONAL ELEMENT
       !|--------------------------------------------------------------------*/
       lenu = 1
       lenl = 0
       jw(ii) = ii
       jwrev(ii) = ii    
       call blank_3D_submat(w,ii,nsubrow,nsubcol) !w(ii) = 0.0
       do j=1,rowz
          jcol = rowj%cols(j) 
          if (jcol < ii) then 
             lenl=lenl+1 
             jpos=lenl
             jw(jpos) = jcol 
             jwrev(jcol) = lenl 
          else if (jcol .eq. ii) then
             jpos=ii
          else 
             jpos = ii+lenu 
             jw(jpos) = jcol 
             jwrev(jcol) = jpos 
             lenu=lenu+1 
          end if
          do iisub=1,nsubrow
             do jjsub=1,nsubcol
                w(iisub,jjsub,jpos)=rowm%submat(iisub,jjsub,j)
!write(6,*)"w(",ii,jpos,")=",w(iisub,jjsub,jpos)
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
             call swap_submat(w,jj,k,nsubrow)
          end if
          !/*---------------------------------------------------------------------
          !|     zero out element in row.
          !|--------------------------------------------------------------------*/
          jwrev(jrow) = -1 
          !/*---------------------------------------------------------------------
          !|     get the multiplier for row to be eliminated (jrow).
          !|--------------------------------------------------------------------*/
          rowm => ilusch%U%pa(jrow) 
          rowj => ilusch%U%pj(jrow) 
          rowz =  ilusch%U%nnzrow(jrow) 
          !fact = w(jj) * rowm(0) 
          do iisub=1,nsubrow
             do jjsub=1,nsubcol
                fact(iisub,jjsub)=0.0
                do kksub=1,nsubrow !assumes square submatrix
                   fact(iisub,jjsub) = fact(iisub,jjsub) + w(iisub,kksub,jj)*rowm%submat(kksub,jjsub,1)
                end do
             end do
          end do

!if(L1Norm(fact,nsubrow)>10.)write(6,*)"2:ilutp ii:",ii," pivot norm(",jrow,")=",L1Norm(fact,nsubrow)

          if (L1Norm(fact,nsubrow) <= drop5 ) then  !/*  DROPPING IN L  */
!write(6,*)"ilutd dropping norm ",L1Norm(fact,nsubrow)," in L row,column(",ii,",",jrow,")",L1Norm(fact,nsubrow),"<=",drop5
             cycle
          end if

          !/*---------------------------------------------------------------------
          !|     combine current row and row jrow
          !|--------------------------------------------------------------------*/
          do k=2,rowz
             jcol = rowj%cols(k) 
             jpos = jwrev(jcol) 
             if (jpos == -1) then
                !/*---------------------------------------------------------------------
                !|     this is a fill-in element
                !|--------------------------------------------------------------------*/
                if (jcol >= ii) then
                   !/*---------------------------------------------------------------------
                   !|     dealing with U
                   !|--------------------------------------------------------------------*/
                   if(lenu < fil_SU) then
                   !idiag=amat%pd(ii)
                   !if(idiag<=0)idiag=1
                   !if(lenu < (fil_SU+(amat%nnzrow(ii)-idiag+1))) then
                      jpos=ii+lenu
                      jw(jpos) = jcol 
                      jwrev(jcol) = jpos 
                      lenu=lenu+1 
                      call blank_3D_submat(w,jpos,nsubrow,nsubcol)
                   end if
                else 
                   !/*---------------------------------------------------------------------
                   !|     dealing  with L
                   !|--------------------------------------------------------------------*/
                   if(lenl < fil_SL) then
                   !idiag=amat%pd(ii)
                   !if(idiag<=0)idiag=1
                   !if(lenl < (fil_SL+idiag-1)) then
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
                         w(iisub,jjsub,jpos) = w(iisub,jjsub,jpos) - fact(iisub,kksub)*rowm%submat(kksub,jjsub,k)
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

       end do !loop jj

       !/*---------------------------------------------------------------------
       !|     reset nonzero indicators
       !|--------------------------------------------------------------------*/
       do j=1,lenl   !/*  L block  */
          jwrev(jw(j)) = -1 
       end do
       do j=1 ,lenu   !/*  U block  */
          jwrev(jw(ii+j-ofst)) = -1 
       end do

       !/*---------------------------------------------------------------------
       !|     done reducing this row, now store L
       !|--------------------------------------------------------------------*/
       !lenl = len > fil5 ? fil5 : len 
       lenl=len
       if(lenl>fil_SL) lenl=fil_SU
       !idiag=amat%pd(ii)
       !if(idiag<=0)idiag=1
       !if(lenl>(fil_SL+idiag-1)) lenl=fil_SL+(idiag-1)

       ilusch%L%nnzrow(ii) = lenl 
!write(6,*)"U%nnzrow(",ii,")=",ilusch%U%nnzrow(ii)
       !/* quick sort */
       if (len > lenl) call qsplit(w,jw,len,lenl,nsubrow)!quick sort alg (sorts based on L1norm of submats)     
       if (lenl > 0) then 
          allocate(ilusch%L%pj(ii)%cols(lenl),STAT=istat)
          allocate(ilusch%L%pa(ii)%submat(nsubrow,nsubcol,lenl),STAT=istat)
          do k=1,lenl
             ilusch%L%pj(ii)%cols(k)=jw(k)
             do iisub=1,nsubrow
                do jjsub=1,nsubcol
                   ilusch%L%pa(ii)%submat(iisub,jjsub,k)=w(iisub,jjsub,k)!copy w into new U elements, U[k]=w[ii+(k-1)]
                end do
             end do
          end do
       endif

       !/*---------------------------------------------------------------------
       !|     store the diagonal element of U
       !|
       !|     dropping in U if size is less than drop1 * diagonal entry
       !|--------------------------------------------------------------------*/

       tnorm=L1Norm2(w,ii,nsubrow)!ceb

       len = 0 
       do j=2,lenu
          jcol=ii+j-ofst
          if ( L1Norm2(w,jcol,nsubrow) > drop6*tnorm ) then !use matrix norm
             len =len+1
             !jw(len) = jw(jcol) !w(len) = w(ii+j)
             jw(ii+len) = jw(jcol) 
             do iisub=1,nsubrow
                do jjsub=1,nsubcol
                   !w(iisub,jjsub,len)=w(iisub,jjsub,jcol)
                   w(iisub,jjsub,ii+len)= w(iisub,jjsub,jcol) !w(ii+len) = w(ii+j) 
                end do
             end do
          else
!           write(6,*) "ilutd  :U Dropping small column element  norm(w(", jcol, "))=", L1Norm2(w,jcol,nsubrow),&
!                " < tol*tnorm(",ii,")=",drop6*tnorm
          end if
       end do

       !lenu = len+1 > fil6 ? fil6 : len+1 
       lenu=len+1
       len=lenu
       if(lenu > fil_SU)lenu=fil_SU
       !idiag=amat%pd(ii)
       !if(idiag<=0)idiag=1
       !if(lenu > (fil_SU+(amat%nnzrow(ii)-idiag+1)))lenu=fil_SU+(amat%nnzrow(ii)-idiag+1)
       ilusch%U%nnzrow(ii) = lenu
       jpos = lenu 

       if (len > jpos) call qsplit(w,jw,len,jpos,nsubrow)!quick sort alg (sorts based on L1norm of submats)

       if (lenu>0) then
          allocate(ilusch%U%pa(ii)%submat(nsubrow,nsubcol,lenu),STAT=istat)
          allocate(ilusch%U%pj(ii)%cols(lenu),STAT=istat)
       end if


       wnorm=L1Norm2(w,ii,nsubrow)
       if(abs(wnorm-DBL_EPSILON) <= DBL_EPSILON*wnorm) then
          !write(6,*) "ilutD: augmenting diag ",ii
          do iisub=1,nsubrow
             w(iisub,iisub,ii) = (0.0001+drop6)*tnorm 
          end do
       end if




       !/*---------------------------------------------------------------------
       !|     copy the rest of U (all but diagonal)
       !|--------------------------------------------------------------------*/
       !memcpy(&ilusch%U%pj(ii)(1), jw, jpos*sizeof(int)) 
       !/* copy the rest of U */   !ceb do we need to deal with diagonal here? 
       do k=2,jpos! (skips diagonal)
          ilusch%U%pj(ii)%cols(k)=jw(ii+k-ofst)!copy jw into new U_ja elements
          ! ilusch%U%pj(ii)%cols(k)=jw(k-ofst)!copy jw into new U_ja elements
         do iisub=1,nsubrow
             do jjsub=1,nsubcol
                !ilusch%U%pa(ii)%submat(iisub,jjsub,k)=w(iisub,jjsub,k-ofst)!copy w into new U elements, U[k]=w[ii+(k-1)]
                ilusch%U%pa(ii)%submat(iisub,jjsub,k)=w(iisub,jjsub,ii+k-ofst)!copy w into new U elements, U[k]=w[ii+(k-1)]
             end do
          end do
       end do



       !ilusch%U%pa(ii)(0) = 1.0 / t 
       !===========================================================
       !//compute L,U for pivot (diagonal of combined LU matrix)
       ! Lij,Uij submatrices refer only to the diagonal block element
       tmp=0.0
       !diagp=iau(ii)
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
             ilusch%U%pa(ii)%submat(jjsub,kksub,1)=0.0
          end do
          ilusch%U%pa(ii)%submat(jjsub,jjsub,1)=1.0
       end do
       
       !//write LU into U[diagp]
       do kksub=1,nsubcol !loop subcols
          do jjsub=1,nsubcol-1 !get D[1,..,neqn-1] ! loop subcols
             d(jjsub)=ilusch%U%pa(ii)%submat(jjsub,kksub,1) !Identity matrix row jj 
             do iisub=1,jjsub-1
                d(jjsub) = d(jjsub) - Lij(iisub,jjsub)*d(iisub)
             end do
             d(jjsub) = d(jjsub)*Lij(jjsub,jjsub)!ceb not sure if this is right
          end do
          
          do jjsub=nsubcol,1,-1 !get alu[diagp,jj,kk] ! loop subcols
             if(jjsub.eq.nsubcol) then     
                !D=D-L*D (left of subdiag)
                do iisub=1,jjsub-1
                   ilusch%U%pa(ii)%submat(jjsub,kksub,1) = ilusch%U%pa(ii)%submat(jjsub,kksub,1)  - &
                        Lij(iisub,jjsub)*d(iisub)
                end do
                ilusch%U%pa(ii)%submat(jjsub,kksub,1) = ilusch%U%pa(ii)%submat(jjsub,kksub,1)*Lij(jjsub,jjsub)
                !D=
             else      
                !D=D-U*D (right of sub diag)
                ilusch%U%pa(ii)%submat(jjsub,kksub,1)=d(jjsub)
                do iisub=jjsub+1,nsubrow
                   ilusch%U%pa(ii)%submat(jjsub,kksub,1) = ilusch%U%pa(ii)%submat(jjsub,kksub,1) - &
                        Uij(iisub,jjsub)*ilusch%U%pa(ii)%submat(iisub,kksub,1)
                end do
             end if
          end do
       end do
       
       !do kksub=1,nsubcol !loop subcols
       !   do jjsub=nsubcol,1,-1 !get D[1,..,neqn-1] ! loop subcols
       !      write(6,*)"U(",ii,ii,")=",ilusch%U%pa(ii)%submat(jjsub,kksub,1)
       !   end do
       !end do

       !U%ja(ii)%cols(1) = jw(ii) !rowj(0) = jw(ii);
       !===========================================================
       ilusch%U%pj(ii)%cols(1) = ii 

    end do
    !/*---------------------------------------------------------------------
    !|     end main loop - now do clean up
    !|--------------------------------------------------------------------*/
    if (rmax > 0) then
       deallocate(jw) 
       deallocate(w) 
       deallocate(jwrev) 
    end if




    do i=1,ilusch%n
!write(6,*)"ilutd row ",i," L%rownnz(i)=",ilusch%L%nnzrow(i)," U%rownnz(i)=",ilusch%U%nnzrow(i)
       do j=1,ilusch%L%nnzrow(i)
          write(188,*) i,ilusch%L%pj(i)%cols(j), &
               ((ilusch%L%pa(i)%submat(iisub,jjsub,j),iisub=1,ilusch%L%nsubrow),jjsub=1,ilusch%L%nsubcol)
       end do
       do j=1,ilusch%U%nnzrow(i)
          write(189,*) i,ilusch%U%pj(i)%cols(j), & 
               ((ilusch%U%pa(i)%submat(iisub,jjsub,j),iisub=1,ilusch%U%nsubrow),jjsub=1,ilusch%U%nsubcol)
       end do
    end do



    !/*---------------------------------------------------------------------
    !|     done  --  correct return
    !|--------------------------------------------------------------------*/
    return!(0) 
9991 continue
    write(6,*)"  Incomprehensible error. Matrix must be wrong.  "
    return!(1) 
    !/*  9992 continue
    !   Memory allocation error.  
    !   return!(2)   */
9995 continue
    write(6,*)"  illegal value for lfil or last entered  "
    return!(5) 
9996 continue
    write(6,*)"  zero row encountered  "
    return!(6) 
    !}
    !/*---------------------------------------------------------------------
    !|     end of ilutNEW
    !|--------------------------------------------------------------------*/
  end subroutine ilutD
  
end module ilutp
