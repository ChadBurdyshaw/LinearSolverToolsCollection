module ilut
  use kinddefs!, only: dp,ja_type,submat_type
  use blas
implicit none
contains
!===================================================================
subroutine BLKILUT(dim,neqn,nnz_mat,ia,ja,iau,A,ILU,droptol,fill,nnz_pc,ia_pc,ja_pc,iau_pc)
  integer, intent(in) ::fill,dim,neqn,nnz_mat
  real(dp), intent(in):: droptol !dropping threshold
  integer, intent(in), dimension(dim+1) :: ia
  integer, intent(in), dimension(nnz_mat) :: ja
  integer, intent(in), dimension(dim) :: iau
  real(dp), intent(in), dimension(neqn,neqn,nnz_mat) :: A
  integer, intent(inout) :: nnz_pc
  integer, intent(inout), allocatable, dimension(:) :: ia_pc
  integer, intent(inout), allocatable, dimension(:) :: ja_pc
  integer, intent(inout), allocatable, dimension(:) :: iau_pc
  real(dp), intent(inout), allocatable, dimension(:,:,:):: ILU
  integer :: ierr,i,j,k,ii,jj,kk,jrow,j1,j2,p,diag,irow
  integer :: nnz_prec,nnz,itmp,rowlength,istat
  integer :: m,n,nsub,iisub,jjsub,kksub,jpos,jcol,jstart,jend,istart,diagp,flag
  integer :: lenl,lenu,len,idiag
  real(dp) :: tmp,tnorm,dt,wnorm
  real(dp) :: Lij(neqn,neqn)
  real(dp) :: Uij(neqn,neqn)
  real(dp) :: tmat(neqn,neqn)
  real(dp) :: pivot(neqn,neqn)
  real(dp) :: s(neqn,neqn)
  real(dp) :: d(neqn)
  integer,  allocatable,dimension(:)     :: jw
  real(dp), allocatable,dimension(:,:,:) :: w

  ! The limitations of Fortran to resize allocated data
  ! requires that we rethink the structure of the L,U data sand indexing
  ! Instead of using ia(nnodes+1) we can use row_nnz(nnodes) to track the number of nonzero nodes per row
  ! Also we can use a pointer to array  ja(nnodes)->cols(row_nnz(irow)) where we can map back to actual column indices
  ! Finally U can be expressed as a pointer to multidimensional array U(nnodes)->submat(neqn,neqn,row_nnz(irow))
  ! this way we can dynamically allocate L,U data one row at a time

  !when L,U have been determined, perhaps these can be combined into our familiar form of ILU and associated ia,ja
  type(matrix_type) :: U
  type(matrix_type) :: L
  type(ja_type), target :: ja_row

  integer :: test
  test=0 !
  n=dim
  nsub=neqn*neqn

  !=======================================================
  L%nnodes=n
  allocate(L%rownnz(n),STAT = istat)						
  allocate(L%ja(n),STAT = istat) !just enough space for a diagonal block matrix			
  allocate(L%row(n),STAT = istat) !just enough space for a diagonal block matrix			

  U%nnodes=n					
  allocate(U%rownnz(n),STAT = istat)						
  allocate(U%ja(n),STAT = istat) !just enough space for a diagonal block matrix			
  allocate(U%row(n),STAT = istat) !just enough space for a diagonal block matrix			
  !data->nnz_mat  = 0;
  nnz_prec = 0
  !=======================================================

  !/* malloc memory for working arrays w and jw */
  allocate(w(neqn,neqn,2*n),STAT = istat)!(size 2*n to store entries for bot L and U)
  allocate(jw(2*n),STAT = istat)

  !ceb nnzrow lists the number of nonzero columns in row
  !ceb  this is used instead of using ia to list the index of j columns in the matrix
  !ceb would like to rewrite to use ia as we are accustomed to use in remainder of code
  
  !parms_Map is;! not sure what this is
  !is = self->is;
  

  !do i=1,n
  !   do j=ia(i),ia(i+1)-1
  !      write(18,*) i,ja(j),iau(i)-ia(i)+1,((A(iisub,jjsub,j),iisub=1,neqn),jjsub=1,neqn)
  !   end do
  !end do



  !ceb note jw[1..n] maps block column index in w[1..n] to original block column index in A
  !         jw[n+1..2n] maps block (n+) original column index in A to block column index in w[1..n]
  !/* initialize nonzero indicator array */
  !//initialize column indexing array entries to -1
  do j=1,n
     jw(n+j)=-1
  end do
  
  dt = droptol!input parameter: dropping tolerance
  ierr = 0
  !m = mat_vcsr->n;!restricted to level of schur complement
  m = n!restricted to level of schur complement
  
  !========================================================================
  !/* beginning of main loop */
  do ii=1,m !loop over all included block rows ==============================================================
     !{  
     !compute average absval (L1 norm) of elements in row : tnorm
     tnorm = 0.0
     j1 = ia(ii) !j1 = 1
     j2 = ia(ii+1)-1 !j2 = mat_vcsr->nnzrow[ii];!lists number of nonzeros in row ii
     !write(6,*) "j1=",j1," j2=",j2," nnz=",nnz_mat
     ! check to see that row is not singular
     do k=j1,j2 !loop over columns in block row
        !tnorm = tnorm + L1Norm2(A,k,neqn,nnz_mat)
        tnorm = tnorm + L1Norm2(A,k,neqn)
     end do   
     ! if row sum is essentially zero, then we have a singular system, exit solver
     if (tnorm < 1.0e-32) then 
        ierr = -3
        write(6,*)"block row ",ii," is singular row, exiting pc solver"
        return !break
     endif
     tnorm = tnorm/real(j2-j1) !divide by number of elements in row*subrow*subcol     

     !write(6,*) "ilut: row ",ii,"tnorm=",tnorm
     
     !/* unpack L-part and U-part of row of A in row array w */
     
     !lenu and lenl record length of block row for U and L matrices
     lenl = 0 !lenl does not include diagonal  
     lenu = 1!second col

     jw(ii) = ii !set index to current row
     jw(n+ii) = ii  ! ?? ceb n+index is work array holding temporary values for L,U row operations       
     call blank_3D_submat(w,ii,neqn,neqn)!,2*n)
     do j=j1,j2!loop over columns in current row
        jcol = ja(j)!ja for orig A  !jcol = mat_vcsr->pj[ii][j];!ja for orig A
        if (jcol < ii) then !if current col < diagonal 
           lenl=lenl+1 !if block element has been filled, we can increment lenl counter
           jpos=lenl
           jw(jpos) = jcol !i.e jpos = lenl    
           jw(n+jcol) = jpos !work array index corresponding to L index in w
        else if(jcol == ii) then !if current col == diagonal
           jpos=ii
        else !if current col > diagonal
           jpos = ii + lenu 
           jw(jpos)   = jcol
           jw(n+jcol) = jpos 
           lenu=lenu+1!if block element has been filled, we can increment lenu counter
        end if
        do iisub=1,neqn
           do jjsub=1,neqn
              w(iisub,jjsub,jpos) = A(iisub,jjsub,j) !w[jpos] = Upart[ii,ii+lenu]   = A[ii,jcol];
           end do
        end do
     end do!end loop over block cols
     ! now w contains all elements of row ii in A



          
       !/*---------------------------------------------------------------------
       !|     eliminate previous rows
       !|--------------------------------------------------------------------*/
     !====================================================================
     !ceb need to limit lenl and lenu in the process of adding elements, otherwise row columns can be swapped and original columns
     ! are then lost when the value of lenl or lenu is truncated to match fill level
     len = 0
     !/* eliminate previous rows */
     !loop over block column entries in row ii of L
     !====================================================================
     jj=0
     do !jj=1,lenl !for (jj = 0; jj < lenl; jj++) ! loop over number of block cols in L
        !{
        jj=jj+1
        if (jj.gt.lenl) exit
        !write(6,*) "current lenl:",lenl
        jrow = jw(jj) !note that jrow should be less than diag
        !write(6,*)"ii=",ii," jj=",jj," jrow=",jrow
        !ceb danger here of reordering out the original column structure??
        !  if lenl gets truncated before writing to U,L then yes!
           
        !/* in order to do the elimination in the correct order we must select the smallest column index among jw[k],
        ! k=jj+1,...,lenl.*/
        k = jj
        !loop over columns in jw to look for smallest col index
        do j=jj+1,lenl 
           if (jw(j) < jrow) then
              jrow = jw(j)
              k = j
           end if
        end do
        
        !ceb need to find a faster way to reorder rows or maybe preorder before entering ilut
        
        ! if smaller col index has been found then swap block elements 
        !  with current jj in both jw and w
        if (k .ne. jj) then
           !write(6,*) "reorder L row(",ii,")  swapping col :jw(",jj,"):",jw(jj), " with jw(",k,"):",jw(k)
           !/* exchange in jw */ ! swapping in L index
           j = jw(jj)
           jw(jj) = jw(k)
           jw(k) = j
           !/* exchange in jr */ ! swapping in U index ?
           jw(n+jrow) = jj !this appears to get axed after leaving if statement, so not needed here?
           jw(n+j) = k
           !/* exchange in w */  !swapping occurs in L row
           call swap_submat(w,jj,k,neqn)

           !do iisub=1,neqn
           !   do jjsub=1,neqn
           !      s(iisub,jjsub)=w(iisub,jjsub,jj)!s = w[jj];
           !      w(iisub,jjsub,jj)=w(iisub,jjsub,k)!w[jj] = w[k];
           !      w(iisub,jjsub,k)=s(iisub,jjsub)!w[k] = s;
           !   end do
           !end do
        end if

        jw(n+jrow) = -1 !removes jrow column index from working portion of jw(i>n) !ceb is this necessary
        
        !/* get the multiplier for row to be eliminated (jrow) */
        do iisub=1,neqn
           do jjsub=1,neqn
              pivot(iisub,jjsub)=0.
              do kksub=1,neqn !ceb rowm[1] is diagonal of U here
                 pivot(iisub,jjsub) = pivot(iisub,jjsub) + w(iisub,kksub,jj)*U%row(jrow)%submat(kksub,jjsub,1)!pivot = w[jj]*rowm[0]; !equivalent L(p,jj)*U(diagj) (pivot)
              end do
           end do
        end do

        !if (A[jj]*U[1]< dropping tolerance) then skip to next row in L
        if (L1Norm(pivot,neqn) <= dt) then !skip to end of loop over rows in L
!write(6,*)"ilut dropping norm ",L1Norm(pivot,neqn)," in L row,column(",ii,",",jrow,")",L1Norm(pivot,neqn),"<=",dt
           cycle!continue
        endif
        
        !write(6,*) "checkpt 36"        
       
        !====================================================================================
        !ceb ideally would like to contain both U and L in same ILU data structure
        !with on ia,ja,iau for indexing
     
        !/* combine current row and row jrow */
        !/* if fill-in element is small then disregard */
        do k=2,U%rownnz(jrow) !for (k = 1; k < nnz; k++)  !loop over all block column elements in U (exclude diagonal)
           !{
           jcol = U%ja(jrow)%cols(k)  !column index in A corresponding to k column index in row U(ii)
           !check for col in row ii that corresponds to col in row jrow
           !if col in row ii does not exist then make one
           jpos = jw(n+jcol)!jpos work index for column in w corresponding to jcol column in A 

           if (jpos == -1) then ! this is a fill-in element 
              if (jcol >= ii) then !block column right of or on diagonal
                 idiag=iau(ii)
                 if(idiag==-1)idiag=ia(ii)
                 if(lenu < (fill+(ia(ii+1)-idiag))) then
                    jpos= ii+lenu !working element index
                    jw(jpos) = jcol !add new element index to jw
                    jw(n+jcol) = jpos
                    lenu=lenu+1 !increment block element count for U
                    call blank_3D_submat(w,jpos,neqn,neqn)!,2*n) 
                 end if
              else !left of diagonal
                 idiag=iau(ii)
                 if(idiag==-1)idiag=ia(ii+1)
                 if(lenl < (fill+(idiag-ia(ii)))) then
                    lenl=lenl+1 !increment block element count for L 
                    jpos=lenl
                    jw(jpos) = jcol !add new element index to jw
                    jw(n+jcol) = jpos
                    call blank_3D_submat(w,jpos,neqn,neqn)!,2*n) 
                    !write(6,*) "row ",ii," adding new column ",jcol," to L.  new lenl=",lenl
                 end if
              endif
           end if
           
           if(jpos.ne.-1) then
              do iisub=1,neqn
                 do jjsub=1,neqn
                    do kksub=1,neqn
                       w(iisub,jjsub,jpos) = w(iisub,jjsub,jpos) - pivot(iisub,kksub)*U%row(jrow)%submat(kksub,jjsub,k)
                    end do
                 end do
              end do
           end if
           
           if (lenu > n .or. lenl > n) then ! should not exceed number of possible rows
              ierr = -1
              goto 999 !done
           end if
           !}
        end do!end loop over block cols U:  k=1,U_rownnz(jrow)
        

        !/* store this pivot element -- from left to right -- no danger
        !	  of overlap with the working elements in L (pivots). */
        len=len+1 
        jw(len) = jrow
        do iisub=1,neqn
           do jjsub=1,neqn
              w(iisub,jjsub,len) = pivot(iisub,jjsub) !w[len]  = ft  = L(p,j)
           end do
        end do
        !}
     end do!end loop over block cols: 1..lenl
     !================================================================== 
     
     !write(6,*) "exited loop over block cols: 1..lenl"
     
     !!/* reset double-pointer to -1 (U-part) */
     !itmp=ii!schur_start
     !if (ii < schur_start) itmp=ii
     do k=1,lenu
        jw(n+jw(ii+(k-ofst))) = -1 !ceb may need to decrement this k index by one to match c code
     end do
          
     !ceb If we truncate lenl here then we may be dropping elements from the original ilu0 structure
     !    this is not desireable, so we need to sort jw such that the original columns exist in the 
     !    new row as well as the largest values from the new columns
     !
     ! note new values in L will be located by jw(jpos) where jpos=lenl >iau(ii)-ia(ii)  and jpos < ii
     ! and new values in U will ve located by jw(jpos) where jpos=ii+lenu > ia(ii+1)-iau(ii) and jpos <n 
     
     !/* update l-mat_vcsr */
     !lenl = len > fill ? fill : len;
     lenl=len
     idiag=iau(ii)
     if(idiag==-1)idiag=ia(ii+1)
     if (lenl > (fill+(idiag-ia(ii)))) lenl = fill+(idiag-ia(ii))
     !if (lenl > fill) lenl = fill
     L%rownnz(ii) = lenl
     
     
     !ceb
     !/* weigh the elements before sorting */
     !do k=1,len !for (k = 0; k < len; k++) 
     !   do iisub=1,neqn
     !      do jjsub=1,neqn
     !         tmat(iisub,jjsub)=w(iisub,jjsub,k) !w[k] here is lower triangular block matrix  ! this is L[ii,k]  !copy w[k]
     !      end do
     !   end do
     !   do iisub=1,neqn
     !      do jjsub=1,neqn
     !         w(iisub,jjsub,k) = tmat(iisub,jjsub)*(dt+w(iisub,jjsub,n+jw(k))) !w[k]*=(dt+w[n+jw[k]])
     !      end do
     !   end do
     !end do
     !ceb

     !/* quick sort */
     if (len > lenl) call qsplit(w,jw,len+ofst,lenl+ofst,neqn)!quick sort alg 
     
     if(lenl>0) then
        allocate(L%ja(ii)%cols(lenl),STAT = istat)! create new L_ja elements
        allocate(L%row(ii)%submat(neqn,neqn,lenl),STAT = istat)!create new L elements 
        
        ! fill L with w for lenl elements
        do k=1,lenl !for (k = 0; k < lenl; k++)
           L%ja(ii)%cols(k)=jw(k)!copy jw into new L_ja elements
           do iisub=1,neqn 
              do jjsub=1,neqn 
                 !use first def if modified above for weight and sorting
                 !ceb L%row(ii)%submat(iisub,jjsub,k) = w(iisub,jjsub,k)/(dt+w(iisub,jjsub,n+jw(k)))!rowm[k] = L[k] = w[k] / (dt + w[n+jw[k]]);  ! is this L[k]/(dt+U[k])?
                 L%row(ii)%submat(iisub,jjsub,k) = w(iisub,jjsub,k)!/(dt+w(iisub,jjsub,n+jw(k)))!rowm[k] = L[k] = w[k] / (dt + w[n+jw[k]]);  ! is this L[k]/(dt+U[k])?
              end do
           end do
        end do
     end if


     
     !================================================================== 
     !/* update u-part */
     !do iisub=1,neqn
     !   do jjsub=1,neqn
     !      tmat(iisub,jjsub)=dt*w(iisub,jjsub,ii)
     !   end do
     !end do

     !write(6,*) "ilut: update U: tnorm(",ii,")= ",L1Norm(tmat,neqn)

     !/*---------------------------------------------------------------------
     !|     apply dropping strategy to U  (entries after diagonal)
     !|--------------------------------------------------------------------*/

     ! for k==1 the following code would just copy w[ii]=w[ii], so we just start at k=2
     wnorm=L1Norm2(w,ii,neqn)
     len = 0      
     do k=2,lenu !for (k = 1; k < lenu; k++) 
        jcol=ii+k-ofst
        !if (abs(w(iisub,iisub,ii+k-ofst)) > abs(dt*w(iisub,iisub,ii))) then
        !if (L1Norm2(w,ii+(k-ofst),neqn,2*n) > L1Norm(tmat,neqn)) then  !ceb if magnitude of offdiag is greater than diag
        if (L1Norm2(w,jcol,neqn) > dt*wnorm) then  !ceb if magnitude of offdiag is greater than diag
           !pre increment len until next subrow is needed
!write(6,*) "ilut: off diag > diag: swapping columns ",ii,"<->",jcol
           len=len+1
           jw(ii+len) = jw(jcol)
           do iisub=1,neqn
              do jjsub=1,neqn
                 w(iisub,jjsub,(ii+len)) = w(iisub,jjsub,jcol)
              end do
           end do
        else
!           write(6,*) "ilut :U Dropping small column element norm(w(", jcol, "))=", L1Norm2(w,jcol,neqn),&
!                " < dt*wnorm(",ii,")=", dt*wnorm     
        end if
     end do
     
     
     lenu=len+1 !lenu = len + 1 > fill ? fill: len + 1;
     idiag=iau(ii)
     if(idiag==-1)idiag=ia(ii)
     if(lenu > (fill+ia(ii+1)-idiag) ) lenu= fill +(ia(ii+1)-idiag)
     !if(lenu > fill) lenu= fill
     U%rownnz(ii) = lenu   
     jpos = lenu !jpos = lenu - 1
     
     !ceb -----------------------------------------------    
     !if (test.eq.1) then ! reorder U row. Not required but seems to improve convergence sometimes???
     !   !ceb order upper triangular part of row if necessary
     !   do jj=1,lenu ! loop over number of block cols in U
     !      jrow = jw((ii)+(jj-1))
     !      k = jj
     !      !loop over columns in jw to look for smaller col index
     !      do j=jj+1,lenu 
     !         if (jw((ii)+(j-1)) < jrow) then
     !            jrow = jw((ii)+(j-1))
     !            k = j
     !         end if
     !      end do
     !      if (k .ne. jj) then
     !         !write(6,*) "reorder U row(",ii,")  swapping col :jw(",jj,"):",jw(ii+(jj-1)), " with jw(",k,"):",jw(ii+(k-1))
     !         !/* exchange in jw */ ! swapping in L index
     !         j = jw((ii)+(jj-1))
     !         jw((ii)+(jj-1)) = jw((ii)+(k-1))
     !         jw((ii)+(k-1)) = j
     !         !/* exchange in jr */ !if these two lines are not removed we end up skipping(dropping) some cols
     !         !jw(n+jrow) = (jj)
     !         !jw(n+j)    = (k)
     !         !/* exchange in w */  !swapping occurs in L row
     !         do iisub=1,neqn
     !            do jjsub=1,neqn
     !               s(iisub,jjsub)=w(iisub,jjsub,(ii)+(jj-1))!s = w[jj];
     !               w(iisub,jjsub,(ii)+(jj-1))=w(iisub,jjsub,(ii)+(k-1))!w[jj] = w[k];
     !               w(iisub,jjsub,(ii)+(k-1))=s(iisub,jjsub)!w[k] = s;
     !            end do
     !         end do
     !      end if
     !   end do
     !   !------------------------------------------------ceb
     !end if
     
     
     ! note that qsplitC uses L1 norm of submatrix to sort here, when iisub has not cycled through all rows yet
     if (jpos < len+ofst) then
        !qsplitC(&w[ii+1], &jw[ii+1], neqn,len, jpos);!quick sort alg
        itmp=ii !ceb does this index need to be incremented by 1 to match c code
        call qsplitC(w,jw,len+ofst,jpos,neqn,itmp)!quick sort alg (sorts based on L1norm of submats)
     end if
     
     if (lenu>0) then
        !write(6,*)"allocating U(",ii,")_submat size =",lenu     
        allocate(U%row(ii)%submat(neqn,neqn,lenu),STAT = istat)!create new U elements
        allocate(U%ja(ii)%cols(lenu),STAT = istat)!create new U_ja elements
     end if

     !/* copy the rest of U */   !ceb do we need to deal with diagonal here? 
     do k=2,jpos! (skips diagonal)
        U%ja(ii)%cols(k)=jw(ii+(k-ofst))!copy jw into new U_ja elements
        do iisub=1,neqn
           do jjsub=1,neqn
              U%row(ii)%submat(iisub,jjsub,k)=w(iisub,jjsub,ii+(k-ofst))!copy w into new U elements, U[k]=w[ii+(k-1)]
           end do
        end do
     end do


     !===========================================================
     !//compute L,U for pivot (diagonal of combined LU matrix)
     ! Lij,Uij submatrices refer only to the diagonal block element
     tmp=0.0
     !diagp=iau(ii)
     do jjsub=1,neqn
        do kksub=1,neqn
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
     do jjsub=1,neqn
        do kksub=1,neqn
           U%row(ii)%submat(jjsub,kksub,1)=0.0
        end do
        U%row(ii)%submat(jjsub,jjsub,1)=1.0
     end do
     
     !//write LU into U[diagp]
     do kksub=1,neqn !loop subcols
        do jjsub=1,neqn-1 !get D[1,..,neqn-1] ! loop subcols
           d(jjsub)=U%row(ii)%submat(jjsub,kksub,1) !Identity matrix row jj 
           do iisub=1,jjsub-1
              d(jjsub) = d(jjsub) - Lij(iisub,jjsub)*d(iisub)
           end do
           d(jjsub) = d(jjsub)*Lij(jjsub,jjsub)!ceb not sure if this is right
        end do
        
        do jjsub=neqn,1,-1 !get alu[diagp,jj,kk] ! loop subcols
           if(jjsub.eq.neqn) then     
              !D=D-L*D (left of subdiag)
              do iisub=1,jjsub-1
                 U%row(ii)%submat(jjsub,kksub,1) = U%row(ii)%submat(jjsub,kksub,1)  - Lij(iisub,jjsub)*d(iisub)
              end do
              U%row(ii)%submat(jjsub,kksub,1) = U%row(ii)%submat(jjsub,kksub,1)*Lij(jjsub,jjsub)
              !D=
           else      
              !D=D-U*D (right of sub diag)
              U%row(ii)%submat(jjsub,kksub,1)=d(jjsub)
              do iisub=jjsub+1,neqn
                 U%row(ii)%submat(jjsub,kksub,1) = U%row(ii)%submat(jjsub,kksub,1) - &
                      Uij(iisub,jjsub)*U%row(ii)%submat(iisub,kksub,1)
              end do
           end if
        end do
     end do
     
     U%ja(ii)%cols(1) = jw(ii) !rowj(0) = jw(ii);
     
     !===========================================================
     
     
     !}     
  end do !loop over included rows
  
  !=================================================================================
  
  
  
  
  
  
  
  !  write(6,*)"U:"
  !    do i=1,n
  !       do iisub=1,neqn
  !          do jjsub=1,neqn
  !             write(6,*) ( "(",i,U%ja(i)%cols(j),")",U%row(i)%submat(iisub,jjsub,j),j=1,U%rownnz(i) )
  !         end do
  !       end do
  !    end do
  !  
  !  write(6,*)"L:"
  !    do i=1,n
  !       do iisub=1,neqn
  !         do jjsub=1,neqn
  !            write(6,*) ( "(",i,L%ja(i)%cols(j),")",L%row(i)%submat(iisub,jjsub,j),j=1,L%rownnz(i) )
  !        end do
  !      end do
  !   end do
  
  
  
  
999 continue !done:
  
  !write(6,*)"after done"
  
  !data->schur_start = is->schur_start
  !data->n = n
  
  flag = 0
  !nnz   = mat_vcsr->nnzrow[m-1];
  !nnz   = ia(m)-1;!ceb may not be right (if start is not ia[1])
  !copy ja into ja_row
  rowlength=ia(m)-ia(m-1)
  allocate(ja_row%cols(rowlength),STAT = istat)
  k=1
  do j=ia(m-1),ia(m)-1
     ja_row%cols(k)=ja(j)
     k = k+1
  end do
  !rowj  = mat_vcsr->pj[m-1];!ja for orig A
  !for (i = 0; i < nnz; i++) 
  do i=1,rowlength
     !if (rowj(i) == n-1) then
     if (ja_row%cols(i) == n) then
        flag = 1
        exit !break;
     end if
  end do
  deallocate(ja_row%cols,STAT = istat)
  
  

  
  ! deallocate/dereference pointers
  deallocate(w,STAT=istat)  
  deallocate(jw,STAT=istat)
  
  !multiply to see if we get back A
  !  do i=1,n
  !     do iisub=1,neqn
  !        do jjsub=1,neqn
  !           !write(6,*) (L%row(i)%submat(iisub,jjsub,j),j=1,L_rownnz(i))
  !       end do
  !     end do
  !  end do
  
  
  !/* compute the total number of nonzeros in pc */
  nnz_pc=0
  do i=1,m
!write(6,*)"ilut row ",i," L%rownnz(i)=",L%rownnz(i)," U%rownnz(i)=",U%rownnz(i)
     nnz_pc = nnz_pc + L%rownnz(i) + U%rownnz(i)
  end do
  !nnz_pc=nnz_pc+1
  write(6,*) "nnz_pc=",nnz_pc
  write(6,*) "nnz_mat=",ia(m+1)-1
  
  !transfer LU data to ILU data structure
  allocate(ILU(neqn,neqn,nnz_pc),STAT=istat)
  allocate(ia_pc(m+1),STAT=istat)
  allocate(iau_pc(m),STAT=istat)
  allocate(ja_pc(nnz_pc),STAT=istat)
  ia_pc(1)=1
  j=0

  do i=1,m !loop over non schur rows
     !add elements from L
     do jcol=1,L%rownnz(i)!+ofst
        j=j+1
        ja_pc(j)=L%ja(i)%cols(jcol)
        do iisub=1,neqn
           do jjsub=1,neqn
              ILU(iisub,jjsub,j)=L%row(i)%submat(iisub,jjsub,jcol)
           end do
        end do
     end do
     !add elements from U     
     do jcol=1,U%rownnz(i)!+ofst
        j=j+1
        ja_pc(j)=U%ja(i)%cols(jcol)   
        if(jcol .eq. 1) then
           iau_pc(i)=j
        end if
        do iisub=1,neqn
           do jjsub=1,neqn
              ILU(iisub,jjsub,j)=U%row(i)%submat(iisub,jjsub,jcol)
           end do
        end do
     end do
     ia_pc(i+1)=j+1
  end do!end loop over non schur rows
  
  
  
  !do p=1,m
  !write(6,*)"row:",p
  !     do j=ia_pc(p),ia_pc(p+1)-1
  !        !if(ja(j).ne.ja_pc(j)) write(6,*)"ja(",j,"):",ja(j),"!=ja_pc(",j,"):",ja_pc(j)
  !        write(6,*)"ja_pc(",j,"):",ja_pc(j)
  !     end do
  !end do 
  
  
  !ceb test output ========
!  write(6,*)""
  !  write(6,*)"ilu:"
  !do p=1,m
  !if(iau(p).ne.iau_pc(p)) write(6,*) "iau(",p,")=",iau(p),"iau_pc(",p,")=",iau_pc(p)
  !if(ia(p).ne.ia_pc(p)) 
  !write(6,*) "ia(",p,")=",ia(p),"ia_pc(",p,")=",ia_pc(p)
  !write(6,*) "ja_pc(iau_pc(",p,"))=",ja_pc(iau_pc(p))
  !write(6,*)"row:",p
  !do j=ia_pc(p),ia_pc(p+1)-1
  !   if(ja(j).ne.ja_pc(j)) write(6,*)"ja(",j,"):",ja(j),"!=ja_pc(",j,"):",ja_pc(j)
  !   write(6,*)"ja_pc(",j,"):",ja_pc(j)
  !end do
  
       !do iisub=1,neqn
       !   do jjsub=1,neqn   
       !      write(6,*) (ILU(iisub,jjsub,j),j=ia_pc(p),ia_pc(p+1)-1)
       !   end do
       !end do
  
  !  do iisub=1,neqn
  !     write(6,*) (ILU(iisub,iisub,iau_pc(p)))
  !  end do   
  !write(6,*)""
  !end do 
  
  !do p=1,nnz_mat
  !   if(ja(p).ne.ja_pc(p)) write(6,*) "ja(",p,")=",ja(p),"ja_pc(",p,")=",ja_pc(p)
  !end do  
  !=========================
  
  deallocate(L%row,STAT=istat)
  deallocate(U%row,STAT=istat)
  deallocate(L%ja,STAT=istat)
  deallocate(U%ja,STAT=istat)
  deallocate(L%rownnz,STAT=istat)
  deallocate(U%rownnz,STAT=istat)

  !write(6,*) "ilut chkpt 11.4"
  return! ierr
  !=============================================================================     
end subroutine BLKILUT

end module ilut
