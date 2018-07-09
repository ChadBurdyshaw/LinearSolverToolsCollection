
!===================================================================
subroutine BLKILUT_SCHUR(dim,neqn,nnz_mat,ia,ja,iau,A,iord,riord,L,U,droptol,fill,schur_start,nnz_pc,ia_pc,ja_pc,iau_pc)
  integer, intent(in) ::fill,dim,neqn,nnz_mat
  integer, intent(inout) ::schur_start
  real(dp), intent(in):: droptol !dropping threshold
  integer, intent(in), dimension(dim+1) :: ia
  integer, intent(in), dimension(nnz_mat) :: ja
  integer, intent(in), dimension(dim) :: iau
  real(dp), intent(in), dimension(neqn,neqn,nnz_mat) :: A
  integer, intent(inout) :: nnz_pc
  integer, intent(inout), allocatable, dimension(:) :: ia_pc
  integer, intent(inout), allocatable, dimension(:) :: ja_pc
  integer, intent(inout), allocatable, dimension(:) :: iau_pc
  !real(dp), intent(inout), allocatable, dimension(:,:,:):: ILU
  integer, intent(in),dimension(dim)::iord,riord
  integer :: ierr,i,j,k,ii,jj,kk,jrow,j1,j2,p,diag,irow
  integer :: nnz_prec,nnz,itmp,rowlength,istat
  integer :: m,n,nsub,iisub,jjsub,kksub,start,jpos,jcol,jstart,jend,istart,diagp,flag
  integer :: lenl,lenu,len
  integer :: offst
  real(dp) :: tmp,tnorm,dt
  real(dp) :: Lij(neqn,neqn)
  real(dp) :: Uij(neqn,neqn)
  real(dp) :: tmat(neqn,neqn)
  real(dp) :: pivot(neqn,neqn)
  real(dp) :: s(neqn,neqn)
  real(dp) :: d(neqn)


  integer,  allocatable,dimension(:)     :: nnz_schur
  integer,  allocatable,dimension(:)     :: jw
  real(dp), allocatable,dimension(:,:,:) :: w

  ! The limitations of Fortran to resize allocated data
  ! requires that we rethink the structure of the L,U data sand indexing
  ! Instead of using ia(nnodes+1) we can use row_nnz(nnodes) to track the number of nonzero nodes per row
  ! Also we can use a pointer to array  ja(nnodes)->cols(row_nnz(irow)) where we can map back to actual column indices
  ! Finally U can be expressed as a pointer to multidimensional array U(nnodes)->submat(neqn,neqn,row_nnz(irow))
  ! this way we can dynamically allocate L,U data one row at a time

  !type :: ja_type
  !   integer, dimension(:), allocatable :: cols ! list of nonzero block cols in row
  !end type ja_type

  !type :: submat_type
  !   real(dp),dimension(:,:,:), allocatable :: submat
  !end type submat_type

  integer,  pointer :: rowj(:)
  real(dp), pointer :: rowm(:,:,:)
  !when L,U have been determined, perhaps these can be combined into our familiar form of ILU and associated ia,ja
  integer :: U_nnodes
  integer,                  allocatable, dimension(:) :: U_rownnz
  type(ja_type),     target,allocatable, dimension(:) :: U_ja
  type(submat_type), target,allocatable, dimension(:) :: U
  integer :: L_nnodes
  integer,                  allocatable, dimension(:) :: L_rownnz
  type(ja_type),     target,allocatable, dimension(:) :: L_ja
  type(submat_type), target,allocatable, dimension(:) :: L
  type(ja_type), target :: ja_row
  integer :: test
  test=0 !
  n=dim
  nsub=neqn*neqn

  !=======================================================
  L_nnodes = n
  allocate(L_rownnz(n),STAT = istat)						
  allocate(L_ja(n),STAT = istat) !just enough space for a diagonal block matrix			
  allocate(L(n),STAT = istat) !just enough space for a diagonal block matrix			

  U_nnodes = n						
  allocate(U_rownnz(n),STAT = istat)						
  allocate(U_ja(n),STAT = istat) !just enough space for a diagonal block matrix			
  allocate(U(n),STAT = istat) !just enough space for a diagonal block matrix			
  
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
  offst=1
  start = 0 !param->start;! offset from normal start of 1					
  if (schur_start == -1) then 
     schur_start = n+1 ! not offset but starting row/col index for schur complement
     m=n
  else
     m=schur_start-1
  endif
  
  !ceb note jw[1..n] maps block column index in w[1..n] to original block column index in A
  !         jw[n+1..2n] maps block (n+) original column index in A to block column index in w[1..n]
  !/* initialize nonzero indicator array */
  !//initialize column indexing array entries to -1
  do jj=1,n
     jw(n+jj)=-1
  end do
  
  dt = droptol!input parameter: dropping tolerance
  ierr = 0
  !m = mat_vcsr->n;!restricted to level of schur complement
  m = n !restricted to level of schur complement ??
  

  !========================================================================
  !/* beginning of main loop */
  do ii=1,m !loop over all included block rows ==============================================================
     !{  
     irow=iord(ii)
     !compute average absval (L1 norm) of elements in row : tnorm
     tnorm = 0.0
     j1 = ia(irow) !j1 = 1
     j2 = ia(irow+1)-1 !j2 = mat_vcsr->nnzrow[ii];!lists number of nonzeros in row ii
     !write(6,*) "j1=",j1," j2=",j2," nnz=",nnz_mat
     ! check to see that row is not singular
     do k=j1,j2 !loop over columns in block row
        tnorm = tnorm + L1Norm2(A,k,neqn,nnz_mat)
     end do    
     ! if row sum is essentially zero, then we have a singular system, exit solver
     if (tnorm < 1.0e-32) then 
        ierr = -3
        write(6,*)"block row ",ii," is singular row, exiting pc solver"
        return !break
     endif
     tnorm = tnorm/real(j2-j1) !divide by number of elements in row*subrow*subcol     
     !write(6,*) "row ",ii,"tnorm=",tnorm
     
     !/* unpack L-part and U-part of row of A in row array w */
     
     !lenu and lenl record length of block row for U and L matrices
     lenu = 0 ! block element index into index jw,row w
     lenl = 0 !lenl does not include diagonal  
  
   if (ii+start < schur_start) then !if current row is less than schur_start row 
        !{
        jw(ii+start) = ii+start !set index to current row
        jw(n+ii+start) = ii+start  ! ?? ceb n+index is work array holding temporary values for L,U row operations       
        call blank_3D_submat(w,ii+start,neqn,neqn,2*n)
        lenu = 1!second col
        do j=j1,j2!loop over columns in current row
           !{
           jcol = riord(ja(j))!ja for orig A  !jcol = mat_vcsr->pj[ii][j];!ja for orig A
           if (jcol < ii+start) then !if current col < diagonal 
              jpos=lenl+offst
              jw(jpos) = jcol !i.e jpos = lenl    
              jw(n+jcol) = jpos !work array index corresponding to L index in w
              lenl=lenl+1 !if block element has been filled, we can increment lenl counter
           else if(jcol == ii+start) then !if current col == diagonal
              jpos=ii+start
           else !if current col > diagonal !this includes elements in W=[L^-1][F]
              jpos = ii + lenu + start 
              jw(jpos)   = jcol
              jw(n+jcol) = jpos 
              lenu=lenu+1!if block element has been filled, we can increment lenu counter
           end if
           do iisub=1,neqn
              do jjsub=1,neqn
                 w(iisub,jjsub,jpos) = A(iisub,jjsub,j) !w[jpos] = Upart[ii,ii+lenu]   = A[ii,jcol];
              end do
           end do
           !}
        end do!end loop over block cols
        ! now w contains all elements of row ii in A
        !}
     else !if current row >= schur_start !storing W portion of factored matrix ===========================
        !{        
        do j=j1,j2 !loop over columns in current row
           !{
           jcol = riord(ja(j)) !jcol = ja(j) for reordered A
           if (jcol < schur_start) then !elements in G
              jpos=lenl+offst
              jw(jpos)   = jcol
              jw(n+jcol) = jpos !work array index points to column in w corresponding to L
              lenl=lenl+1!if block element has been filled, we can increment lenl counter
           else !store cols related to schur complement 
              jpos       =  schur_start + lenu! should this be lenu+1 ?
              jw(jpos)   = jcol
              jw(n+jcol) = jpos !work array index points to column in w corresponding to U
              lenu=lenu+1!if block element has been filled, we can increment lenu counter
           end if
           do iisub=1,neqn
              do jjsub=1,neqn ! loop over subcols              
                 w(iisub,jjsub,jpos) = A(iisub,jjsub,j) !w[jpos]    = A[p,j];
              end do
           end do
           !}
        end do!loop over block cols
        !}
     end if ! end if row>schur_start 
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
           do iisub=1,neqn
              do jjsub=1,neqn
                 s(iisub,jjsub)=w(iisub,jjsub,jj)!s = w[jj];
                 w(iisub,jjsub,jj)=w(iisub,jjsub,k)!w[jj] = w[k];
                 w(iisub,jjsub,k)=s(iisub,jjsub)!w[k] = s;
              end do
           end do
        end if
        jw(n+jrow) = -1 !removes jrow column index from working portion of jw(i>n) !ceb is this necessary
            

        !/* get the multiplier for row to be eliminated (jrow) */
        do iisub=1,neqn
           do jjsub=1,neqn
              pivot(iisub,jjsub)=0.
              do kksub=1,neqn !ceb rowm[1] is diagonal of U here
                 pivot(iisub,jjsub) = pivot(iisub,jjsub) + w(iisub,kksub,jj)*U(jrow)%submat(kksub,jjsub,1)!pivot = w[jj]*rowm[0]; !equivalent L(p,jj)*U(diagj) (pivot)
              end do
           end do
        end do

        !if (A[jj]*U[1]< dropping tolerance) then skip to next row in L
        if (L1Norm(pivot,neqn) <= dt) then !skip to end of loop over rows in L
        !write(6,*)"row ",ii," skip block col ",jrow," (dropping) and go to next"
           cycle!continue
        endif
        
!write(6,*) "checkpt 36"        
       
        !====================================================================================
        !ceb ideally would like to contain both U and L in same ILU data structure
        !with on ia,ja,iau for indexing
     
        !/* combine current row and row jrow */
        !/* if fill-in element is small then disregard */
        if (ii+start < schur_start) then !!/* ilut */========================================================================
           !{
           do k=2,U_rownnz(jrow) !loop over all block column elements in U (exclude diagonal)
              !{
              jcol = U_ja(jrow)%cols(k)  !k column index in row jrow of U
              !check for col in row ii that corresponds to col in row jrow
              !if col in row ii does not exist then make one
              jpos = jw(n+jcol)!jpos work index for column in w corresponding to jcol column in A 
              if (jpos == -1) then ! this is a fill-in element 
                 !{
                 if (jcol >= ii+start) then !block column right of or on diagonal 
                    if(lenu < (fill+(ia(irow+1)-iau(irow)))) then
                       jpos= ii+start+lenu !working element index
                       lenu=lenu+1 !increment block element count for U
                       jw(jpos) = jcol !add new element index to jw
                       jw(n+jcol) = jpos
                       call blank_3D_submat(w,jpos,neqn,neqn,2*n) 
                    end if
                 else !left of diagonal
                    if(lenl < (fill+(iau(irow)-ia(irow)))) then
                       jpos=lenl+offst
                       lenl=lenl+1 !increment block element count for L 
                       jw(jpos) = jcol !add new element index to jw
                       jw(n+jcol) = jpos
                       call blank_3D_submat(w,jpos,neqn,neqn,2*n) 
                       !write(6,*) "row ",ii," adding new column ",jcol," to L.  new lenl=",lenl
                    end if
                 endif
                 !}
              end if
              
              if(jpos.ne.-1) then
                 do iisub=1,neqn
                    do jjsub=1,neqn
                       do kksub=1,neqn
                          w(iisub,jjsub,jpos) = w(iisub,jjsub,jpos) - pivot(iisub,kksub)*U(jrow)%submat(kksub,jjsub,k)
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
           !}
        else !/* partial ilut */ ! apply to W G and S portions of factored matrix ===================================================
           !{
           do k=2,U_rownnz(jrow) !loop over all block column elements in U (exclude diagonal)
              !{
              jcol = U_ja(jrow)%cols(k) 
              jpos = jw(n+jcol)!jpos work index for column in w corresponding to jcol column in A 
              if (jpos == -1) then !/* this is a fill-in element */
                 !{    
                 !ceb no fill requirements to constrain lenl lenu here ?            
                 if (jcol >= schur_start) then !both row and column are > schur start so we are in S matrix
                    !if(lenu < (fill+(ia(irow+1)-iau(irow)))) then
                    !/* dealing with upper part */
                    jpos = schur_start+lenu
                    lenu=lenu+1 !increment block element count for U
                    jw(jpos) = jcol
                    jw(n+jcol) = jpos
                    call blank_3D_submat(w,jpos,neqn,neqn,2*n)
                    !end if 
                 else ! jcol < schur start !we are in G = [E][U^-1] matrix
                    !if(lenl < (fill+(iau(irow)-ia(irow)))) then
                    !/* dealing with lower part */
                    jpos=lenl+offst!ceb
                    lenl=lenl+1 !increment block element count for L
                    jw(jpos) = jcol
                    jw(n+jcol) = jpos
                    call blank_3D_submat(w,jpos,neqn,neqn,2*n)
                    !endif 
                 end if
                 !}
              end if

              if (jpos.ne.-1) then
                 do iisub=1,neqn
                    do jjsub=1,neqn
                       do kksub=1,neqn
                          w(iisub,jjsub,jpos) = w(iisub,jjsub,jpos) - pivot(iisub,kksub)*U(jrow)%submat(kksub,jjsub,k)
                       end do
                    end do
                 end do
              end if
              
              if (lenu > n .or. lenl > n) then ! should not exceed number of possible rows
                 ierr = -1
                 goto 999 !done
              end if
              !}
           end do !end loop over block cols
           !}
        end if  !if (ii+start < schur_start) else !==================================================================
     
        !/* store this pivot element -- from left to right -- no danger
        !	  of overlap with the working elements in L (pivots). */
        do iisub=1,neqn
           do jjsub=1,neqn
              w(iisub,jjsub,len+offst) = pivot(iisub,jjsub) !w[len]  = ft  = L(p,j)
           end do
        end do
                
        jw(len+offst) = jrow
        len=len+1 
        !}
     end do!end loop over block cols: 1..lenl
     !================================================================== 

!write(6,*) "exited loop over block cols: 1..lenl"

     !!/* reset double-pointer to -1 (U-part) */
     itmp=schur_start
     if (ii+start < schur_start) itmp=ii+start
     do k=1,lenu
        jw(n+jw(itmp+(k-offst))) = -1 !ceb may need to decrement this k index by one to match c code
     end do


!ceb If we truncate lenl here then we may be dropping elements from the original ilu0 structure
!    this is not desireable, so we need to sort jw such that the original columns exist in the 
!    new row as well as the largest values from the new columns
!
! note new values in L will be located by jw(jpos) where jpos=lenl >iau(ii)-ia(ii)  and jpos < ii+start
! and new values in U will ve located by jw(jpos) where jpos=ii+start+lenu > ia(ii+1)-iau(ii) and jpos <n 

     !/* update l-mat_vcsr */
     !lenl = len > fill ? fill : len;
     if (len > (fill+(iau(irow)-ia(irow)))) then
        lenl = fill+(iau(irow)-ia(irow))
     else
       lenl = len
     end if
     L_rownnz(ii+start) = lenl


     !/* weigh the elements before sorting */
     do k=1,len
        !{
        do iisub=1,neqn
           do jjsub=1,neqn
              tmat(iisub,jjsub)=w(iisub,jjsub,k) !w[k] here is lower triangular block matrix  ! this is L[ii,k]  !copy w[k]
           end do
        end do
        do iisub=1,neqn
           do jjsub=1,neqn
              w(iisub,jjsub,k) = tmat(iisub,jjsub)*(dt+w(iisub,jjsub,n+jw(k))) !w[k]*=(dt+w[n+jw[k]])
           end do
        end do
        !}
     end do   

     !/* quick sort */
     if (len > lenl) then ! sort row w
        itmp=1
        call qsplitC(w,jw,len+offst,lenl+offst,neqn,itmp)!quick sort alg 
     end if

     !write(6,*)"allocating L( start:",start,",+ii:",ii,") size =",lenl     
     !     allocate(L_ja(start+ii)%cols(lenl+offst),STAT = istat)! create new L_ja elements
     allocate(L_ja(start+ii)%cols(lenl),STAT = istat)! create new L_ja elements
     !     allocate(L(start+ii)%submat(neqn,neqn,lenl+offst),STAT = istat)!create new L elements 
     allocate(L(start+ii)%submat(neqn,neqn,lenl),STAT = istat)!create new L elements 
     
     ! fill L with w for lenl elements
     do k=1,lenl !for (k = 0; k < lenl; k++)
        L_ja(start+ii)%cols(k)=jw(k)!copy jw into new L_ja elements
        do iisub=1,neqn 
           do jjsub=1,neqn 
              !use first def if modified above for weight and sorting
              L(start+ii)%submat(iisub,jjsub,k) = w(iisub,jjsub,k)/(dt+w(iisub,jjsub,n+jw(k)))!rowm[k] = L[k] = w[k] / (dt + w[n+jw[k]]);  ! is this L[k]/(dt+U[k])?
           end do
        end do
     end do
     

     !================================================================== 
     !/* update u-part */
     if (ii+start < schur_start) then
        !{ 
        len = 0      
        do iisub=1,neqn
           do jjsub=1,neqn
              tmat(iisub,jjsub)=dt*w(iisub,jjsub,ii+start)
           end do
        end do
        
        ! for k==1 the following code would just copy w[ii]=w[ii], so we just start at k=2
        do k=2,lenu 
           !if (abs(w(iisub,iisub,ii+k-offst)) > abs(dt*w(iisub,iisub,ii))) then
           if (L1Norm2(w,ii+start+(k-offst),neqn,2*n) > L1Norm(tmat,neqn)) then  !ceb if magnitude of offdiag is greater than diag
              !pre increment len until next subrow is needed
              !write(6,*) "off diag > diag: swapping columns ",ii,"<->",ii+(k-offst)
              len=len+1
              do iisub=1,neqn
                 do jjsub=1,neqn
                    w(iisub,jjsub,(ii+start+len)) = w(iisub,jjsub,ii+start+(k-offst))
                 end do
              end do
              jw(ii+start+len) = jw(ii+start+(k-offst))
           end if
        end do
        
        !lenu = len + 1 > fill ? fill: len + 1;
        if(len+1 > (fill+ia(irow+1)-iau(irow)) ) then
           lenu= fill +(ia(irow+1)-iau(irow))
        else
           lenu=len+1
        end if

        !ceb -----------------------------------------------

        ! note that qsplitC uses L1 norm of submatrix to sort here, when iisub has not cycled through all rows yet
        jpos = (lenu-1)+offst !jpos = lenu - 1
        if (len > (lenu-1)) then   !if (len > (jpos-offst)) then
           itmp=ii+start+1 !ceb does this index need to be incremented by 1 to match c code
           call qsplitC(w,jw,len+offst,jpos,neqn,itmp)!quick sort alg (sorts based on L1norm of submats)
        end if
              
        U_rownnz(ii+start) = lenu   
        !write(6,*)"allocating U(",ii,")_submat size =",lenu     
        allocate(U(ii+start)%submat(neqn,neqn,lenu),STAT = istat)!create new U elements
        allocate(U_ja(ii+start)%cols(lenu),STAT = istat)!create new U_ja elements
        !/* copy the rest of U */   !ceb do we need to deal with diagonal here? 
        do k=2,jpos! (skips diagonal)
           U_ja(ii+start)%cols(k)=jw(ii+start+(k-offst))!copy jw into new U_ja elements
           do iisub=1,neqn
              do jjsub=1,neqn
                 U(ii+start)%submat(iisub,jjsub,k)=w(iisub,jjsub,ii+start+(k-offst))!copy w into new U elements, U[k]=w[ii+(k-1)]
              end do
           end do
        end do
        



        !===========================================================
        !//compute L,U for pivot (diagonal of combined LU matrix)
        ! Lij,Uij submatrices refer only to the diagonal block element
        tmp=0.0
        diagp=iau(irow)
        do jjsub=1,neqn
           do kksub=1,neqn
              if (jjsub.eq.kksub) then!//set diagonal (sub diagonal)
                 tmp=w(jjsub,jjsub,ii+start)!A[p,p](jj,jj)
                 do iisub=1,jjsub-1
                    tmp = tmp - Lij(iisub,jjsub)*Uij(jjsub,iisub) ! A[p,p](jj,jj)=A[p,p](jj,jj)- L[pj](ii,jj)*U[pj](jj,ii)
                 end do
                 Lij(kksub,kksub)= 1.0/tmp 
              else if(kksub<jjsub) then!//set lower
                 tmp = w(jjsub,kksub,ii+start)!A[p,p](jjsub,kksub) L 
                 do iisub=1,kksub-1
                    tmp = tmp - Lij(iisub,jjsub)*Uij(kksub,iisub)
                 end do
                 Lij(kksub,jjsub)=tmp
              else if(kksub>jjsub) then !//set upper
                 tmp = w(jjsub,kksub,ii+start)!A[p,p](jjsub,kksub) U
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
              U(ii)%submat(jjsub,kksub,1)=0.0
           end do
           U(ii)%submat(jjsub,jjsub,1)=1.0
        end do
                
        !//write LU into U[diagp]
        do kksub=1,neqn !loop subcols
           !{  
           do jjsub=1,neqn-1 !get D[1,..,neqn-1] ! loop subcols
              d(jjsub)=U(ii)%submat(jjsub,kksub,1) !Identity matrix row jj 
              do iisub=1,jjsub-1
                 d(jjsub) = d(jjsub) - Lij(iisub,jjsub)*d(iisub)
              end do
              d(jjsub) = d(jjsub)*Lij(jjsub,jjsub)!
           end do
           
           do jjsub=neqn,1,-1 !get alu[diagp,jj,kk] ! loop subcols
              if(jjsub.eq.neqn) then     
                 !Djk=Djk-Lij*Dii (left of subdiag)
                 do iisub=1,jjsub-1
                    U(ii)%submat(jjsub,kksub,1) = U(ii)%submat(jjsub,kksub,1)  - Lij(iisub,jjsub)*d(iisub)
                 end do
                 !Djk=Djk*Ljj
                 U(ii)%submat(jjsub,kksub,1) = U(ii)%submat(jjsub,kksub,1)*Lij(jjsub,jjsub)
              else      
                 !Djk=Djk-Uij*Dik (right of sub diag)
                 U(ii)%submat(jjsub,kksub,1)=d(jjsub)
                 do iisub=jjsub+1,neqn
                    U(ii)%submat(jjsub,kksub,1) = U(ii)%submat(jjsub,kksub,1) - Uij(iisub,jjsub)*U(ii)%submat(iisub,kksub,1)
                 end do
              end if
           end do
           !}
        end do
        
        U_ja(ii+start)%cols(1) = jw(ii+start) !rowj(0) = jw(ii+start);
        
        !===========================================================
        !}
     else if (ii >= schur_start) then !store upper part of schur complement S=A_len+1=[C]-[E][B^-1][F]=[C]-[G][W]
        U_rownnz(ii+start) = lenu
        allocate(U_ja(start+ii)%cols(lenu+offst),STAT = istat)!create new Uja elements
        allocate(U(start+ii)%submat(neqn,neqn,lenu+offst),STAT = istat)!create new U elements     

        do k=1,lenu
           U_ja(start+ii)%cols(k)=jw(schur_start+(k-offst))!copy new Uja elements
           do iisub=1,neqn        
              do jjsub=1,neqn
                 U(start+ii)%submat(iisub,jjsub,k) = w(iisub,jjsub,schur_start+(k-offst)) !copy new U elements
              end do
           end do
        end do
     endif !ii+start < schur_start
     !}     
  end do !loop over included rows
  
  !=================================================================================
  





!  write(6,*)"U:"
!    do i=1,n
!       do iisub=1,neqn
!          do jjsub=1,neqn
!             write(6,*) ( "(",i,U_ja(i)%cols(j),")",U(i)%submat(iisub,jjsub,j),j=1,U_rownnz(i) )
!         end do
!       end do
!    end do
!  
!  write(6,*)"L:"
!    do i=1,n
!       do iisub=1,neqn
!         do jjsub=1,neqn
!            write(6,*) ( "(",i,L_ja(i)%cols(j),")",L(i)%submat(iisub,jjsub,j),j=1,L_rownnz(i) )
!        end do
!      end do
!   end do
  
999 continue !done:
  
  write(6,*)"after done"
  
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
     ja_row%cols(k)=riord(ja(j))
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
  
  
  !===============================================
  !ceb   create schur complement matrix/indices
  !actually we are adding the extra columns and rows that make up our off diagonal blocks in L. i.e. G matrix
  if (flag==1) then
     if (((n+1) - schur_start) > 0) then ! if n>=schur_start
        !{
        !PARMS_NEWARRAY(data->nnzschur, n-schur_start);
        allocate(nnz_schur((n+1)-schur_start),STAT = istat)!number of rows in schur complement
        do i = schur_start,n !loop over included rows
           !{
           lenl = L_rownnz(i)!number of columns in G row i
           k = 1
           do j=1,lenl !loop over length of L row i !loop over columns in G
              if (L_ja(i)%cols(j) < schur_start) then ! if col is in G
                 k=k+1
              endif
           end do
           
           nnz_schur(1+i-schur_start) = k !number of cols in G row i
           
           j1 = 1
           j2 = k
           do j=1,lenl !loop over length of L row i !loop over row columns in G row i
              !{
              if (rowj(j) < schur_start) then !col is in G
                 jcol=j1
                 j1=j1+1
              else !rowj(j) >= schur_start !col is in schur
                 jcol=j2
                 j2=j2+1
              endif
              jw(jcol)  = L_ja(i)%cols(j)
              do iisub=1,neqn
                 do jjsub=1,neqn
                    w(iisub,jjsub,jcol) = L(i)%submat(iisub,jjsub,j) !store L in correct part of working array
                 end do
              end do
              !}
           end do

           do j=1,lenl !loop over length of L row i !loop over row columns in G row i             
              L_ja(i)%cols(j)=jw(j) !assign column index to row i of G
              do iisub=1,neqn
                 do jjsub=1,neqn
                    L(i)%submat(iisub,jjsub,j)=w(iisub,jjsub,j)!only columns < schur_start are assigned back to G     
                 end do
              end do
           end do
           !}
        end do !i=schur_start .. n
        !}  
     end if
  end if
  
  ! deallocate/dereference pointers
  deallocate(w,STAT=istat)  
  deallocate(jw,STAT=istat)
  
  !multiply to see if we get back A
  !  do i=1,n
  !     do iisub=1,neqn
  !        do jjsub=1,neqn
  !           !write(6,*) (U(i)%submat(iisub,jjsub,j),j=1,U_rownnz(i))
  !           !write(6,*) (L(i)%submat(iisub,jjsub,j),j=1,L_rownnz(i))
  !       end do
  !     end do
  !  end do
  
  !cant get rid of schur yet. need to fill ILU with W,G,A_lev+1  
  !deallocate(nnz_schur,STAT=istat)
  !/* compute the number of nonzeros in matrix */
  !ceb this should not have changed from original
  !for (i = 0; i < m; i++) 
  !do i=1,m
  !   !data->nnz_mat += mat_vcsr->nnzrow[i];
  !   nnz_mat += nnzrow[i];
  !end do
  !nnz_mat=ja(m+1) !ceb
  
  
  
  
  
  !/* compute the total number of nonzeros in pc */
  nnz_pc=0
  do i=1,m
     nnz_pc = nnz_pc + L_rownnz(i) + U_rownnz(i)
  end do
  !nnz_pc=nnz_pc+1
  write(6,*) "nnz_pc=",nnz_pc
  !write(6,*) "nnz_mat=",ia(iord(m)+1)-1
  
  !transfer LU data to ILU data structure
  allocate(ILU(neqn,neqn,nnz_pc),STAT=istat)
  allocate(ia_pc(m+1),STAT=istat)
  allocate(iau_pc(m),STAT=istat)
  allocate(ja_pc(nnz_pc),STAT=istat)
  ia_pc(1)=1
  j=0

  do i=1,m !loop over non schur rows
     !add elements from L
     do jcol=1,L_rownnz(i)!+offst
        j=j+1
        ja_pc(j)=L_ja(i)%cols(jcol)
        do iisub=1,neqn
           do jjsub=1,neqn
              ILU(iisub,jjsub,j)=L(i)%submat(iisub,jjsub,jcol)
           end do
        end do
     end do
     !add elements from U     
     do jcol=1,U_rownnz(i)!+offst
        j=j+1
        ja_pc(j)=U_ja(i)%cols(jcol)   
        if(jcol .eq. 1) then
           iau_pc(i)=j
        end if
        do iisub=1,neqn
           do jjsub=1,neqn
              ILU(iisub,jjsub,j)=U(i)%submat(iisub,jjsub,jcol)
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

!     do iisub=1,neqn
!        do jjsub=1,neqn   
!           write(6,*) (ILU(iisub,jjsub,j),j=ia_pc(p),ia_pc(p+1)-1)
!        end do
!     end do

   !  do iisub=1,neqn
   !     write(6,*) (ILU(iisub,iisub,iau_pc(p)))
   !  end do   
!write(6,*)""
!end do 

  !do p=1,nnz_mat
  !   if(ja(p).ne.ja_pc(p)) write(6,*) "ja(",p,")=",ja(p),"ja_pc(",p,")=",ja_pc(p)
  !end do  
  !=========================
  
  deallocate(L,STAT=istat)
  deallocate(U,STAT=istat)
  deallocate(L_ja,STAT=istat)
  deallocate(U_ja,STAT=istat)
  !write(6,*) "ilut chkpt 11.4"
  return! ierr
!=============================================================================     
end subroutine BLKILUT_SCHUR
