module armslib
  use kinddefs !, only : dp
  use blas
  use matsol_lib
  !use ilut
  !use pilu_lib
implicit none
contains


subroutine indepSetReorder(A,ia,ja,n,nsub,nnz,bsize,iord,riord,nnod,tol,nbnd) 
  !/*--------------------------------------------------------------------- 
  !| greedy algorithm for independent set ordering -- 
  !|----------------------------------------------------------------------
  !|     Input parameters:
  !|     -----------------
  !|     (mat)  =  matrix in SparRow format
  !|     
  !|     bsize  =  integer (input) the target size of each block.
  !|               each block is of size >= bsize. 
  !|
  !|     w      =  weight factors for the selection of the elements in the
  !|               independent set. If w(i) is small i will be left for the
  !|               vertex cover set. 
  !|
  !|     tol    =  a tolerance for excluding a row from independent set.
  !|
  !|     nbnd   =  number of interior variables.
  !|
  !|     Output parameters:
  !|     ------------------ 
  !|     iord   = permutation array corresponding to the independent set 
  !|     ordering.  Row number i will become row number iord[i] in 
  !|     permuted matrix.
  !|     riord   = reversed ordering of iord
  !|     
  !|     nnod   = (output) number of elements in the independent set.
  !|     
  !|----------------------------------------------------------------------- 
  !|     the algorithm searches nodes in lexicographic order and groups
  !|     the (BSIZE-1) nearest nodes of the current to form a block of
  !|     size BSIZE. The current algorithm does not use values of the matrix.
  !|---------------------------------------------------------------------*/ 
  integer, intent(in) :: n,nsub,nnz,bsize,nbnd
  real(dp),intent(in) :: A(nsub,nsub,nnz)
  integer,intent(in) ::ia(n+1),ja(nnz)
  real(dp),intent(in)::tol
  
  integer,intent(inout)::nnod
  integer,intent(inout)::iord(n)
  integer,intent(inout)::riord(n)
  !/*   local variables   */
  integer :: i,jstart,jend,nod, jcount, lastlev, begin, last0, last, nback, mid
  integer :: j1, j2, jcol, inod, jnod, j, k, jcount0, begin0
  integer prog,istat
  real(dp),allocatable,dimension(:) :: w
  integer,allocatable,dimension(:)::iaT
  integer,allocatable,dimension(:)::jaT
  !allocate memory
  !/*-----------------------------------------------------------------------*/
  allocate(w(n),STAT=istat)
  allocate(iaT(n+1),STAT=istat)
  allocate(jaT(nnz),STAT=istat)
  !/*  	 call weights to compute the weights for  input matrix.. */
  
  !ceb  create iaT,jaT  index for matT
  !ceb a transpose of mat is needed for case where mat is nonsymmetric
  call SparTran_mod(ia,ja,iaT,jaT,n,nnz) !transpose matrix mat and store in matT
  call weightsC_mod(A,ia,ja,n,nsub,nnz,w)
  
  !/*---------------------------------------------------------------------- 
  !| scan all nodes first to eliminate those not satisfying DD criterion 
  !+----------------------------------------------------------------------*/
  !nback = n-1
  nback = n
  !nod = 1
  
  do j=1, n
     iord(j) = -1
  end do
  
  !add all rows with index greater than nbnd to complement set
  do j = n, nbnd+1, -1
     call add2com(nback, j, iord, riord)
  end do
  
  !add rows with diagonal dominance less than tolerance to complement set
  do j=1, nbnd
     if (w(j) < tol) then
        call add2com(nback, j, iord, riord)
        !nod=nod+1!not sure why this is incremented since its not used here and gets reset below
     end if
  end do
  
  !last = -1
  last = 1
  nod=0
  !for (nod=0; nod<n; nod++)
  !===================================================================
  do 
     nod=nod+1
     if(nod>n) exit
     !while (iord[nod] .ne. -1)   
     do !loop through nodes to find row that has not yet been added to a set
        !if row(nod) has not been added to a set
        if (iord(nod) .eq. -1) exit
        nod=nod+1
        !if (nod >= n) goto 150 !exit and begin reorder of nodes
        if (nod > n) goto 150 !exit and begin reorder of nodes
     end do
     
     !/*-------------------- initialize level-set - contains nod (only)*/
     call add2is(last, nod, iord, riord)
     begin   = last
     begin0  = begin
     lastlev = begin
     jcount  = 1
          
     !/*----------------------------------------------------------------------
     !|     put all the nearest neighbor nodes of the current node into
     !|     the block until the number is BSIZE.
     !|---------------------------------------------------------------------*/
     prog = 1
     !while (jcount < bsize && prog) 
     do
        if((jcount >= bsize) .or. (prog.eq.0)) then
           exit
        end if
        !ceb jnode may be a bit larger than bsize since it is accumulated
        !    below during traversal of the row/colum of mat/matT and does not
        !    check against bsize intil row/col traversal has finished
        !/*--------------------   traverse all the current level-set   */
        last0 = last
        jcount0 = jcount
        do inod=begin,last0 
           jnod = riord(inod)
           !/*--------------------   This assumes A is not symmetric.   */
           !set up so that this operation is performed once for mat and once for matT    
           !mat
           !loop over columns in row jnod
           do j=ia(jnod),ia(jnod+1)-1
              jcol = ja(j)
              if (iord(jcol) == -1) then
                 !call add2is(last, jcol, iord, riord,n)
                 call add2is(last, jcol, iord, riord)
                 jcount=jcount+1
              end if
           end do
           !matT
           !loop over rows in column jnod
           do j=iaT(jnod),iaT(jnod+1)-1
              jcol = jaT(j)
              if (iord(jcol) == -1) then
                 !call add2is(last, jcol, iord, riord,n)
                 call add2is(last, jcol, iord, riord)
                 jcount=jcount+1
              end if
           end do
        end do
        !prog = jcount > jcount0 ? 1 : 0
        prog=0
        if(jcount>jcount0) prog=1
        lastlev = begin
        begin = last0+1
     end do
          
     !/*-----------------------------------------------------------------------
     !| the neighbors of elements of last level go to the complement   
     !| gmat loop over original matrix and its transpose 
     !+-----------------------------------------------------------------------*/ 
     !gmat = mat
     !do k=1,2 !set up so that this operation is performed once for mat and once for matT

     !mat
     do inod=lastlev,last  
        jnod = riord(inod)
        do j=ia(jnod),ia(jnod+1)-1
           jcol = ja(j)
           if (iord(jcol) == -1) then
              !call add2com(nback, jcol, iord, riord,n)
              call add2com(nback, jcol, iord, riord)
           end if
        end do
     end do
     !matT
     do inod=lastlev,last  
        jnod = riord(inod)
        do j=iaT(jnod),iaT(jnod+1)-1
           jcol = jaT(j)
           if (iord(jcol) == -1) then
              !call add2com(nback, jcol, iord, riord,n)
              call add2com(nback, jcol, iord, riord)
           end if
        end do
     end do
     !/*   reverse ordering for this level   */
     mid = (begin0+last) / 2
     do inod=begin0, mid
        j = last - inod + begin0
        jnod = riord(inod)
        riord(inod) = riord(j)
        riord(j) = jnod
     end do
     
  end do !nod=1, n   !/*--------------------------------------------------
    
  !/*-------------------- relabel nodes of vertex cover   */
150 continue!  label50:
    
  nnod = last
  j1 = nnod
  do j2=nnod+1,n 
     if (iord(riord(j2)) .ne. -1) then
        j1=j1+1
        if (j1 .ne. j2) then 
           j = riord(j2)
           riord(j2) = riord(j1)
           riord(j1) = j
        end if
     end if
  end do
  
  !/*-------------------- obtain reverse permutation array   */
  do j=1,n
     iord(riord(j)) = j
  end do
  
  !ceb output test
  !do j=1,n  
  !   write(6,*)"iord(",j,")=",iord(j)
  !end do
  
  !nnod=nnod+1!=last+1  why inc nnode here?
  
  write(6,*)"indepSetReorder: nnod=",nnod
  
  !do k = 1,n
  !   i=iord(k)
  !   jstart = ia(i)
  !   jend   = ia(i+1) - 1
  !   do j = jstart,jend
  !      jcol = iord(ja(j))
  !      write(46,'(2(i5,1x))')i,jcol
  !  end do
  !end do
  
  !write(6,*) "inpedsetReorder chkpt 8"
  !cleanCS(matT) 
  !if(allocated(riord)) deallocate(riord,STAT=istat)
  !write(6,*) "inpedsetReorder chkpt 9"
  if(allocated(w))deallocate(w,STAT=istat)
  !write(6,*) "inpedsetReorder chkpt 10"
  if(allocated(iaT))deallocate(iaT,STAT=istat)
  !write(6,*) "inpedsetReorder chkpt 11"
  if(allocated(jaT))deallocate(jaT,STAT=istat)
  !write(6,*) "inpedsetReorder chkpt 12"
  return
  !}=====================================================================
end subroutine indepSetReorder





!int PQperm(csptr mat, int bsize, int *Pord, int *Qord, int *nnod, double tol, int nbnd) 
!ceb greedy matching set algorithm
subroutine PQperm( mat,  bsize, Pord, Qord, nnod,  tol,  nbnd) 
  !{
  !/*--------------------------------------------------------------------- 
  !| algorithm for nonsymmetric  block selection - 
  !|----------------------------------------------------------------------
  !|     Input parameters:
  !|     -----------------
  !|     (mat)  =  matrix in SparRow format
  !|     
  !|     tol    =  a tolerance for excluding a row from B block 
  !|
  !|     bsize not used here - it is used in arms2..
  !|
  !|     nbnd   =  number of interior variables.
  !|
  !|     Output parameters:
  !|     ------------------ 
  !|     Pord   = row permutation array.  Row number i will become row number 
  !|		Pord[i] in permuted matrix. (Old to new labels) 
  !|     Qord   = column permutation array.  Column number j will become 
  !|		row number Qord[j] in permuted matrix.
  !|             [destination lists] - i.e., old to new arrays= 
  !|     
  !|     nnod   = number of elements in the B-block 
  !|     
  !|---------------------------------------------------------------------*/ 
  type(cs_type),pointer :: mat
  integer,intent(in) :: bsize,nbnd
  integer,intent(inout)::nnod
  integer,intent(inout),pointer,dimension(:) :: Pord,Qord
  real(dp),intent(in) :: tol
  !/*--------------------   local variables   */
  integer :: i,j,ii,k,col,jj,rnz,nzi,n,count,numnode
  integer, dimension(:),allocatable :: icor,jcor
  real(dp) :: aij,rn
  type(ja_type), pointer  :: row
  type(submat_type), pointer  :: mrow
  integer :: istat,nsubrow,itmp,job
  !/*-----------------------------------------------------------------------*/  
write(6,*)"calling pqperm"

  nsubrow=mat%nsubrow
  n=mat%n
  do j=1,n
     Pord(j) = -1
     Qord(j) = -1 
  end do

  allocate(icor(n),STAT=istat)
  allocate(jcor(n),STAT=istat)
  numnode = 0
  count = 0

  !do i=1,mat%n
  !   do j=1,mat%nnzrow(i)
  !!      write(74,*) i, mat%pj(i)%cols(j)
  !      write(74,*) i, mat%pj(i)%cols(j),((mat%pa(i)%submat(ii,jj,j),ii=1,mat%nsubrow),jj=1,mat%nsubcol)
  !   end do
  !end do

  !/*-------------------- wDiag selects candidate entries in a sorted oder */
  job = 1 !allow swapping of matrix entries

  !preselect ordered set based on diagonal dominance weight calculation
  !order by decreasing weight (best to worst)
  call preSel(mat, icor, jcor, job, tol, count, nbnd) !filter out poor rows

write(6,*)"pqperm: checkpt 1 count=",count

  !do i=1,mat%n
  !   do j=1,mat%nnzrow(i)
  !!      write(94,*) i, mat%pj(i)%cols(j)
  !      write(75,*) i, mat%pj(i)%cols(j),((mat%pa(i)%submat(ii,jj,j),ii=1,mat%nsubrow),jj=1,mat%nsubcol)
  !   end do
  !end do

  !do i = 1,count
  !   write(6,*)"icor:jcor(",i,")=",icor(i),jcor(i)
  !end do

  !/*-------------------- add entries one by one to diagnl */
  !/* needs recoding so as to scan rows only once instead of 2 */
  do i = 1,count !loop over independent row set 

     ii = icor(i)
     jj = jcor(i)
     !write(6,*)""
     !write(6,*)"pqperm: checkpt 1.1 check element(",ii,jj,") for diagonal placment, Qord(",jj,")=",Qord(jj) 

     if (Qord(jj) .ne. -1) cycle !skip col index if already ordered
     !/*-------------------- */
     row => mat%pj(ii) 
     mrow => mat%pa(ii)
     nzi = mat%nnzrow(ii)
     !/*-------------------- rnz = already assigned cols (either in L or F) */
     itmp=1
     rn = L1Norm2(mrow%submat,itmp,nsubrow)
     rnz = (nzi-1)

     !write(6,*)"pqperm: checkpt 2 row:",i," rnz=",rnz," rn(",ii,"1)=",rn

     !alg 4.3 Matching for Augmented triangular B with dynamic averages
     do k=1,nzi  !loop over columns
        aij = L1Norm2(mrow%submat,k,nsubrow)!should aij be a scalar or matrix
        col = row%cols(k)
        if (Qord(col) >0 ) then !element already assigned
           rn = rn - aij !subtract already assigned col norm from row norm 
           rnz=rnz-1 
           !write(6,*)"pqperm: checkpt 2.5 row:",i," rnz=",rnz," rn(1)=",rn
        else if (Qord(col) == -2) then !if rejected col
           rnz=rnz-1
        endif
     end do

     if (rn < 0.0) cycle !if col 1 of row i is not dominant skip to next row
     !else add to independent set

     !Add i,j to independent set
     numnode=numnode+1 
     Pord(ii) = numnode
     Qord(jj) = numnode

     !/*-------------------- acceptance test among others */    
     do k=1,nzi !loop over columns
        col = row%cols(k)
        if (Qord(col) .ne. -1) cycle !if already added then skip col
        aij = L1Norm2(mrow%submat,k,nsubrow)
        if (rnz*aij > rn) then !reject col ?
           Qord(col) = -2
        else 
           rn = rn - aij
        end if
        rnz=rnz-1
     end do
     
  end do

  !  /*-------------------- number of B nodes */
  nnod = numnode 
  !/*--------------------------------------------------
  !|    end-main-loop - complete permutation arrays
  !|-------------------------------------------------*/


  do i=1,n
     if (Pord(i) < 0) then 
        numnode=numnode+1
        Pord(i) = numnode !fill arbitrary Pord
     end if
  end do

  if (numnode .ne. n) then
     write(6,*)"PQperm: counting error type 1"
     return 
  end if

  numnode = nnod!reset numnode

  do j=1,n
     if (Qord(j) < 0) then 
        numnode=numnode+1
        Qord(j) = numnode !fill arbitrary Qord
     end if
  end do
  
  !/*--------------------              */ 
  if (numnode .ne. n) then
     write(6,*)"PQperm: counting error type 2"
     return
  end if

  !ceb reverse pord
 do j=1,n
    icor(j)=Pord(j)
    jcor(j)=Qord(j)
    !write(6,*)"Pord(",j,")=",Pord(j)
 end do


 do j=1,n
    Pord(icor(j))=j
    Qord(jcor(j))=j
 end do
 !do j=1,n
 !   write(6,*)"Pord(",j,")=",Pord(j)
 !   write(6,*)"Qord(",j,")=",Qord(j)
 !end do

  !/*-------------------- debugging - check permutations */
  !  /* 
  !     printf(" checking P  and Q  :    -->  \n") ;
  !     check_perm(n, Pord) ;
  !     check_perm(n, Qord) ;
  !  */
  !/*--------------------------------------------------
  !|  clean up before returning
  !|-------------------------------------------------*/
  deallocate(icor,STAT=istat)
  deallocate(jcor,STAT=istat)
  return! 0;
  !}
end subroutine PQperm
!/*---------------------------------------------------------------------
!|--------------------------------------------------------------------*/






!int preSel(csptr mat, int *icor, int *jcor, int job, double tol, int *count, int nbnd)  
subroutine preSel(mat, icor, jcor, job, tol, count, nbnd)  
!{
!/*---------------------------------------------------------------------
!| does a preselection of possible diagonal entries. will return a list
!| of "count" bi-indices representing "good" entries to be selected as 
!| diagonal elements -- the order is important (from best to
!| to worst). The list is in the form (icor(ii), jcor(ii)) 
!|
!|      ON ENTRY: 
!       mat   = matrix in csptr format 
!|       tol   = tolerance used for selecting best rows 
!|       job   = indicates whether or not to permute the max entry in 
!|               each row to first position 
!|        NOTE: CAN RECODE BY HAVING JCOR CARRY THE INDEX IN ROW[I] WHERE
!|              MAX IS LOCATED.. 
!|       
!|      ON RETURN: 
!|       icor  = list of row indices of entries selected 
!|       jcor  = list of column indices of entries selected 
!|       count = number of entries selected (size of B block) 
!|--------------------------------------------------------------------*/
  type(cs_type),pointer :: mat
  integer,dimension(:) :: icor,jcor
  integer,intent(inout) ::count
  integer :: job,nbnd
  real(dp) tol

  !internal vars
  integer :: i, k, kmax, col, jmax, countL
  integer,dimension(:),allocatable :: nzi
  real(dp) :: rownorm,t,tmax,wmax
  real(dp),dimension(:),allocatable :: weight
  type(ja_type),pointer ::jcol
  type(submat_type),pointer ::mrow 
  integer :: nsubrow,nsubcol,istat,iisub,jjsub
 
!/*---------------------begin ------------------------------------------------*/
  nsubrow=mat%nsubrow
  nsubcol=mat%nsubcol

  allocate(weight(nbnd),STAT=istat)
  allocate(nzi(nbnd),STAT=istat)
!write(6,*)"presel: checkpt 1"
  !if ( weight==NULL) return 1  

!  /*-------------------- compute max entry for each row */
  wmax = 0.0
  do i=1,nbnd !loop over rows

     jcol => mat%pj(i)
     mrow => mat%pa(i)
     tmax = 0.0 
     kmax = 1 
     rownorm = 0.0
     nzi(i) = 0

     if(mat%nnzrow(i)<=0) write(6,*)"presel: checkpt 2 mat%nnzrow(",i,")=",mat%nnzrow(i)

     do k = 1,mat%nnzrow(i) !loop over cols
        col = jcol%cols(k)
        if (col <= nbnd) then
           nzi(i)=nzi(i)+1
           t = L1Norm2(mrow%submat,k,nsubrow)
           if(abs(t-DBL_EPSILON) > DBL_EPSILON*t) then !if norm is sufficiently nonzero
              rownorm = rownorm + t !add to row sum
              if (tmax < t) then !check for max t column
                 tmax = t
                 kmax = k
              end if
           end if
        end if
     end do

!write(6,*)"presel 3: row ",i," maxcol(",kmax,")=",tmax
     jmax = jcol%cols(kmax) !column index of max t element from row i
     jcor(i) = jmax
     if (job.ne.0 .and. kmax.ne.1) then !shift max col element to first col index
        !matrix swap
        call swap_submat(mrow%submat,1,kmax,nsubrow)
        jcol%cols(kmax) = jcol%cols(1)
        jcol%cols(1) = jmax
     end if
     !/*-------------------- save max diag. dominance ratio  */
     t = tmax / rownorm
     if (wmax < t)  wmax = t 
     weight(i) = t
     !    /* remove!! ALREADY ASSIGNED  */
     !jcor(i) = jmax !ceb removed
  end do

  !/*-------------------- now select rows according to tol */
  countL = 0
  do i=1,nbnd
     t = weight(i)
     col = jcor(i)
     if (t < wmax*tol) cycle ! if row weight less than maxweight*tol add to B set list
     countL=countL+1
     weight(countL) =  t / nzi(i) 
     icor(countL) = i 
     jcor(countL) = col
  end do

  !/*-------------------- sort rows ----------------------------
  !pass in weights array, get new icor,jcor
  call qsortR2I(weight, icor, jcor, 1, countL)
  count = countL;
  deallocate(weight,STAT=istat)
  deallocate(nzi,STAT=istat)
  return
end subroutine preSel
!/*---------------------------------------------------------------------
!|---- end of preSel ---------------------------------------------------
!|--------------------------------------------------------------------*/











subroutine indsetC(mat, bsize, iord, nnod, tol, nbnd) 
  !/*--------------------------------------------------------------------- 
  !| greedy algorithm for independent set ordering -- 
  !|----------------------------------------------------------------------
  !|     Input parameters:
  !|     -----------------
  !|     (mat)  =  matrix in SparRow format
  !|     
  !|     bsize  =  integer (input) the target size of each block.
  !|               each block is of size >= bsize. 
  !|
  !|     w      =  weight factors for the selection of the elements in the
  !|               independent set. If w(i) is small i will be left for the
  !|               vertex cover set. 
  !|
  !|     tol    =  a tolerance for excluding a row from independent set.
  !|
  !|     nbnd   =  number of interior variables.
  !|
  !|     Output parameters:
  !|     ------------------ 
  !|     iord   = permutation array corresponding to the independent set 
  !|     ordering.  Row number i will become row number iord[i] in 
  !|     permuted matrix.
  !|     
  !|     nnod   = (output) number of elements in the independent set.
  !|     
  !|----------------------------------------------------------------------- 
  !|     the algorithm searches nodes in lexicographic order and groups
  !|     the (BSIZE-1) nearest nodes of the current to form a block of
  !|     size BSIZE. The current algorithm does not use values of the matrix.
  !|---------------------------------------------------------------------*/ 

  type(cs_type),pointer ::mat
  integer,intent(in) :: bsize, nbnd
  integer,intent(inout) :: nnod
  integer,intent(inout),dimension(:):: iord
  real(dp),intent(in)::tol

  !/*   local variables   */
  integer:: nod, jcount, lastlev, begin, last0, last, nback, mid,istat
  integer::  j1, j2, jcol, inod, jnod, i,j,k, iisub,jjsub,jcount0, begin0,prog,n,pos
  integer,dimension(:),allocatable :: riord
  type(ja_type),pointer::rowj
  real(dp),dimension(:),allocatable::w
  type(cs_type),pointer ::matT
  type(cs_type),pointer::gmat

  nullify(matT)
  nullify(gmat)
  nullify(rowj)

  n=mat%n
!write(6,*),"indsetC checkpt 0"
  !/*-----------------------------------------------------------------------*/
  allocate(riord(n),STAT=istat)
  allocate(w(n),STAT=istat)

  !/!*  	 weights to compute the weights for  input matrix.. */

  !do i=1,mat%n
  !   do j=1,mat%nnzrow(i)
  !      write(110,*) i, mat%pj(i)%cols(j)!,((mat%pa(i)%submat(iisub,jjsub,j),iisub=1,mat%nsubrow),jjsub=1,mat%nsubcol)
  !   end do
  !end do

  call setupCS(matT,mat%n,mat%nsubrow,mat%nsubcol,1)!creates new matrix matT 
  call SparTran(mat, matT, 1, 0)

  !do i=1,matT%n
  !   do j=1,matT%nnzrow(i)
  !      write(111,*) i, matT%pj(i)%cols(j)!,((matT%pa(i)%submat(iisub,jjsub,j),iisub=1,matT%nsubrow),jjsub=1,matT%nsubcol)
  !   end do
  !end do

  !call SparTran(matT, mat, 1, 1) 
  call weightsC(mat, w) 



  !/*---------------------------------------------------------------------- 
  !| scan all nodes first to eliminate those not satisfying DD criterion 
  !+----------------------------------------------------------------------*/
  nback = n
  nod = 1
  !initialize
  do j=1,n
     iord(j) = -1
  end do
  
  do j = n,nbnd+1,-1
     call add2com(nback, j, iord, riord)!decrements nback
  end do
  
  do j=1, nbnd
     if (w(j) < tol) then
        !write(6,*),"w[",j,"]=",w(j)
        call add2com(nback, j, iord, riord)
        nod=nod+1
     end if
  end do
  
  last = 0
  nod=0
  do
     nod=nod+1
     if(nod>n) exit

     do
        if(iord(nod) .eq. -1) exit
        nod=nod+1
        if(nod>n) goto 50
     end do
     
     !/*-------------------- initialize level-set - contains nod (only)*/
     call add2is(last, nod, iord, riord)

     begin   = last
     begin0  = begin 
     lastlev = begin
     jcount  = 1
     !/*----------------------------------------------------------------------
     !|     put all the nearest neighbor nodes of the current node into
     !|     the block until the number is BSIZE.
     !|---------------------------------------------------------------------*/
     prog = 1
     !while (jcount < bsize && prog.ne.0) {
     do
        !if(jcount .gt. bsize .or. prog.eq.0) then
        if(jcount .gt. bsize) then
          write(6,*)"exiting indsetc: jcount==bsize"
           exit
        else if(prog.eq.0) then
           write(6,*)"exiting indsetc: prog==0"
           exit
        end if
        !/*--------------------   traverse all the current level-set   */
        last0 = last
        jcount0 = jcount
        
        do inod=begin,last0
           jnod = riord(inod) 
           !/*--------------------   This assumes A is not symmetric.   */
           gmat => mat 
           do k=1,2 
              rowj => gmat%pj(jnod)
              do j=1,gmat%nnzrow(jnod) 
                 jcol = rowj%cols(j)
                 if (iord(jcol) .eq. -1 ) then	
                    call add2is(last, jcol, iord, riord)!increments last
                    jcount=jcount+1
                 end if
              end do
              gmat => matT	
           end do
        end do

        prog=0
        if(jcount>jcount0)prog=1
        lastlev = begin
        begin = last0+1
     end do!end while
     
     !/*-----------------------------------------------------------------------
     !| the neighbors of elements of last level go to the complement   
     !| gmat loop over original matrix and its transpose 
     !+-----------------------------------------------------------------------*/ 
     gmat => mat 
     do k=1,2 
        do inod=lastlev,last  	
           jnod = riord(inod) 
           rowj => gmat%pj(jnod)
           do j=1,gmat%nnzrow(jnod)
               jcol = rowj%cols(j)
              if (iord(jcol) .eq. -1) call add2com(nback, jcol, iord, riord)
           end do
        end do
        gmat => matT 	
     end do

     !/*   reverse ordering for this level   */
     mid = (begin0+last) / 2
     do inod=begin0,mid 
        j = last - inod + begin0
        jnod = riord(inod)
        riord(inod) = riord(j)
        riord(j) = jnod
     end do
  end do!loop over nod

  !/*--------------------------------------------------
  !|  end-main-loop
  !|-------------------------------------------------*/


  !/*-------------------- relabel nodes of vertex cover   */
50 continue
  nnod = last
  j1 = nnod
  do j2=nnod+1,n
     if (iord(riord(j2)) > -1) then
        j1=j1+1
        if (j1 .ne. j2) then
           j = riord(j2)
           riord(j2) = riord(j1)
           riord(j1) = j
        end if
     end if
  end do

  !/*-------------------- obtain reverse permutation array   */
  do j=1,n
     iord(riord(j)) = j
  !   !write(46,*) "iord(",riord(j),")=",iord(riord(j))
  end do

!write(6,*),"indsetC checkpt 5: nnod=",nnod
!do i=1,mat%n
!   do j=1,mat%nnzrow(i)
!      write(46,*) iord(i),iord(mat%pj(i)%cols(j))
!   end do
!end do
!write(6,*),"indsetC checkpt 5"
!do i=1,mat%n
!   do j=1,mat%nnzrow(i)
!      write(47,*) riord(iord(i)),riord(iord(mat%pj(i)%cols(j)))
!   end do
!end do
!write(6,*),"indsetC checkpt 6"

  call cleanCS(matT) !implement this
  deallocate(riord,STAT=istat)
  deallocate(w,STAT=istat)
  return
end subroutine indsetC
!/*---------------------------------------------------------------------
!|-----end-of-indsetC---------------------------------------------------
!|--------------------------------------------------------------------*/





!/*----------------------------------------------------------------------
!|   adds element nod to independent set
!|---------------------------------------------------------------------*/
subroutine add2is(last, nod, iord, riord)
  integer,intent(in)::nod
  integer,intent(inout)::last 
  integer,intent(inout),dimension(:)::iord, riord
  last=last+1
!write(6,*) "adding row ",nod," to independent set ind ",last
  iord(nod) = last
  riord(last) = nod
end subroutine add2is

!/*----------------------------------------------------------------------
!|   adds element nod to schur complement set
!|---------------------------------------------------------------------*/
subroutine add2com(nback, nod, iord, riord) 
  integer,intent(in)::nod
  integer,intent(inout):: nback
  integer,intent(inout), dimension(:):: iord, riord
!write(6,*) "adding row ",nod," to complement set ind ",nback
  iord(nod) = nback
  riord(nback) = nod
  nback=nback-1
end subroutine add2com
!=======================================================================================================


subroutine copy_to_CS(A,ia,iau,ja,amat)
  !assumes amat exists
  integer,intent(in),dimension(:)::ia,iau,ja
  real(dp),intent(in),dimension(:,:,:)::A
  type(cs_type),pointer::amat
  integer:: istat,i,j,jcol,iisub,jjsub
  allocate(amat%pd(amat%n),STAT=istat)
  do i=1,amat%n
     amat%nnzrow(i)=ia(i+1)-ia(i)
     amat%pd(i)=iau(i)-ia(i)+1 
     !amat%pd(i)=iau(i)-ia(i)+1 !ceb make this orig nnzrow 
!write(6,*) "copy_to_cs: amat%nnzrow(",i,")=",amat%nnzrow(i)," amat%pd(",i,")=",amat%pd(i)
     allocate(amat%pj(i)%cols(amat%nnzrow(i)),STAT=istat)     
     allocate(amat%pa(i)%submat(amat%nsubrow,amat%nsubcol,amat%nnzrow(i)),STAT=istat)     
     jcol=0
     do j=ia(i),ia(i+1)-1
        jcol=jcol+1
        amat%pj(i)%cols(jcol)=ja(j)
        do iisub=1,amat%nsubrow
           do jjsub=1,amat%nsubcol
              amat%pa(i)%submat(iisub,jjsub,jcol)=A(iisub,jjsub,j)
           end do
        end do
     end do
  end do
  return
end subroutine copy_to_CS
!|--------------------------------------------------------------------*/







subroutine setupCS(mat,len,nsubr,nsubc,job)
  !/*----------------------------------------------------------------------
  !| Initialize SparRow structs.
  !|----------------------------------------------------------------------
  !| on entry:
  !|==========
  !| ( mat )  =  Pointer to a SparRow struct.
  !|     len   =  size of matrix
  !|     job   =  0: pattern only
  !|              1: data and pattern
  !|
  !| On return:
  !===========
  !|
  !|  mat->n
  !|      ->*nnzrow
  !|      ->**ja
  !|      ->**ma
  !|
  !| integer value returned:
  !|             0   --> successful return.
  !|             1   --> memory allocation error.
  !|--------------------------------------------------------------------*/
  type(cs_type),pointer::mat
  integer::len,job,nsubr,nsubc
  integer:: istat

!write(6,*)"setupCS checkpt 0"
  !if(associated(mat) .and. allocated(mat))
  call cleanCS(mat) 
!write(6,*)"setupCS checkpt 1"
  allocate(mat,STAT=istat)
!write(6,*)"setupCS checkpt 2"
  mat%n = len
  mat%nsubrow=nsubr
  mat%nsubcol=nsubc
!write(6,*)"setupCS checkpt 0"
  allocate(mat%nnzrow(len),STAT=istat)
  allocate(mat%pd(len),STAT=istat)
!write(6,*)"setupCS checkpt 4"
  allocate(mat%pj(len),STAT=istat)
!write(6,*)"setupCS checkpt 5"
  if( job .eq. 1 ) then 
     allocate(mat%pa(len),STAT=istat)
!write(6,*)"setupCS checkpt 5.1"
  else
     !nullify(mat->pa)! = NULL
  end if
  return
end subroutine setupCS
!/*---------------------------------------------------------------------
!|     end of setupCS
!|--------------------------------------------------------------------*/





!=======================================================================================================
subroutine cleanCS(mat)
  !/*----------------------------------------------------------------------
  !| Free up memory allocated for SparRow structs.
  !|----------------------------------------------------------------------
  !| on entry:
  !|==========
  !| ( mat )  =  Pointer to a SparRow struct.
  !|--------------------------------------------------------------------*/
  !type(csptr),intent(inout)::mat
  type(cs_type),pointer::mat
  integer:: i,istat
!write(6,*)"cleanCS checkpt 0"
  if (.not.associated(mat)) return! 0;
  !if ((.not.associated(mat)) .or. (.not.allocated(mat))) return! 0;
!write(6,*)"cleanCS checkpt 1"
  if (mat%n < 1) return 
!write(6,*)"cleanCS checkpt 2 mat%n=",mat%n
  
  !if(allocated(mat%nnzrow))deallocate(mat%nnzrow,STAT=istat)
!write(6,*)"cleanCS checkpt 2.1"
     
!write(6,*)"cleanCS checkpt 3.1 size(mat%pj)=",size(mat%pj)
!write(6,*)"cleanCS checkpt 3.2 size(mat%pa)=",size(mat%pa)

  if(size(mat%pa)>0 .and. size(mat%pj)>0) then

     do i=1,mat%n
!        write(6,*)"cleanCS checkpt 3"
        
        if (allocated(mat%pa(i)%submat) .and. size(mat%pa(i)%submat)>0) then
!           write(6,*)"cleanCS checkpt 3.3 size(mat%pa(",i,")%submat=",size(mat%pa(i)%submat),")"
           deallocate(mat%pa(i)%submat,STAT=istat)
        end if
     
!        write(6,*)"cleanCS checkpt 4"
     
        if (allocated(mat%pj(i)%cols) .and. size(mat%pj(i)%cols)>0) then
!           write(6,*)"cleanCS checkpt 4.3 size(mat%pj(",i,")%cols)=",size(mat%pj(i)%cols),")"
           deallocate(mat%pj(i)%cols,STAT=istat)
        end if
!        write(6,*)"cleanCS checkpt 5"
     end do
  
!     write(6,*)"cleanCS checkpt 6"
     deallocate(mat%pa,STAT=istat)
!     write(6,*)"cleanCS checkpt 7"
     deallocate(mat%pj,STAT=istat)
!     write(6,*)"cleanCS checkpt 8"
     if(size(mat%nnzrow)>0) then
        deallocate(mat%nnzrow,STAT=istat)
        deallocate(mat%pd,STAT=istat)
     end if

  end if
  
!  write(6,*)"cleanCS checkpt 9"
  nullify(mat)!deallocate(amat)
!  write(6,*)"cleanCS checkpt 10"
  return 
end subroutine cleanCS
!/*---------------------------------------------------------------------
!|     end of cleanCS
!|--------------------------------------------------------------------*/
!=======================================================================================================









!int csSplit4(csptr amat, int bsize, int csize, csptr B, csptr F, csptr E, csptr C)
subroutine csSplit4(amat, bsize, csize, B, F, E, C)
  !{
  !/*---------------------------------------------------------------------
  !| Convert permuted csrmat struct to PerMat4 struct 
  !|                - matrix already permuted
  !|----------------------------------------------------------------------
  !| on entry:
  !|========== 
  !| ( amat )  =  Matrix stored in SparRow format.
  !|              Internal pointers (and associated memory) destroyed before
  !|              return.
  !|
  !| On return:
  !|===========
  !|
  !| B, E, F, C = 4 blocks in 
  !| 
  !|          | B   F |      
  !|   Amat = |       | 
  !|          | E   C | 
  !| 
  !|
  !|       integer value returned:
  !|             0   --> successful return.
  !|             1   --> memory allocation error.
  !|--------------------------------------------------------------------*/
  type(cs_type),pointer:: amat,B,F,E,C
  integer :: bsize,csize
  integer :: i,j,k, j1, numr, numl, ind, newj, rowz,nsubrow,nsubcol
  integer :: iisub,jjsub,kksub,istat,itmp
  type(ja_type),pointer::rowj
  type(ja_type) :: new1j,new2j
  type(submat_type),pointer::rowm
  type(submat_type)::new1m,new2m

  nullify(B)
  nullify(F)
  nullify(E)
  nullify(C)

  nsubrow=amat%nsubrow
  nsubcol=amat%nsubcol

  !/*---------------------------------------------------------------------
  !|     Sort the matrix and separate into   |  B  F  |
  !|                                         |        |
  !|                                         |  E  C  |
  !|--------------------------------------------------------------------*/
!write(6,*)"cssplit4: checkpt 0"
  call setupCS(B,bsize,amat%nsubrow,amat%nsubcol,1) 
  call setupCS(F,bsize,amat%nsubrow,amat%nsubcol,1)
  call setupCS(E,csize,amat%nsubrow,amat%nsubcol,1)
  call setupCS(C,csize,amat%nsubrow,amat%nsubcol,1)
!write(6,*)"cssplit4: checkpt 1"
  allocate(new1j%cols(bsize),STAT=istat)
  allocate(new2j%cols(csize),STAT=istat)
  allocate(new1m%submat(nsubrow,nsubcol,bsize),STAT=istat)! = (int *) Malloc(bsize*sizeof(int), "csSplit4:1" );
  allocate(new2m%submat(nsubrow,nsubcol,csize),STAT=istat)
!write(6,*)"cssplit4: checkpt 2"

  !ceb
  allocate(B%pd(bsize),STAT=istat)
  allocate(C%pd(csize),STAT=istat)
  !ceb  

  !/*    B and F blocks */ 
  do j=1,bsize
     numl = 0
     numr = 0
     rowz = amat%nnzrow(j)
     rowj => amat%pj(j)
     rowm => amat%pa(j)
     do j1=1,rowz
        !if (rowj%cols(j1)<bsize) then
        if (rowj%cols(j1).le.bsize) then
           numl = numl+1
        else 
           numr = numr+1
        end if
     end do

     B%nnzrow(j) = numl
     F%nnzrow(j) = numr
!write(6,*)"cssplit4: checkpt 3"

     if (numl.gt.0) then
        allocate(B%pj(j)%cols(numl),STAT=istat)! = (int *) Malloc(numl*sizeof(int), "csSplit4:5" );
        allocate(B%pa(j)%submat(nsubrow,nsubcol,numl),STAT=istat)! = (double *) Malloc(numl*sizeof(double), "csSplit4:6" );
     end if

     if (numr.gt.0) then
        allocate(F%pj(j)%cols(numr),STAT=istat) ! = (int *) Malloc(numr*sizeof(int), "csSplit4:7" );
        allocate(F%pa(j)%submat(nsubrow,nsubcol,numr),STAT=istat)! = (double *) Malloc(numr*sizeof(double), "csSplit4:8" );
     end if


     numl = 0
     numr = 0

     do j1=1,rowz
        newj = rowj%cols(j1)
        !if (newj<bsize) then
        if (newj.le.bsize) then
           numl = numl+1
           new1j%cols(numl) = newj
           !copy submat
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 new1m%submat(iisub,jjsub,numl) = rowm%submat(iisub,jjsub,j1)
              end do
           end do
        else 
           numr = numr+1
           new2j%cols(numr) = newj - bsize
           !copy submat
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 new2m%submat(iisub,jjsub,numr) = rowm%submat(iisub,jjsub,j1)
              end do
           end do
        end if
     end do
     
     B%pd(j)=numl!ceb make this orig nnzrow
     do k=1,numl
        B%pj(j)%cols(k)=new1j%cols(k)
     !   if(j==B%pj(j)%cols(k)) B%pd(j)=k!ceb
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              B%pa(j)%submat(iisub,jjsub,k)= new1m%submat(iisub,jjsub,k)
           end do
        end do
     end do
     do k=1,numr
        F%pj(j)%cols(k)=new2j%cols(k)
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              F%pa(j)%submat(iisub,jjsub,k)= new2m%submat(iisub,jjsub,k)
           end do
        end do
     end do

  end do
!write(6,*)"cssplit4: checkpt 5"
  
  !/!*    E and C blocks */
  do j=1,csize !{
     numl = 0
     numr = 0
     ind = bsize + j
     rowz = amat%nnzrow(ind)
     rowj => amat%pj(ind)
     rowm => amat%pa(ind)
     do j1=1,rowz
        !if (rowj%cols(j1)<bsize) then
        if (rowj%cols(j1).le.bsize) then
           numl=numl+1
        else 
           numr=numr+1
        end if
     end do
     E%nnzrow(j) = numl
     C%nnzrow(j) = numr
     if (numl>0) then
	allocate(E%pj(j)%cols(numl),STAT=istat)
        allocate(E%pa(j)%submat(nsubrow,nsubcol,numl),STAT=istat)
     end if

     if (numr>0) then
	allocate(C%pj(j)%cols(numr),STAT=istat)
        allocate(C%pa(j)%submat(nsubrow,nsubcol,numr),STAT=istat)
     end if
     numl = 0
     numr = 0
     do j1=1,rowz 
	newj = rowj%cols(j1)
        !if (newj<bsize) then
        if (newj.le.bsize) then
           numl=numl+1
           new1j%cols(numl) = newj
           !copy submat
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 new1m%submat(iisub,jjsub,numl) = rowm%submat(iisub,jjsub,j1)
              end do
           end do
        else 
           numr=numr+1
           new2j%cols(numr) = newj - bsize;
           !copy submat
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 new2m%submat(iisub,jjsub,numr) = rowm%submat(iisub,jjsub,j1)
              end do
           end do
	end if
     end do


     !memcpy E
     do k=1,numl
        E%pj(j)%cols(k)=new1j%cols(k)
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              E%pa(j)%submat(iisub,jjsub,k)= new1m%submat(iisub,jjsub,k)
           end do
        end do
     end do
     !memcpy C

     C%pd(j)=numr !ceb make this orig nnzrow
     do k=1,numr
        C%pj(j)%cols(k)=new2j%cols(k)
     !   if(j==C%pj(j)%cols(k)) C%pd(j)=k!ceb
        do iisub=1,nsubrow
           do jjsub=1,nsubcol
              C%pa(j)%submat(iisub,jjsub,k)= new2m%submat(iisub,jjsub,k)
           end do
        end do
     end do
  end do
!write(6,*)"cssplit4: checkpt 6"
  nullify(rowj)
  nullify(rowm)
  deallocate(new1j%cols,STAT=istat)
  deallocate(new2j%cols,STAT=istat)
  deallocate(new1m%submat,STAT=istat)
  deallocate(new2m%submat,STAT=istat)
!write(6,*)"cssplit4: checkpt 7"

  return !0;
111 continue
  return !1;
end subroutine csSplit4
!/*---------------------------------------------------------------------
!|     end of csSplit4
!|--------------------------------------------------------------------*/




subroutine createP4 (amat) 
  !/*----------------------------------------------------------------------
  !| initialize PerMat4 struct given the F, E, blocks.  
  !|----------------------------------------------------------------------
  !type(p4_type),intent(inout)::amat
  type(p4_type),pointer::amat
  allocate(amat)!?

  amat%nsubrow=0!F%nsubrow
  amat%nsubcol=0!F%nsubcol
  amat%n = 0!n
  amat%nB =0! Bn 
  amat%symperm=0
  amat%mylev=0
  nullify(amat%prev)
  nullify(amat%next)
  nullify(amat%L)
  nullify(amat%U)
  nullify(amat%E)
  nullify(amat%F)
  nullify(amat%rperm)
  nullify(amat%perm)
  nullify(amat%D1)
  nullify(amat%D2)
  nullify(amat%wk)
  return
end subroutine createP4
!/*---------------------------------------------------------------------
!|     end of setupP4 
!|--------------------------------------------------------------------*/

subroutine setupP4 (amat, Bn, Cn, F, E, lev) 
  !/*----------------------------------------------------------------------
  !| initialize PerMat4 struct given the F, E, blocks.  
  !|----------------------------------------------------------------------
  !| on entry:
  !|==========
  !| ( amat )  =  Pointer to a PerMat4 struct.
  !|     Bn    =  size of B block
  !|     Cn    =  size of C block
  !|     F, E  = the two blocks to be assigned to srtruct - without the
  !|
  !| On return:
  !|===========
  !|
  !|  amat->L                for each block: amat->M->n
  !|      ->U                                       ->nnzrow
  !|      ->E                                       ->pj
  !|      ->F                                       ->pa
  !|      ->perm
  !|      ->rperm       (if meth[1] > 0)
  !|      ->D1          (if meth[2] > 0)
  !|      ->D2          (if meth[3] > 0)
  !|
  !|  Scaling arrays are initialized to 1.0.
  !|
  !|       integer value returned:
  !|             0   --> successful return.
  !|             1   --> memory allocation error.
  !|--------------------------------------------------------------------*/
  !type(p4ptr),intent(inout)::amat
  type(p4_type),pointer::amat
  type(cs_type),pointer::F,E
  integer,intent(in) :: Bn,Cn
  integer,intent(in), optional :: lev
  !internal vars
  integer :: n,istat,nsubrow,nsubcol

  if(present(lev)) then
     amat%mylev=lev
  else
     amat%mylev=0
  end if

  n = Bn+Cn  !/* size n */
  amat%nsubrow=F%nsubrow
  amat%nsubcol=F%nsubcol
  amat%n = n
  amat%nB = Bn 
  !/*   assign space for wk -- note that this is only done at 1st level
  !   at other levels, copy pointer of wk from previous level */
!write(6,*)"setupP4: checkpt 1"
  !ceb need to determine how pointer is points to null
  if (.not.associated(amat%prev)) then!/* wk has 2 * n entries now */ 
     allocate(amat%wk(amat%nsubrow,2*n),STAT=istat)
  else 
     amat%wk => amat%prev%wk !ceb so far amat%wk is not a pointer
  end if
!write(6,*)"setupP4: checkpt 2 setting up L: Bn=",Bn
  !/*-------------------- L and U */ 
  call setupCS(amat%L,Bn,amat%nsubrow,amat%nsubcol,1)
!write(6,*)"setupP4: checkpt 3 setting up U: Bn=",Bn
  !/*    fprintf(stdout,"  -- BN %d   Cn   %d \n", Bn,Cn);  */
  call setupCS(amat%U,Bn,amat%nsubrow,amat%nsubcol,1)
!write(6,*)"setupP4: checkpt 4"
  allocate(amat%F,STAT=istat)
  allocate(amat%E,STAT=istat)
!write(6,*)"setupP4: checkpt 5"
  amat%F => F !F passed in is pointer
  amat%E => E !E passed in is pointer
!write(6,*)"setupP4: checkpt 6"
  nullify(amat%rperm)
  nullify(amat%perm)
  nullify(amat%D1)
  nullify(amat%D2)
  !nullify(amat%prev)
  !nullify(amat%next)
  return
end subroutine setupP4
!/*---------------------------------------------------------------------
!|     end of setupP4 
!|--------------------------------------------------------------------*/




subroutine cleanP4(amat)
  !/*----------------------------------------------------------------------
  !| Free up memory allocated for Per4Mat structs.
  !|----------------------------------------------------------------------
  !| on entry:
  !|==========
  !| ( amat )  =  Pointer to a Per4Mat struct.
  !|--------------------------------------------------------------------*/
  type(p4_type),pointer::amat
  integer:: istat
  if (.not.associated(amat)) return! 0;
  
  if (amat%n < 1) return
  
  if (associated(amat%perm)) then
     deallocate(amat%perm,STAT=istat) 
     nullify(amat%perm)
  end if
  
  if (amat%symperm.eq.0) then 
     if (associated(amat%rperm)) then
        deallocate(amat%rperm,STAT=istat) 
        nullify(amat%rperm)
     end if
  end if
  
  if (associated(amat%F)) then
     call cleanCS(amat%F) 
     nullify(amat%F)
  end if
  if (associated(amat%E)) then
     call cleanCS(amat%E) 
     nullify(amat%E)
  end if
  if (associated(amat%L)) then
     call cleanCS(amat%L)
     nullify(amat%L)
  end if
  if (associated(amat%U)) then
     call cleanCS(amat%U)
     nullify(amat%U)
  end if
  
  if (.not.associated(amat%prev)) then
     if (associated(amat%wk)) then
        deallocate(amat%wk,STAT=istat)  
        nullify(amat%wk)
     end if
  end if
  if (associated(amat%D1)) then
     deallocate(amat%D1,STAT=istat)
     nullify(amat%D1)
  end if
  if (associated(amat%D2)) then 
     deallocate(amat%D2,STAT=istat)
     nullify(amat%D2)
  end if
  nullify(amat)
  return
end subroutine cleanP4
!/*---------------------------------------------------------------------
!|     end of cleanP4
!|--------------------------------------------------------------------*/

subroutine createILUT (amat) 
  !/*----------------------------------------------------------------------
  !type(ilut_type),intent(inout)::amat
  type(ilut_type),pointer::amat
  allocate(amat)
  amat%nsubrow=0!F%nsubrow
  amat%nsubcol=0!F%nsubcol
  nullify(amat%C)
  nullify(amat%L)
  nullify(amat%U)
  nullify(amat%rperm)
  nullify(amat%perm)
  nullify(amat%perm2)
  nullify(amat%D1)
  nullify(amat%D2)
  nullify(amat%wk)
  return
end subroutine createILUT

!int setupILUT(ilutptr amat, int len)
subroutine setupILUT(amat,len,nsubrow,nsubcol)
  !/*----------------------------------------------------------------------
  !| Allocate pointers for ILUTfac structs.
  !|----------------------------------------------------------------------
  !| on entry:
  !|==========
  !| ( amat )  =  Pointer to a ILUTfac struct.
  !|     len   =  size of L U  blocks
  !|
  !| On return:
  !|===========
  !|
  !|  amat->L                for each block: amat->M->n
  !|      ->U                                       ->nnzrow
  !|                                                ->pj
  !|                                                ->pa
  !|      ->rperm       (if meth[0] > 0)
  !|      ->perm2       (if meth[1] > 0)
  !|      ->D1          (if meth[2] > 0)
  !|      ->D2          (if meth[3] > 0)
  !|
  !|  Permutation arrays are initialized to the identity.
  !|  Scaling arrays are initialized to 1.0.
  !|
  !|       integer value returned:
  !|             0   --> successful return.
  !|             1   --> memory allocation error.
  !|--------------------------------------------------------------------*/
  type(ilut_type),pointer :: amat
  integer :: len,nsubrow,nsubcol,istat
  if(.not.associated(amat)) allocate(amat,STAT=istat)
  amat%n = len
  amat%nsubrow=nsubrow
  amat%nsubcol=nsubcol
  allocate(amat%wk(nsubrow,2*len),STAT=istat)
  allocate(amat%L,STAT=istat)
  amat%L%n=0
!write(6,*)"setupILUT checkpt 1"
  call setupCS(amat%L,len,amat%nsubrow,amat%nsubcol,1)!) return 1;
!write(6,*)"setupILUT checkpt 2"
  allocate(amat%U,STAT=istat)
  amat%U%n=0
!write(6,*)"setupILUT checkpt 3"
  call setupCS(amat%U,len,amat%nsubrow,amat%nsubcol,1)!) return 1;
!write(6,*)"setupILUT checkpt 4"
  return !0;
  !/*---------------------------------------------------------------------
end subroutine setupILUT
  !    |---------------------------------------------------------------------*/



!int cleanILUT(ilutptr amat, int indic)
subroutine cleanILUT(amat, indic)
  !/*----------------------------------------------------------------------
  !| Free up memory allocated for IluSpar structs.
  !|----------------------------------------------------------------------
  !| on entry:
  !|==========
  !| ( amat )  =  Pointer to a IluSpar struct.
  !|  indic    = indicator for number of levels.  indic=0 -> zero level.
  !|--------------------------------------------------------------------*/
  type(ilut_type),pointer :: amat
  integer :: indic,istat
  
  if(.not.associated(amat)) return

  if (associated(amat%wk)) then
     deallocate(amat%wk,STAT=istat)
     nullify(amat%wk)
  end if

  if(associated(amat%L)) call cleanCS(amat%L)
  if(associated(amat%U)) call cleanCS(amat%U)
  if (indic.eq.1) call cleanCS(amat%C)  
  !/*-------------------- nonsymmetric permutation */
  if (associated(amat%rperm)) then
     deallocate(amat%rperm,STAT=istat)
     nullify(amat%rperm)
  end if

  if (associated(amat%perm)) then
     deallocate(amat%perm,STAT=istat) 
     nullify(amat%perm)
  end if
  !/*-------------------- ilutp permutation */
  if (associated(amat%perm2)) then
     deallocate(amat%perm2,STAT=istat)
     nullify(amat%perm2)
  end if
  !/*-------------------- diagonal scalings */
  if (associated(amat%D1)) then
     deallocate(amat%D1,STAT=istat)
     nullify(amat%D1)
  end if
  if (associated(amat%D2)) then
     deallocate(amat%D2,STAT=istat)
     nullify(amat%D2)
  end if
  nullify(amat)
  return !0
  !/*---------------------------------------------------------------------
end subroutine cleanILUT
  !    |---------------------------------------------------------------------*/





subroutine roscalC(mata, diag, nrm)
  !/*---------------------------------------------------------------------
  !|
  !| This routine scales each row of mata so that the norm is 1.
  !|
  !|----------------------------------------------------------------------
  !| on entry:
  !| mata  = the matrix (in SparRow form)
  !| nrm   = type of norm
  !|          0 (\infty),  1 or 2
  !|
  !| on return
  !| diag  = diag[j] = 1/norm(row[j])
  !|
  !|     0 --> normal return
  !|     j --> row j is a zero row
  !|--------------------------------------------------------------------*/
  type(cs_type),target :: mata
  real(dp),dimension(:,:)::diag
  integer :: nrm
  !/*   local variables    */
  integer :: i, k,iisub,kksub,nsubrow,nsubcol
  !double *kr, scal;
  type(submat_type),pointer::kr
  real(dp)::scal
  nsubrow=mata%nsubrow
  nsubcol=mata%nsubcol
write(6,*)"calling roscalC"  
  do i=1,mata%n
     kr => mata%pa(i)

     do iisub =1,nsubrow!loop subrows

        scal = 0.0
        if (nrm .eq. 0) then ! infinity norm
           !find largest diag val in subcols k,kk
           do k=1,mata%nnzrow(i)!loop columns
              do kksub=1,nsubcol!loop subcols
                 if (abs(kr%submat(iisub,kksub,k)) > scal) scal = abs(kr%submat(iisub,iisub,k))
              end do
           end do
        else if (nrm .eq. 1) then !L1 norm
           do k=1,mata%nnzrow(i)!loop columns
              do kksub=1,nsubcol!loop subcols
                 scal = scal + abs(kr%submat(iisub,kksub,k))
              end do
           end do
        else   !/* nrm = 2 */ !L2 norm
           do k=1,mata%nnzrow(i)!loop columns
              do kksub=1,nsubcol!loop subcols
                 scal = scal + kr%submat(iisub,kksub,k)*kr%submat(iisub,kksub,k)
              end do
           end do
        end if
        if (nrm .eq. 2) scal = sqrt(scal);
        if ( abs(scal-DBL_EPSILON) .le. DBL_EPSILON*abs(scal) ) then
           scal = 1.0
           !/* YS. return i+1; */
        else 
           scal = 1.0 / scal 
        end if

        diag(iisub,i) = scal
!write(6,*)"roscal: diag(",iisub,i,")=",diag(iisub,i)
        do k=1,mata%nnzrow(i)!loop columns
           do kksub=1,nsubcol!loop subcols
              kr%submat(iisub,kksub,k) = kr%submat(iisub,kksub,k) * scal
           end do
        end do

     end do

  end do
  return !0;
end subroutine roscalC
!/*---------------end of roscalC-----------------------------------------
!----------------------------------------------------------------------*/


subroutine coscalC(mata, diag, nrm)
  !/*---------------------------------------------------------------------
  !|
  !| This routine scales each column of mata so that the norm is 1.
  !|
  !|----------------------------------------------------------------------
  !| on entry:
  !| mata  = the matrix (in SparRow form)
  !| nrm   = type of norm
  !|          0 (\infty),  1 or 2
  !|
  !| on return
  !| diag  = diag[j] = 1/norm(row[j])
  !|
  !|     0 --> normal return
  !|     j --> column j is a zero column
  !|--------------------------------------------------------------------*/
  type(cs_type),target :: mata
  real(dp),pointer,dimension(:,:) :: diag !should this have subrow dimensions as well
  integer :: nrm
  !/*   local variables    */
  integer :: i, j, k,iisub,jjsub,kksub,nsubrow,nsubcol
  type(submat_type),pointer :: kr
  type(ja_type),pointer :: ki
  nsubrow=mata%nsubrow
  nsubcol=mata%nsubcol
write(6,*)"calling coscalC"  

  do i=1,mata%n
     do iisub=1,nsubrow
        diag(iisub,i) = 0.0
     end do
  end do
  !/*---------------------------------------
  !|   compute the norm of each column
  !|--------------------------------------*/
  do i=1, mata%n
     kr => mata%pa(i)
     ki => mata%pj(i)
     do iisub=1,nsubrow
        if (nrm .eq. 0) then
           do k=1, mata%nnzrow(i)
              j = ki%cols(k)
              do kksub=1,nsubcol       
                 if (abs(kr%submat(iisub,kksub,k)) > diag(iisub,j)) diag(iisub,j) = abs(kr%submat(iisub,kksub,k))
              end do
           end do
        else if (nrm .eq. 1) then
           do k=1,mata%nnzrow(i)
              do kksub=1,nsubcol
                 diag(iisub,ki%cols(k)) = diag(iisub,ki%cols(k))+ abs(kr%submat(iisub,kksub,k))
              end do
           end do
        else   !/*  nrm = 2 */
           do k=1,mata%nnzrow(i)!loop cols
              do kksub=1,nsubcol
                 diag(iisub,ki%cols(k)) = diag(iisub,ki%cols(k)) + kr%submat(iisub,kksub,k)*kr%submat(iisub,kksub,k) 
              end do
           end do
        end if
     end do
  end do
  if (nrm .eq. 2) then
     do i=1,mata%n!loop rows
        do iisub=1,nsubrow!loop subrows
           diag(iisub,i) = sqrt(diag(iisub,i))
        end do
     end do
  end if
  !/*---------------------------------------
  !|   invert
  !|--------------------------------------*/
  do i=1,mata%n!loop rows
     do iisub=1,nsubrow!loop subrows
        if (abs(diag(iisub,i)-DBL_EPSILON) .le. DBL_EPSILON*diag(iisub,i)) then
           !/* return i+1;*/
           diag(iisub,i) = 1.0 
        else 
           diag(iisub,i) = 1.0 / diag(iisub,i)
        end if
!write(6,*)"coscal: diag(",iisub,i,")=",diag(iisub,i),1.0/diag(iisub,i)
     end do
  end do
  !/*---------------------------------------
  !|   C = A * D
  !|--------------------------------------*/
  do i=1,mata%n!loop rows
     kr => mata%pa(i)
     ki => mata%pj(i)
     do iisub=1,nsubrow
        do k=1,mata%nnzrow(i)!loop cols
           do kksub=1,nsubrow
              kr%submat(iisub,kksub,k) = kr%submat(iisub,kksub,k) * diag(kksub,ki%cols(k))
           end do
        end do
     end do
  end do

!do i=1,mata%n
!   do j=1,mata%nnzrow(i)!loop cols
!      write(61,*) i,mata%pj(i)%cols(j),((mata%pa(i)%submat(iisub,jjsub,j),iisub=1,mata%nsubrow),jjsub=1,mata%nsubcol)
!   end do
!end do

  return! 0;
end subroutine coscalC
!/*---------------end of coscalC-----------------------------------------
!----------------------------------------------------------------------*/










!void setup_arms (arms Levmat) 
subroutine setup_arms (Levmat) 
  type(arms_data_type),pointer::Levmat
  !allocate(Levmat%ilus,STAT=istat)!replace with createILUT ?
  !allocate(Levmat%levmat,STAT=istat)
  return
end subroutine setup_arms
!|--------------------------------------------------------------------*/



subroutine cleanARMS(ArmsPre) 
  type(arms_data_type),pointer::ArmsPre
  type(p4_type),pointer:: amat,levc,levn
  type(ilut_type),pointer:: cmat
  integer:: indic,istat
  amat=>ArmsPre%levmat
  cmat=>ArmsPre%ilus
  !/*----------------------------------------------------------------------
  !| Free up memory allocated for entire ARMS preconditioner.
  !|----------------------------------------------------------------------
  !| on entry:
  !|==========
  !| ( amat )  =  Pointer to a Per4Mat struct.
  !| ( cmat )  =  Pointer to a IluSpar struct.
  !|--------------------------------------------------------------------*/
  !/* case when nlev == 0 */  
  !int indic=(amat->nB != 0) ;
  indic=0
  if(amat%nB .ne. 0) indic=1
  !    /*  && amat->next !=NULL) ; */
  allocate(levc,STAT=istat)
  allocate(levn,STAT=istat)
  levc => amat 
  if (indic.eq.1) then 
    !while (levc)
     !if (cleanP4(levc)) return(1);
     call cleanP4(levc) 
     levn = levc%next
     nullify(levc);
     levc => levn	
  else 	
     if (associated(amat)) then
        call cleanP4(amat)
        nullify(amat)
     end if
  end if
  call cleanILUT(cmat,indic); 

  !if (associated(cmat)) then
  !  cleanP4(cmat);	
  !  nullify(cmat)
  !end if
  deallocate(ArmsPre,STAT=istat)
  nullify(ArmsPre)
  return! 0;
end subroutine cleanARMS
!/*---------------------------------------------------------------------
!|     end of cleanARMS 
!|--------------------------------------------------------------------*/


!static int arms_free_vcsr(parms_Operator *self)
subroutine arms_free_vcsr(self)
  type(parms_Operator),pointer::self
  type(arms_data_type),pointer :: data!again why cant I just pass self
  data => self%data
  call cleanARMS(data)
  return!
end subroutine arms_free_vcsr
!|---------------------------------------------------------------------*/



subroutine parms_arms_nnz(self, nnz_mat, nnz_pc)
  type(parms_Operator),target::self
  type(arms_data_type),pointer::data
  integer,pointer::nnz_mat,nnz_pc
  data = self%data
  nnz_mat => data%nnz_mat
  nnz_pc  => data%nnz_prec
end subroutine parms_arms_nnz
!|---------------------------------------------------------------------*/



subroutine parms_arms_lsol_vcsr(self, y, x)
  !
  type(parms_Operator)::self
  real(dp),pointer,dimension(:,:)::x,y
  type(arms_data_type),pointer:: data
  type(p4_type),pointer:: levmat
  type(ilut_type),pointer:: ilus
  integer ::schur_start, nlev,i,iisub

  data => self%data
  levmat => data%levmat
  ilus   => data%ilus
  schur_start = data%schur_start

  if (data%ipar(1) .eq. 0) then
    call Lsolp(schur_start, ilus%L, y, x); 
    return 
  end if
  nlev = data%nlev
  do i=1,data%n
     do iisub=1,data%nsubrow
        x(iisub,i)=y(iisub,i)
     end do
  end do
  levmat =>  Lvsol2(x, nlev, levmat, ilus, 0)!ceb
  return !0;
end subroutine parms_arms_lsol_vcsr
!|---------------------------------------------------------------------*/




!|---------------------------------------------------------------------*/
subroutine parms_arms_sol_vcsr(op, y, x) 
  type(parms_Operator),pointer::op
  real(dp),pointer,dimension(:,:)::x,y
  integer i,j,iisub
  type(arms_data_type),pointer:: data
  data => op%data
!write(6,*)"parms_arms_sol_vcsr: checkpt 1"!ceb
  !set x=b

  do i=1,data%n
     do iisub=1,data%nsubrow
        x(iisub,i)=y(iisub,i)
     end do
  end do

!write(6,*)"parms_arms_sol_vcsr: checkpt 3"!ceb
  call armsol2(x, data)
!write(6,*)"parms_arms_sol_vcsr: checkpt 4"!ceb
  return!
end subroutine parms_arms_sol_vcsr
!|---------------------------------------------------------------------*/




subroutine parms_arms_invs_vcsr(self,  y, x) 
  !
  type(parms_Operator),target:: self
  real(dp),pointer,dimension(:,:)::x,y
  type(arms_data_type),pointer:: data
  type(ilut_type),pointer:: ilus
  integer:: i,iisub,schur_start
  data => self%data
  schur_start = data%schur_start
  ilus   => data%ilus
  !if (.not.allocated(y) .or. .not.allocated(x)) return
  if (data%ipar(1) .eq. 0) then
     call invsp(schur_start, ilus, y, x)
  else 
     do i=1,ilus%n
        do iisub=1,data%nsubrow
           x(iisub,i)=y(iisub,i)
        end do
     end do
     call SchLsol(ilus, x);
     call SchUsol(ilus, x);
  end if
  return!
end subroutine parms_arms_invs_vcsr
  !    |---------------------------------------------------------------------*/


!static int parms_arms_ascend_vcsr(parms_Operator self, FLOAT *y, FLOAT *x)  
subroutine parms_arms_ascend_vcsr(self, y, x)  
  !
  type(parms_Operator),target:: self
  real(dp),pointer,dimension(:,:)::x,y

  type(arms_data_type),pointer:: data
  type(p4_type),pointer:: levmat,last
  type(ilut_type),pointer:: ilus
  integer :: schur_start, nloc, lenB, first
  real(dp),pointer,dimension(:,:):: shift_x

  nullify(levmat)
  nullify(last)

  data => self%data;
  levmat => data%levmat
  ilus   => data%ilus
  schur_start = data%schur_start;
  
  if (data%ipar(1) .eq. 0) then
    call Usolp(schur_start, ilus%U, y, x)
    return!
  end if

  do !while (associated(levmat)) 
     if(.not.associated(levmat)) goto 4000
     last => levmat
     levmat => levmat%next
  enddo
4000 continue

  levmat => last
  nloc=levmat%n
  lenB=levmat%nB
  first = data%n - nloc 

  !/*-------------------- last level                                 */
  first = first + lenB 

  !/*-------------------- other levels                               */
  do 
     if(.not.associated(levmat)) goto 5000
     nloc = levmat%n 
     first = first - levmat%nB;
     shift_x => x(1:,first:)
     if (levmat%n.gt.0) call ascend(levmat, shift_x,shift_x)! need to pass x,y at particular index
     !/*-------------------- right scaling */
     if (associated(levmat%D2)) then!or if allocated? 
write(6,*) "parms_arms_ascend_vcsr: calling dscale"
        call dscale(nloc, levmat%nsubrow, levmat%D2, shift_x, shift_x)! need to pass x,y at particular index
     end if
     levmat => levmat%prev
  end do
5000 continue
  return! 0;
end subroutine parms_arms_ascend_vcsr
  !    |---------------------------------------------------------------------*/


 
integer function parms_arms_getssize_vcsr(self)
  type(parms_Operator),target:: self
  !type(arms_data_type),pointer:: data
  !data => self%data
  parms_arms_getssize_vcsr=self%data%schur_start
  return 
end function parms_arms_getssize_vcsr
  !    |---------------------------------------------------------------------*/


subroutine parms_OperatorCreate(self)
  type(parms_Operator),pointer::self
  integer :: istat
  if(.not.associated(self)) allocate(self,STAT=istat)!not sure why we have to create newOpt and point to it
  !allocate(newOpt)
  !newOpt%ref = 1
  !PARMS_NEW0((newOpt)->ops);
  !newOpt%ops%apply  = 0
  !newOpt%ops%lsol   = 0
  !newOpt%ops%invs   = 0
  !newOpt%ops%ascend = 0
  !create an arms_data_type
  call create_arms_data_type(self%data)
  return!
end subroutine parms_OperatorCreate
  !    |---------------------------------------------------------------------*/

subroutine print_parms_Operator(op)
  type(parms_Operator),pointer::op
write(6,*)"parms_Operator_data"
  if(associated(op)) then
     call print_arms_data_type(op%data)
  end if

  return
end subroutine print_parms_Operator

  !    |---------------------------------------------------------------------*/


subroutine create_arms_data_type(data)
  type(arms_data_type),pointer :: data
  integer::istat
  allocate(data,STAT=istat)
  data%n=0
  data%nsubrow=0
  data%nsubcol=0
  data%nlev=0
  data%schur_start=0
  data%ind=0
  data%nnz_mat=0
  data%nnz_prec=0
  nullify(data%ilus)
  nullify(data%levmat)
end subroutine create_arms_data_type
  !    |---------------------------------------------------------------------*/
subroutine print_arms_data_type(data)
  type(arms_data_type),pointer :: data
  if(associated(data)) then
     write(6,*)"arms_data_type"
     write(6,*)"n=",data%n
     write(6,*)"nsubrow=",data%nsubrow
     write(6,*)"nsubcol=",data%nsubcol
     write(6,*)"nlev=",data%nlev
     write(6,*)"schur_start=",data%schur_start
     write(6,*)"ind=",data%ind
     write(6,*)"nnz_mat=",data%nnz_mat
     write(6,*)"nnz_prec=",data%nnz_prec
     !if(associated(data%ilus))!print_ILUT_data(data%ilus)
     !if(associated(data%levmat))!print_dataP4_data(data%levmat)
  end if
end subroutine print_arms_data_type
  !    |---------------------------------------------------------------------*/





subroutine create_FactParam(param, &
                            l_fil,u_fil,eu_fil,lf_fil,s_fil,sl_fil,su_fil, &
                            droptol,n,nsubrow,nsubcol,nlev,bsize,nkrylov_dirs,mreord, &
                            rscal,cscal,rscal_S,cscal_S,pq_S)
  !need to set this up to read from file
  type(parms_FactParam),intent(inout) :: param
  integer,intent(in)::n,nsubrow,nsubcol,nlev,bsize,nkrylov_dirs,mreord,rscal,cscal,rscal_S,cscal_S
  integer,intent(in)::l_fil,u_fil,eu_fil,lf_fil,s_fil,sl_fil,su_fil,pq_S
  real(dp),intent(in)::droptol

  integer::i

  !allocate(param)
  param%n=n
  param%nsubrow=nsubrow
  param%nsubcol=nsubcol
  param%start=1
  param%schur_start=n
  param%mc=0!multicoloring

!write(6,*)"lfil=",lfil
  param%lfil(1)= l_fil !amount of fill-in kept  in L [B].
  param%lfil(2)= u_fil !amount of fill-in kept  in U [B].
  param%lfil(3)=eu_fil !amount of fill-in kept  in E U^-1 
  param%lfil(4)=lf_fil !amount of fill-in kept  in L^-1 F 
  param%lfil(5)= s_fil !amount of fill-in kept  in S     
  param%lfil(6)=sl_fil !amount of fill-in kept  in L [S]  
  param%lfil(7)=su_fil !amount of fill-in kept  in U [S]  

  param%ipar(1)= nlev!number of levels (reduction processes)
  param%ipar(2)= mreord!level-reordering option to be used.
  !                if ipar[2]==0 ARMS uses a block independent set ordering
  !                  with a sort of reverse cutill Mc Kee ordering to build 
  !                  the blocks. This yields a symmetric ordering. 
  !                  in this case the reordering routine called is indsetC
  !                if ipar[2] == 1, then a nonsymmetric ordering is used.
  !                  In this case, the B matrix is constructed to be as
  !                  diagonally dominant as possible and as sparse as possble.
  !                  in this case the reordering routine called is ddPQ.
  param%ipar(3)= bsize ! minimum size of independent set blocks.(not valid if using non symmetric permutations:mreord=1 )
  !param%ipar(4)=iout !if (iout > 0) statistics on the run are printed to FILE *ft
  param%ipar(5)= nkrylov_dirs  !Krylov subspace dimension for last level 
  !                 ipar[5] == 0 means only backward/forward solve is performed on last level. 
  !param%ipar(6)=! maximum # iterations on last level
  param%ipar(7)= 0!ipar[7-10] NOT used 
  param%ipar(8)= 0
  param%ipar(9)= 0
  param%ipar(10)=0
  !       ipar[11-14] == meth[1:4] = method flags for interlevel blocks
  !       ipar[15-18] == meth[1:4] = method flags for last level block -
  !         with the following meaning  
  !           meth[1] permutations of rows  0:no 1: yes. affects rperm
  !             NOT USED IN THIS VERSION ** enter 0.. Data: rperm
  !           meth[2] permutations of columns 0:no 1: yes. So far this is 
  !             USED ONLY FOR LAST BLOCK [ILUTP instead of ILUT]. (so
  !             ipar[12] does no matter - enter zero). If ipar[16] is one
  !             then ILUTP will be used instead of ILUT. Permutation data
  !             stored in: perm2.  
  !           meth[3] diag. row scaling. 0:no 1:yes. Data: D1
  !           meth[4] diag. column scaling. 0:no 1:yes. Data: D2
  !             similarly for meth[15], ..., meth[18] all transformations
  !             related to parametres in meth[*] (permutation,
  !             scaling,..) are applied to the matrix before processing it   

  !all zeros means no permutations and no scaling

  !parameters relating to independent set
  param%ipar(11)=0!pq perm (not used)
  param%ipar(12)=0!not used
  param%ipar(13)=rscal!meth[3] diag. row scaling. 0:no 1:yes. Data: D1
  param%ipar(14)=cscal!cscal!meth[4] diag. column scaling. 0:no 1:yes. Data: D2

  !parameters relating to schur complement
  param%ipar(15)=pq_S !pqperm (second permutation performed on schur complement)
  param%ipar(16)=0!not used
  param%ipar(17)=rscal_S!meth[3] diag. row scaling. 0:no 1:yes. Data: D1
  param%ipar(18)=cscal_S!meth[4] diag. column scaling. 0:no 1:yes. Data: D2

  ! droptol[1] = threshold for dropping  in L [B]. See piluNEW.c:
  ! droptol[2] = threshold for dropping  in U [B].
  ! droptol[3] = threshold for dropping  in L^{-1} F 
  ! droptol[4] = threshold for dropping  in E U^{-1} 
  ! droptol[5] = threshold for dropping  in Schur complement
  ! droptol[6] = threshold for dropping  in L in last block [see ilutpC.c]  
  ! droptol[7] = threshold for dropping  in U in last block [see ilutpC.c] 
!write(6,*)"droptol=",droptol

  do i=1,7
     param%droptol(i)=droptol !currently set up to have one value for all matrices
  end do

  !param%pfgpar(1)=0!?
  !param%pfgpar(2)=0

  !    tolind   = tolerance parameter used by the indset function. 
  !      a row is not accepted into the independent set if the *relative*
  !      diagonal tolerance is below tolind. see indset function for
  !      details. Good values are between 0.05 and 0.5 -- larger values
  !      tend to be better for harder problems.
  !param%tolind=0.45
  param%tolind=0.3
  param%isalloc=1
  return
end subroutine create_FactParam
  !    |---------------------------------------------------------------------*/







subroutine create_parms_Mat(Mat)
  type(parms_Mat)::Mat
  !integer:: ref
  !!  parms_Mat_ops ops;
  !!  void        *data;
  !integer :: isserial;!BOOL
  !integer :: issetup;!BOOL
  !integer :: isperm;!BOOL
  !integer :: isalloc;!BOOL
  !!MATTYPE     type;
  !!PCILUTYPE   ilutype;
  !integer :: m
  !integer :: n
  !integer :: MM
  !integer :: NN
  !type(parms_Map):: is
  call create_parms_Map(Mat%is)
end subroutine create_parms_Mat
  !    |---------------------------------------------------------------------*/

subroutine create_parms_Map(Map)
  type(parms_Map)::Map
  !integer :: ref
  !parms_Table table;		!/**< contains pair (gindex,lindex) */
  !MPI_Comm    comm;		!/**< MPI communicator */
  !integer         pid;		!/**< processor ID */
  !integer         npro;		!/**< number of processors */
  !Map%lsize=!integer         lsize;		!/**< number of local variables */
  !integer         gsize;		!/**< number of global variables */
  !/*! numbering style used.
  ! *
  ! *  \f{tabular}{lcl} 
  ! *   FORTRAN &-& 1 \\
  ! *   C       &-& 0 
  ! *   \f}
  ! */
  Map%start=1!integer :: start
  !   The number of variables associated with each vertex.
  !Map%dof=!integer :: dof!if square submat this would be nsubrow	
  !/*! style of labelling variables u_i,v_i associated with the i-th vertex.
  ! * - NONINTERLACED u_1,v_1,u_2,v_2.
  ! * - INTERLACED u_1,u_2,\cdots,v_1,v_2.
  !VARSTYPE    vtype;	
  !/*! array of size lsize, stores local variables in global labels */
  !integer,dimension(:) :: lvars		
  Map%isserial=1!integer:: isserial	!BOOL	!//!< data layout on one processor?
  !/*! variables are permuted or not.
  ! *  -true: interior variables labelled first followed by interface variables */
  !Map%isperm= !integer:: isperm !BOOL		
  !Map%isvecperm= !integer:: isvecperm !BOOL	!/* used to check if vector has been permuted in gmres */	
  !Map%ispermalloc=0 !integer:: ispermalloc !BOOL	!/* permutation array allocated? */
  !Map%perm !integer,dimension(:),allocatable :: perm		!//!< permutation array of size lsize.
  !Map%iperm !integer,dimension(:),allocatable :: iperm		!//!< inverse permutation array. 
  !Map%nint= !integer:: nint		!//!< number of interior variables.
  !Map%ninf= !integer:: ninf		!//!< number of interface variables.
  !  \param schur_start  start of the local Schur complement which
  !  may be lesser than nbnd for the matrix with unsymmetric pattern 
  Map%schur_start=1 !integer:: schur_start	

  !integer:: ninf_send!	//!< number of variables to be sent.
  ! param vsend  array of size ninf_send which stores variables in local indices 
  !integer,dimension(:) :: vsend		
  ! param vstable  stores variables to be sent to adjacent processors in pairs
  ! 	 (local_index, index_in_vsend) 
  !parms_Table vstable
end subroutine create_parms_Map
  !    |---------------------------------------------------------------------*/

end module armslib
  !    |---------------------------------------------------------------------*/
  !    |---------------------------------------------------------------------*/
