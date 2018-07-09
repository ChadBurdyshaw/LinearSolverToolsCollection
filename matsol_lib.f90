module matsol_lib
  use kinddefs
  use blas
contains
!#include <stdio.h>
!#include <stdlib.h>
!#if defined(C99)
!#include <tgmath.h>
!#else
!#include <math.h>
!#endif 
!#include "protos.h"

!subroutine diag_scal( vbsptr vbmat )
!!{
!!  !----------------------------------------------------------------------------
!!   * Diagonal scaling:
!!   * For the matrix with block diagonals D1, D2, ..., Dp :
!!   *       D1 x  x  x
!!   *       x  D2 x  x
!!   *  A =  x  x ... x
!!   *       x  x  x  Dp
!!   * simply take the block diagonal matrix
!!   *       D1 0  0  0
!!   *       0  D2 0  0
!!   *  D =  0  0 ... 0
!!   *       0  0  0  Dp
!!   * invert D and do A := inv(D)*A
!!   * then the diagonal blocks of A are now identities.
!!   *----------------------------------------------------------------------------
!!   * Parameters
!!   *----------------------------------------------------------------------------
!!   * on entry:
!!   * =========
!!   * vbmat    = block matrix stored in VBSparRow format -- see globheads.h for
!!   *            details on format, the block sizes might be different
!!   * on return:
!!   * ==========
!!   * vbmat    = inv(D)*vbmat
!!   *--------------------------------------------------------------------------*/
!  integer i, j, k, dim, sz, size, ierr = 0, col 
!  real(dp) one=1.0, zero=0.0  
!  integer nnzrow, n = vbmat%n, *bsz = vbmat%bsz, *ja 
!  integer bufsz = sizeof(real(dp))*MAX_BLOCK_SIZE*MAX_BLOCK_SIZE 
!  BData *ba, *D, buf 
!  D = (BData *)Malloc( sizeof(BData)*n, "diag_scal" ) 
!  buf = (BData)Malloc( bufsz, "diag_scal" ) 
!  for( i = 0  i < n  i++ ) {
!    nnzrow = vbmat%nnzrow(i) 
!    ja = vbmat%ja(i) 
!    for( j = 0  j < nnzrow  j++ ) {
!      if( ja(j) .ne. i ) continue 
!      dim = B_DIM( bsz, i ) 
!      size = sizeof(real(dp))*dim*dim 
!      D(i) = (BData)Malloc( size, "diag_scal" ) 
!      memcpy(D(i), vbmat%ba(i)(j), size ) 
!      ierr = invSVD( dim, D(i) ) 
!      if( ierr .ne. 0 ) {
!	for( k = 0  k < i  k++ ) free( D(k) ) 
!	free( D ) 
!	fprintf( stderr, "error: Singular diagonal block...\n" ) 
!	return -2 
!      }
!    }
!  }
!  for( i = 0  i < n  i++ ) {
!    dim = B_DIM( bsz, i ) 
!    nnzrow = vbmat%nnzrow(i) 
!    ja = vbmat%ja(i) 
!    ba = vbmat%ba(i) 
!    for( j = 0  j < nnzrow  j++ ) {
!      col = ja(j) 
!      sz = B_DIM( bsz, col ) 
!      GGEMM ("n","n", dim, sz, dim, one, D(i), dim, ba(j), dim,
!	     zero,buf, dim)      
!      copyBData( dim, sz, ba(j), buf, 0 ) 
!    }
!  }
!  vbmat%D = D 
!  free( buf ) 
!  return 0 
!}

!integer diagvec( vbsptr vbmat, BData x, BData y )
!!{
!  !---------------------------------------------------------------------
!!    | This function does y = inv(D) x, where D is diagonals of A.
!!    |----------------------------------------------------------------------
!!    | on entry:
!!    | vbmat = the matrix (in BSparRow form)
!!    | x     = a vector
!!    |
!!    | on return
!!    | y     = the product inv(D) * x
!!    |--------------------------------------------------------------------*/
!  integer i, n = vbmat%n, *bsz = vbmat%bsz, dim, sz = 1 
!  real(dp) zero=0.0, one = 1.0 
!  BData *D = vbmat%D 
!  do (i = 0  i < n  i++ ) {
!    dim = B_DIM( bsz, i ) 
!    DGEMM ("n","n", dim, sz, dim, one,D(i), dim,x+bsz(i), dim,
!	   zero, y+bsz(i), dim)    
!  }
!  return 0 
!!}


subroutine matvec(mata, x, y)  
  !---------------------------------------------------------------------
  !  | This function does the matrix vector product y = A x.
  !  |----------------------------------------------------------------------
  !  | on entry:
  !  | mata  = the matrix (in SparRow form)
  !  | x     = a vector
  !  |
  !  | on return
  !  | y     = the product A * x
  !  |--------------------------------------------------------------------*/
  type(cs_type),pointer::mata
  real(dp),pointer,dimension(:,:)::x,y
  !   local variables    */
  integer :: i, k,iisub,kksub
  type(ja_type),pointer::ki 
  type(submat_type),pointer::kr  
  do i=1,mata%n
     kr => mata%pa(i) 
     ki => mata%pj(i) 
     do iisub=1,mata%nsubrow
        y(iisub,i) = 0.0 
        do k=1,mata%nnzrow(i)
           do kksub=1,mata%nsubcol
              !write(6,*)"matvec: kr%submat(",iisub,kksub,k,")=",kr%submat(iisub,kksub,k)
              !write(6,*)"matvec: x(",kksub,ki%cols(k),")=",x(kksub,ki%cols(k))
              y(iisub,i) = y(iisub,i) + kr%submat(iisub,kksub,k) * x(kksub,ki%cols(k))
           end do
        end do
!write(6,*)"matvec: y(",iisub,i,")=",y(iisub,i)
     end do
  end do
  nullify(kr)
  nullify(ki)
  return 
end subroutine matvec
!----------------------------------------------------------------------*/

subroutine matvec_prll(pcomm,mata, x, y)  
  !---------------------------------------------------------------------
  !  | This function does the matrix vector product y = A x.
  !  |----------------------------------------------------------------------
  !  | on entry:
  !  | mata  = the matrix (in SparRow form)
  !  | x     = a vector
  !  |
  !  | on return
  !  | y     = the product A * x
  !  |--------------------------------------------------------------------*/
  !type(PComm_type),intent(inout)::pcomm
  integer :: pcomm
  type(cs_type),pointer::mata
  real(dp),pointer,dimension(:,:)::x,y
  !   local variables    */
  integer :: i, k,iisub,kksub
  type(ja_type),pointer::ki 
  type(submat_type),pointer::kr  
  do i=1,mata%n
     kr => mata%pa(i) 
     ki => mata%pj(i) 
     do iisub=1,mata%nsubrow
        y(iisub,i) = 0.0 
        do k=1,mata%nnzrow(i)
           do kksub=1,mata%nsubcol
              !write(6,*)"matvec: kr%submat(",iisub,kksub,k,")=",kr%submat(iisub,kksub,k)
              !write(6,*)"matvec: x(",kksub,ki%cols(k),")=",x(kksub,ki%cols(k))
              y(iisub,i) = y(iisub,i) + kr%submat(iisub,kksub,k) * x(kksub,ki%cols(k))
           end do
        end do
!write(6,*)"matvec: y(",iisub,i,")=",y(iisub,i)
     end do
  end do


  !Aswap()

  nullify(kr)
  nullify(ki)
  return 
end subroutine matvec_prll
!----------------------------------------------------------------------*/


!subroutine vbmatvec(vbmat, x, y)
!  !-------------------- matrix -- vector product in VB format */
!  type(vbs_type),pointer::vbmat
!  real(dp),dimension(:,:)::x,y
!  !   local variables    */
!  integer ::i, j,z,iisub,jjsub,nnzrow,col,inc,dim,sz,nBs,nBsj,n
!  integer,pointer,dimension(:)::bsz
!  type(ja_type),pointer::ja 
!  type(submat_type),pointer::ba
!  real(dp)::one
!  inc = 1 
!  n = vbmat%n
!  bsz => vbmat%bsz 
!  one=1.0 
!  !BData *ba 
!  do i = 1,n
!     nBs = bsz(i) 
!     dim = B_DIM(bsz,i) 
!     do j = 1,dim 
!        y(nBs+j) = 0 
!     end do
!     nnzrow = vbmat%nnzrow(i) 
!     ja => vbmat%ja(i) 
!     ba => vbmat%ba(i) 
!     do j = 1,nnzrow
!        col = ja(j) 
!        nBsj = bsz(col) 
!        sz = B_DIM(bsz,col) 
!        !-------------------- operation:  y = Block*x + y */
!        call DGEMV ("n", dim, sz, one, ba(j), dim, &x(nBsj), inc, one, &y(nBs), inc) 
!     end do
!  end do
!end subroutine vbmatvec



subroutine invsp(start, ilusch, y, x)
  !---------------------------------------------------------------------
  !  |
  !  | This routine does the backward solve U x = y, where U is upper or
  !  | bottom part of local upper triangular matrix
  !  |
  !  | Can be done in place.
  !  |
  !  | Zhongze Li, Aug. 17th, 2001
  !  |
  !  |----------------------------------------------------------------------
  !  | on entry:
  !  | start = the index of the first component
  !  | n     = one ore than the index owned by the processor 
  !  | y     = a vector
  !  | ilusch  = the LU matrix as provided from the ILU routines.
  !  |
  !  | on return
  !  | x     = the product U^{-1} * y
  !  |
  !  |---------------------------------------------------------------------*/
  integer :: start   !ceb start may need to be zero else 1 
! probably need start to be zero, depends on how start is passed in as to what to do here
  type(ilut_type),pointer::ilusch
  real(dp),pointer,dimension(:,:)::x,y
  !   local variables    */
  integer ::i, k,iisub,kksub,n
  type(ja_type),pointer::ki 
  type(submat_type),pointer::kr
  real(dp),dimension(ilusch%nsubrow)::t

  n = ilusch%L%n 

  do i=start,n 
     do iisub=1,ilusch%nsubrow
        t(iisub) = y(iisub,i-start) 
        if ( ilusch%L%nnzrow(i) .gt. 0 ) then
           kr => ilusch%L%pa(i) 
           ki => ilusch%L%pj(i) 
           do k=1,ilusch%L%nnzrow(i)
              !if (ki%cols(k) .ge. start .and. ki%cols(k) .lt. n) then
              if (ki%cols(k) .ge. start .and. ki%cols(k) .lt. n) then !should this be ki(k).le.n
                 do kksub=1,ilusch%nsubcol
                    t(iisub) = t(iisub)- kr%submat(iisub,kksub,k)*y(kksub,ki%cols(k)-start) 
                 end do
              end if
           end do
        end if
        x(iisub,i-start) = t(iisub)
     end do
  end do

  !ceb  do i=n-1  i>= start  i-- 
  do i=n,start,-1 
     kr => ilusch%U%pa(i) 
     ki => ilusch%U%pj(i) 
     do iisub=1,ilusch%nsubrow
        t(iisub)  = x(iisub,i-start) 
        !do (k=1  k<ilusch%U%nnzrow(i)  k++)
        do k=2,ilusch%U%nnzrow(i)
           do kksub=1,ilusch%nsubcol
              t(iisub) =t(iisub)- kr%submat(iisub,kksub,k) * x(kksub,ki%cols(k)-start) 
           end do
        end do
        x(iisub,i-start) = t(iisub)*kr%submat(iisub,iisub,1) 
     end do
  end do
  
  nullify(ki)
  nullify(kr)
  return
end subroutine invsp
!----------------------------------------------------------------------*/



subroutine arms_Lsol(mata, b, x)
  !---------------------------------------------------------------------
  !  | This function does the forward solve L x = b.
  !  | Can be done in place.
  !  |----------------------------------------------------------------------
  !  | on entry:
  !  | mata  = the matrix (in SparRow form)
  !  | b     = a vector
  !  |
  !  | on return
  !  | x     = the solution of L x = b 
  !  |--------------------------------------------------------------------*/
  type(cs_type),pointer::mata
  real(dp),pointer,dimension(:,:):: x,b
  !   local variables    */
  integer i, k,iisub,kksub 
  type(ja_type),pointer::ki 
  type(submat_type),pointer::kr
  !write(6,*)"arms_Lsol: checkpt 0"
  do i=1,mata%n
     kr => mata%pa(i) 
     ki => mata%pj(i) 
     do iisub=1,mata%nsubrow
        !write(6,*)"1: arms_Lsol: init b(",iisub,i,")=",b(iisub,i)
        x(iisub,i) = b(iisub,i) 
        !ceb cols(k) should always be < i (if L rows are ordered) so x does not have to be initialized here
        do k=1,mata%nnzrow(i)
           !if(ki%cols(k)>i) write(6,*) "arms_Lsol: error:k:",ki%cols(k)," > i:",i !ceb
           do kksub=1,mata%nsubcol
              !write(6,*)"2:   arms_Lsol: x(",kksub,ki%cols(k),")=",x(kksub,ki%cols(k))
              !write(6,*)"2.5: arms_Lsol: kr%submat(",iisub,kksub,k,")=",kr%submat(iisub,kksub,k)
              x(iisub,i) = x(iisub,i) - kr%submat(iisub,kksub,k)*x(kksub,ki%cols(k))
           end do
        end do
     end do
     !do iisub=1,mata%nsubrow
     !write(6,*)"3: arms_Lsol: final x(",iisub,i,")=",x(iisub,i)
     !end do
  end do
  !write(6,*)"arms_Lsol: checkpt 10"
  nullify(ki)
  nullify(kr)
  return
end subroutine arms_Lsol
!---------------end of arms_Lsol-----------------------------------------
!----------------------------------------------------------------------*/



subroutine arms_Usol(mata, b, x)
  !---------------------------------------------------------------------
  !  | This function does the backward solve U x = b.
  !  | Can be done in place.
  !  |----------------------------------------------------------------------
  !  | on entry:
  !  | mata  = the matrix (in SparRow form)
  !  | b    = a vector
  !  |
  !  | on return
  !  | x     = the solution of U * x = b 
  !  |
  !  |---------------------------------------------------------------------*/
  type(cs_type),pointer::mata
  real(dp),pointer,dimension(:,:):: x,b
  !   local variables    */
  integer i, j,iisub,jjsub 
  real(dp)::t(mata%nsubrow)
  type(ja_type),pointer::ki 
  type(submat_type),pointer::kr
!write(6,*)"arms_Usol: checkpt 0"

  !do i=1,mata%n
  !   do iisub=1,mata%nsubrow
  !      write(6,*)"1: arms_Usol: init b(",iisub,i,")=",b(iisub,i)
  !   end do
  !end do

!write(6,*)""
!write(6,*)"arms_Usol: checkpt 1: levmat%U%n=",mata%n," size(x)=",size(x)/mata%nsubrow," size(b)=",size(b)/mata%nsubrow
  do i=mata%n,1,-1
     kr => mata%pa(i) 
     ki => mata%pj(i) 
     do iisub=1,mata%nsubrow
        x(iisub,i) = b(iisub,i) 
!write(6,*)"1.5:   arms_Usol: init: x(",iisub,i,")=",x(iisub,i)
        do j=2,mata%nnzrow(i)!loop off diagonals
           do jjsub=1,mata%nsubcol
!write(6,*)"2:   arms_Usol: x(",jjsub,ki%cols(j),")=",x(jjsub,ki%cols(j))
!write(6,*)"2.5: arms_Usol: kr%submat(",iisub,jjsub,j,")=",kr%submat(iisub,jjsub,j)
              x(iisub,i) = x(iisub,i) - kr%submat(iisub,jjsub,j) * x(jjsub,ki%cols(j)) 
           end do
        end do
        !x(iisub,i) = x(iisub,i) * kr%submat(iisub,iisub,1) 
!write(6,*)"arms_Usol: x(",iisub,i,")=",x(iisub,i)
     end do
     !//x[i]=D[i]*x[i]
     do iisub=1,mata%nsubrow
        t(iisub)=0.0
        do jjsub=1,mata%nsubcol		
           t(iisub) = t(iisub) + kr%submat(iisub,jjsub,1)*x(jjsub,i)
        end do
     end do

     do iisub=1,mata%nsubrow
        x(iisub,i)=t(iisub)
!write(6,*)"arms_Usol:final  x(",iisub,i,")=",x(iisub,i)
     end do
  end do

!write(6,*)"arms_Usol: checkpt 10"
  nullify(ki)
  nullify(kr)
  return
end subroutine arms_Usol
!----------------end of arms_Usol--------------------------------------
!----------------------------------------------------------------------



subroutine descend(levmat, x, wk)
  !---------------------------------------------------------------------
  !| This function does the (block) forward elimination in ARMS
  !|                       new       old
  !|     |            |  |     |    |    |
  !|     | L        0 |  | wx1 |    | x1 |
  !|     |            |  |     | =  |    | 
  !|     | EU^{-1}  I |  | wx2 |    | x2 |
  !|     |            |  |     |    |    |
  !| x used and not touched -- or can be the same as wk.
  !|--------------------------------------------------------------------
  type(p4_type),pointer::levmat
  real(dp),pointer,dimension(:,:):: x,wk
  !  local variables   
  integer ::j,jjsub,lenB,permj
  real(dp),pointer,dimension(:,:):: shift_work,shift_wk

  lenB=levmat%nB
  !------------------------------------------------------
  !|   apply row permutation P to rhs 
  !|-----------------------------------------------------
  do j=1,levmat%n
     permj=j !row permutation
     if(associated(levmat%rperm)) permj=levmat%rperm(j)
     do jjsub=1,levmat%nsubrow
        levmat%wk(jjsub,j) = x(jjsub,permj)  
        !write(6,*)"descend: permuted b(",jjsub,j,")=",levmat%wk(jjsub,j)
     end do
  end do
  
  !write(6,*)"  descend: checkpt 1: solve [L]g_p=g"!ceb   g_p(nB) = [L^-1] g
  call arms_Lsol(levmat%L, levmat%wk, wk)  
  
  !write(6,*)"  descend: checkpt 2: solve [U] h_p=g_p"!ceb  wk(nB)= [U^-1]g_p
  call arms_Usol(levmat%U, wk, levmat%wk)  
  !-------------------- compute x(lenb:.) = x (lenb:.) - E * work(1) */
  ! i.e  h_p(nC) = h - E[U^-1] g_p
  shift_work => levmat%wk(1:,lenB+1:)  
  shift_wk => wk(1:,lenB+1:)  
  !write(6,*)"  descend: checkpt 3: solve h_p(nC) = h - E[U^-1] f_p"!ceb  h_p(nC) = h - E[U^-1] f_p
  call matvecz (levmat%E, levmat%wk, shift_work, shift_wk)
  
  !do j=1,levmat%n
  !   do jjsub=1,levmat%nsubrow
  !      write(6,*)"descend: final h_p(",jjsub,j,")=",wk(jjsub,j)
  !   end do
  !end do

  !write(6,*)"descend: checkpt 10"!ceb
  nullify(shift_work)
  nullify(shift_wk)
  return 
end subroutine descend
!----end-of-descend-----------------------------------------------------
!|----------------------------------------------------------------------




!|----------------------------------------------------------------------
subroutine ascend (levmat, x, wk) 
  !---------------------------------------------------------------------
  !  | This function does the (block) backward substitution: 
  !  |
  !  |     |            |  |     |    |    |
  !  |     | U  L^{-1}F |  | wk1 |    | x1 |
  !  |     |            |  |     | =  |    |
  !  |     | 0       S  |  | wk2 |    | x2 |  <<-- x2 already computed.
  !  |     |            |  |     |    |    |       and we need x1
  !  |
  !  |    with x2 = S^{-1} wk2 (assumed to have been computed ) 
  !  |-------------------------------------------------------------------
  type(p4_type),pointer::levmat
  real(dp),pointer,dimension(:,:):: x,wk
  !--------------------  local variables  */
  integer ::j,jjsub,len,lenB
  integer,dimension(:),pointer::qperm
  real(dp),pointer,dimension(:,:):: shift_x
  real(dp),pointer,dimension(:,:):: shift_wk
!write(6,*)"ascend: checkpt 0"!ceb

  len=levmat%n
  lenB=levmat%nB
  qperm=>levmat%perm 
  !-------------------- copy x onto wk */  
  if(lenB.ne.len) then
!     write(6,*)"  ascend: checkpt 1: x2 start=",lenB+ofst!ceb
     shift_x=>x(1:,lenB+1:)
     !write(6,*) "compute z_p= [F]z"
     call matvec(levmat%F, shift_x, levmat%wk)! work = F * x_2 
     !write(6,*) "compute z_pp=[L^-1]z_p"
     call arms_Lsol(levmat%L, levmat%wk, levmat%wk)!  work = L \ work
!write(6,*) "compute g_pp=g_p-[L^-1][F]z"
     do j=1,lenB               !  wk1 = wk1 - work 
        do jjsub=1,levmat%nsubrow   
           levmat%wk(jjsub,j) = x(jjsub,j) - levmat%wk(jjsub,j)
           !write(6,*)"ascend: gpp(",jjsub,j,")=",levmat%wk(jjsub,j)         
        end do
     end do
  end if

!write(6,*)"Back sub: y=[U^-1]g_pp"
  call arms_Usol(levmat%U, levmat%wk, levmat%wk)!wk1 = U \ wk1 

!ceb load nC portion of wK with solution vector
  !memcpy(&work(lenB),&x(lenB),(len-lenB)*sizeof(real(dp)))
  !do j=0,(len-lenB)-1              
  do j=1,(len-lenB)             
     do jjsub=1,levmat%nsubrow   
        levmat%wk(jjsub,lenB+j) = x(jjsub,lenB+j)
     end do
  end do

  !---------------------------------------
  !  |   apply reverse permutation
  !  |------------------------------------
  do j=1,len
     do jjsub=1,levmat%nsubrow
        wk(jjsub,qperm(j)) = levmat%wk(jjsub,j) 
     end do
  end do

!write(6,*)"ascend: checkpt 10"!ceb
   
  nullify(qperm)
  nullify(shift_x)
  return
end subroutine ascend
!----end-of-ascend----------------------------------------------------
!|--------------------------------------------------------------------
!|--------------------------------------------------------------------




subroutine matvecz(mata, x, y, z) 
  !---------------------------------------------------------------------
  !  | This function does the matrix vector  z = y - A x.
  !  |------------------------------------------------------------------
  !  | on entry:
  !  | mata  = the matrix (in SparRow form)
  !  | x, y   = two input vector
  !  |
  !  | on return
  !  | z    = the result:  y - A * x
  !  | z-location must be different from that of x 
  !  | i.e., y and x are used but not modified.
  !  |-------------------------------------------------------------------
  !   local variables
  type(cs_type),pointer:: mata
  real(dp),pointer,dimension(:,:)::x,y,z
  type(ja_type),pointer :: ki
  type(submat_type),pointer::kr
  integer i, k,iisub,kksub
  real(dp),dimension(mata%nsubrow)::t
!write(6,*)"matvecz: checkpt 0"!ceb
  do i=1,mata%n
     kr => mata%pa(i) 
     ki => mata%pj(i) 
     do iisub=1,mata%nsubrow
        t(iisub) = y(iisub,i)  
        do k=1,mata%nnzrow(i)
           do kksub=1,mata%nsubcol
              t(iisub) =t(iisub)- kr%submat(iisub,kksub,k) * x(kksub,ki%cols(k))
           end do
        end do
        z(iisub,i) = t(iisub)
!        write(6,*)"matvecz: z(",iisub,i,")=",z(iisub,i)
     end do
  end do
  nullify(ki)
  nullify(kr)
  return
end subroutine matvecz
!---------------end of matvecz----------------------------------------
! *-------------------------------------------------------------------


function Lvsol2(x, nlev, levmat, ilusch, flag) result (returnmat)
  ! Macro L-solve -- corresponds to left (L) part of arms
  !   |  preconditioning operation -- 
  !   |  on entry : 
  !   |   x =  right- hand side to be operated on by the preconditioner
  !   |  on return : x is overwritten
  !   |   x =  output result of operation 
  !   |  
  !   |  Note : in-place operation -- b and x can occupy the same space..
  !   | --------------------------------------------------------------------*
  type(p4_type),pointer:: levmat
  type(ilut_type),pointer:: ilusch
  real(dp),pointer,dimension(:,:)::x
  integer ::nlev,flag,ilev
  real(dp),pointer,dimension(:,:):: shift_x
  type(p4_type),pointer::returnmat
  !-------------------- local variables  
  integer ::first,lenB,lenC,nloc
  type(p4_type),pointer:: last!=levmat  

  nloc=levmat%n
  lenB =levmat%nB 
  last=>levmat
  !-------------------- take care of  special cases :  nlev==0 -% lusol  
!write(6,*)"Lvsol2: checkpt 0"!ceb
  first = 1!ceb 
  !-------------------- descend                                      
  !while (levmat) { 
!write(6,*)"  Lvsol2: checkpt 5.5"!ceb
  ilev=0
  do!loop over all nlev
     if (.not.associated(levmat) .or. lenB.eq.0) exit 
!write(6,*)"  Lvsol2: checkpt 6 ilev=",ilev!ceb
     nloc=levmat%n 
     lenB =levmat%nB 

     !-------------------- left scaling  
!write(6,*)"  6.5 Lvsol2:  first=",first,"  lenB=",lenB
     shift_x=>x(1:,first:)
     if (associated(levmat%D1)) then
!write(6,*)"Lvsol2: calling row scale for levmat"
        call dscale(nloc,levmat%nsubrow,levmat%D1, shift_x, shift_x)
     end if

     !--------------------  RESTRICTION/ DESCENT OPERATION  */
!write(6,*)"  Lvsol2: checkpt 7 "!ceb
     if (lenB.ne.0) call descend(levmat, shift_x, shift_x) !ceb

     if(lenB.eq.nloc) then
        !write(6,*)"exiting Lvsol2: no schur complement"
        returnmat=>last
        return
     end if
     ilev=ilev+1

!write(6,*)"  Lvsol2: checkpt 8: going to next depth level: ilev=",ilev!ceb
     if(ilev.le.nlev) first = first + lenB  
     last => levmat 

     levmat => levmat%next !get next level B,F

     !---------------------------------------------------------------------
     !  | next level 
     !--------------------------------------------------------------------
  end do



  !write(6,*)"Lvsol2:checkpt 9: at deepest level ",ilev,": first=",first                                
  shift_x=>x(1:,first:)
  !write(6,*)"  Lvsol2: checkpt 9.5:
  if (flag.ne.0) call SchLsol(ilusch,shift_x) !ceb
  !write(6,*)"Lvsol2: checkpt 10"!ceb

  returnmat=>last
  nullify(last)
  nullify(shift_x)
  return
end function Lvsol2
!|---------------------------------------------------------------------



subroutine Uvsol2(x, nlev, n, levmat, ilusch) 
  ! Macro U-solve -- corresponds to left (L) part of arms
  !   |  preconditioning operation -- 
  !   |  on entry : 
  !   |  b  =  right- hand side to be operated on by the preconditioner
  !   |  on return  = x has been overwritten =
  !   |  x  =  output result of operation 
  !   |  
  !   |  Note : in-place operation -- b and x  can occupy the same space..
  !   | --------------------------------------------------------------------*
  type(p4_type),pointer::levmat
  type(ilut_type),pointer::ilusch
  real(dp),pointer,dimension(:,:) ::x
  integer :: nlev,n,ilev
  !-------------------- local variables  
  real(dp),pointer,dimension(:,:):: shift_x
  integer :: nloc, lenB, first  
  !-------------------- work array                                        
  !-------------------- take care of  special cases :  nlev==0 --> lusol  
  !-------------------- case of zero levels                             
!write(6,*)"Uvsol2: checkpt 0"!ceb
  !-------------------- general case                               
  nloc=levmat%n  
  lenB=levmat%nB  
  first = n - nloc + 1 !ceb 
!write(6,*)"Uvsol2: checkpt 3 n=",n," nloc=",nloc," n-nloc= first(init)=",first," lenB=",lenB
  !-------------------- last level                                 
  first = first + lenB  
!write(6,*)"  Uvsol2: checkpt 4: lenB=",lenB," first_init+lenB=first=",first!ceb
  shift_x=>x(1:,first:)
  if (lenB.ne.nloc .or.lenB.eq.0) then
!write(6,*)"  Uvsol2: checkpt 4.5: schur U solve"
     call SchUsol(ilusch, shift_x)  !ceb
  end if
  !-------------------- other levels                               
  ilev=nlev!-1
  do
!write(6,*)"  Uvsol2: checkpt 5 ilev=",ilev!ceb
     if(.not.associated(levmat) .or. levmat%nB.eq.0) exit

     nloc = levmat%n  
     first = first - levmat%nB 
     ilev=ilev-1
!write(6,*)"Uvsol2: checkpt 6 first=",first," nB=",levmat%nB!ceb
     shift_x=>x(1:,first:)
     if (levmat%n.ne.0) call ascend(levmat, shift_x,shift_x) 
!write(6,*)"Uvsol2: checkpt 7"!ceb
     !-------------------- right scaling 
     if (associated(levmat%D2)) then
!write(6,*) "Uvsol2: col scaling for levmat"
        call dscale(nloc,levmat%nsubrow,levmat%D2,shift_x,shift_x) 
     end if

!write(6,*)"Uvsol2: checkpt 9"!ceb
     levmat => levmat%prev  !get previous B,F block
  end do

!write(6,*)"Uvsol2: checkpt 10"!ceb
  nullify(shift_x)
  return 
  !--------------------  PROLONGATION/ ASCENT OPERATION */
end subroutine Uvsol2
!|---------------------------------------------------------------------







subroutine armsol2(x, Prec) 
  !     |  combined preconditioning operation -- combines the
  !     |  left and right actions. 
  !     | 
  !     |  on entry : 
  !     |   x =  right- hand side to be operated on by the preconditioner
  !     |  on return : x is overwritten - 
  !     |   x =  output result of operation 
  !     |  
  !     |  Note : in-place operation -- b and x can occupy the same space..
  !     | --------------------------------------------------------------------*/ 
  type(arms_data_type),pointer ::Prec
  real(dp),pointer,dimension(:,:)::x
  !  !-------------------- local variables  */
  type(p4_type),pointer:: levmat ! = Prec%levmat 
  type(ilut_type),pointer:: ilusch ! = Prec%ilus 
  integer ::nlev ! = Prec%nlev 
  integer ::n,i,j,test ! = levmat%n  
  type(p4_type),pointer:: last
  real(dp),dimension(:,:),pointer::wk!ceb

  levmat=>Prec%levmat
  ilusch=>Prec%ilus
  nlev=Prec%nlev
  n=levmat%n
  test=1
!write(6,*)""!ceb

  if (ilusch%n>0) then!=================================
!write(6,*)"armsol2: checkpt 0"!ceb
     last => Lvsol2(x, nlev, levmat, ilusch, 1)
     call Uvsol2(x, nlev, n, last, ilusch)
     !do i=1,ilusch%n
     !   do j=1,ilusch%nsubrow
     !      write(6,*)"armsol: final x(",j,i,")=",x(j,i)
     !   end do
     !end do
  end if!===========================================

  nullify(last)
  nullify(ilusch)
  nullify(levmat)
  return! 0  
end subroutine armsol2
!|---------------------------------------------------------------------*/ 





subroutine SchLsol(ilusch, y) 
  !  !---------------------------------------------------------------------
  !    |  Forward solve for Schur complement part = 
  !    |----------------------------------------------------------------------
  !    | on entry:
  !    | ilusch  = the LU matrix as provided from the ILU functions.
  !    | y       = the right-hand-side vector
  !    |
  !    | on return
  !    | y       = solution of LU x = y. (overwritten) 
  !    |---------------------------------------------------------------------*/
  type(ilut_type),pointer :: ilusch
  real(dp),pointer,dimension(:,:):: y
  !!-------------------- local variables                        */
  integer :: n,j,i
  n=ilusch%n

!write(6,*)"SchLsol: checkpt 0"!ceb
  !!-------------------- begin: right scaling                          */
  if (associated(ilusch%D1)) then
!write(6,*) "schLsol: calling schur dscale row scaling"
     call dscale(n,ilusch%nsubrow, ilusch%D1, y, y) 
  end if
  !!-------------------- ONE SIDED ROW PERMS */
  if (associated(ilusch%rperm)) then !reverses ordering
     do j=1,n
!if(j .ne. ilusch%rperm(j))write(6,*)"  SchLsol: row rperm j:",j,"=",ilusch%rperm(j)
        do jjsub=1,ilusch%nsubrow
           ilusch%wk(jjsub,j) = y(jjsub,ilusch%rperm(j))
      !write(6,*)"SchLsol: b(",jjsub,ilusch%rperm(j),")=",ilusch%wk(jjsub,j)
!write(6,*)"SchLsol: b(",jjsub,j,")=",ilusch%wk(jjsub,j)
        end do
     end do
     !!--------------------  L solve proper */
!write(6,*)"  SchLsol: checkpt 5: permuted Lsolve Schur complement "!ceb
     call arms_Lsol(ilusch%L, ilusch%wk, y)  
  else 
!write(6,*)"  SchLsol: checkpt 7: unpermuted Lsolve Schur complement "!ceb
     call arms_Lsol(ilusch%L, y, y) 
  end if

  return
end subroutine SchLsol
!!---------------end of SchLsol---------------------------------------
!  ----------------------------------------------------------------------*/
!|---------------------------------------------------------------------*/ 



subroutine SchUsol(ilusch, y) 
  !---------------------------------------------------------------------
  !  | U-solve for Schur complement  - 
  !  |----------------------------------------------------------------------
  !  | on entry:
  !  | ilusch  = the LU matrix as provided from the ILU functions.
  !  | y       = the right-hand-side vector
  !  |
  !  | on return 
  !  | y       = solution of U x = y. (overwritten on y) 
  !  |----------------------------------------------------------------------*/
  type(ilut_type),pointer::ilusch
  real(dp),pointer,dimension(:,:)::y

  integer j,jjsub,permj

!write(6,*)"SchUsol: checkpt 0"!ceb
  ! -------------------- begin by U-solving */
  !-------------------- CASE: column pivoting  used (as in ILUTP) */
  if (associated(ilusch%perm2))then
     call arms_Usol(ilusch%U, y, y) 
!write(6,*)"  SchUsol: checkpt 1: column pivoting Usolve permute solution to original numbering "!ceb
     do j=1,ilusch%n
!if(j .ne. ilusch%perm2(j))write(6,*)"  SchUsol: local column pivoting j:",j,"=",ilusch%perm2(j)
        permj=ilusch%perm2(j)
        do jjsub=1,ilusch%nsubrow
           ilusch%wk(jjsub,permj) = y(jjsub,j) 
!write(6,*)"1: SchUsol: wk(",jjsub,permj,")=",y(jjsub,j)
        end do
     end do
  else
     !-------------------- CASE: no column pivoting  used                   */
!write(6,*)"  SchUsol: checkpt 3: unpermuted local cols Usolve Schur complement "!ceb
     call arms_Usol(ilusch%U, y, ilusch%wk) 
  end if


  !-------------------- generic permutation                              */
  if (associated(ilusch%perm)) then
!write(6,*)"  SchUsol: checkpt 4: applying generic permutation to z"!ceb
     do j=1,ilusch%n
!if(j .ne. ilusch%perm(j))write(6,*)"  SchUsol: generic permutation (symmetric?) j:",j,"=",ilusch%perm(j)
        do jjsub=1,ilusch%nsubrow
           !y(jjsub,j) = ilusch%wk(jjsub,ilusch%perm(j)) 
           y(jjsub,ilusch%perm(j)) = ilusch%wk(jjsub,j) 
!write(6,*)"2: SchUsol: qr perm z(",jjsub,j,")=",y(jjsub,j)
        end do
     end do
  else
     do j=1,ilusch%n
        do jjsub=1,ilusch%nsubrow
           y(jjsub,j)=ilusch%wk(jjsub,j)
!write(6,*)"3: SchUsol: prescal z(",jjsub,j,")=",ilusch%wk(jjsub,j) 
        end do
     end do
  end if


!write(6,*)"SchUsol: checkpt 8"!ceb
  !-------------------- case when diagonal scaling is done on columns    */ 
  if (associated(ilusch%D2)) then
!write(6,*) "schUsol: calling dscale col scaling for schur complement"
     call dscale(ilusch%n,ilusch%nsubrow,ilusch%D2,y, y) 
  end if

!ceb
  !do j=1,ilusch%n
  !   do jjsub=1,ilusch%nsubrow
  !write(6,*)"SchUsol: final unscaled z(",jjsub,j,")=",y(jjsub,j)         
  !   end do
  !end do

!write(6,*)"SchUsol: checkpt 10"!ceb
  return
end subroutine SchUsol
!---------------end of SchUsol---------------------------------------
!----------------------------------------------------------------------*/





!subroutine luinv(n, nsub,a, x, y)
!  !--------------------------------------------------------
!  ! *    does the operation y = inv(a) * x
!  ! *    where a has already been factored by Gauss.
!  ! *    LUy = x
!  ! *------------------------------------------------------*/
!  integer ::n,nsub
!  real(dp),dimension(:,:)::a,x,y
!  integer ::i, j,iisub,jjsub
!  integer:: bsA, bsB!block offsets 
!  real(dp),dimension(:) ::sum 
!  ! Ly0 = x -- use Lsol ? */   
!  do i = 1,n 
!     do iisub=1,nsub 
!        sum(iisub) = x(iisub,i) 
!        bsA = i - n 
!        do  j = 1,i 
!           do jjsub=1,nsub
!              bsA = bsA + n 
!              sum(iisub) = sum(iisub) - a(jjsub,bsA) * y(jjsub,j)  ! a(i,j) * y(j) */
!           end do
!        end do
!        y(iisub,i) = sum(iisub) 
!     end do
!  end do
!
!  ! Uy = y0 */
!  bsB = i * n 
!  do i = n,1,-1 
!     do iisub=1,nsub 
!        sum(iisub) = y(iisub,i) 
!        bsB = bsB - n 
!        bsA = i+bsB 
!        do j = i+1,n 
!           do jjsub=1,nsub
!              bsA = bsA + n 
!              sum(iisub) = sum(iisub) - a(jjsub,bsA) * y(jjsub,j)  ! a(i,j) * y(j) */
!           end do
!        end do
!        y(iisub,i) = sum(iisub) * a(iisub,bsB+i)  ! a(i,i) */
!     end do
!  end do
!end subroutine luinv




subroutine lusolC(y, x, lu)
  !----------------------------------------------------------------------
  ! *    performs a forward followed by a backward solve
  ! *    for LU matrix as produced by iluk
  ! *    y  = right-hand-side 
  ! *    x  = solution on return 
  ! *    lu = LU matrix as produced by iluk. 
  ! *--------------------------------------------------------------------*/
  type(ilu_type),pointer::lu
  real(dp),pointer,dimension(:,:)::x,y

  integer::n,i,j,iisub,jjsub,nnzrow
  type(ja_type),pointer::ja
  type(submat_type),pointer::pa
  type(cs_type),pointer::L,U
  real(dp),pointer,dimension(:,:)::D 

  n = lu%n
  L => lu%L 
  U => lu%U 
  D => lu%D 

  ! Block L solve */
  do i = 1,n 
     nnzrow = L%nnzrow(i) 
     ja => L%pj(i) 
     do iisub=1,lu%nsubrow
        x(iisub,i) = y(iisub,i) 
        do j = 1,nnzrow
           do jjsub=1,lu%nsubcol        
              x(iisub,i) = x(iisub,i) - x(jjsub,ja%cols(j)) * L%pa(i)%submat(iisub,jjsub,j) 
           end do
        end do
     end do
  end do
  ! Block -- U solve */
  do i = n,1,-1
     nnzrow = U%nnzrow(i)
     ja = U%pj(i) 
     do iisub=1,lu%nsubrow
        do  j = 1,nnzrow
           do jjsub=1,lu%nsubcol    
              x(iisub,i) = x(iisub,i) - x(jjsub,ja%cols(j)) * U%pa(i)%submat(iisub,jjsub,j) 
           end do
        end do
        x(iisub,i) = x(iisub,i) * D(iisub,i)
     end do
  end do

  nullify(L)
  nullify(U)
  nullify(D)
  nullify(pa)
  nullify(ja)
  return   
end subroutine lusolC
!|---------------------------------------------------------------------*/ 




!integer vblusolC( real(dp) *y, real(dp) *x, vbiluptr lu)  
!!{
!  !----------------------------------------------------------------------
!  ! *    performs a forward followed by a backward block solve
!  ! *    for LU matrix as produced by VBILUT
!  ! *    y  = right-hand-side 
!  ! *    x  = solution on return 
!  ! *    lu = LU matrix as produced by VBILUT
!  ! *
!  ! *    note: lu%bf is used to store vector
!  ! *--------------------------------------------------------------------*/
!  integer n = lu%n, *bsz = lu%bsz, i, j, bi, icol, dim, sz 
!  integer nnzrow, nBs, nID, *ja, inc = 1, OPT 
!  real(dp) *data, alpha = -1.0, beta = 1.0, alpha2 = 1.0, beta2 = 0.0 
!  vbsptr L, U 
!  BData *D, *ba 
!
!  L = lu%L 
!  U = lu%U 
!  D = lu%D 
!  OPT = lu%DiagOpt 
!  ! Block L solve */
!  for( i = 0  i < n  i++ ) {
!    dim = B_DIM(bsz,i) 
!    nBs = bsz(i) 
!    for( j = 0  j < dim  j++ ) {
!      nID = nBs + j 
!      x(nID) = y(nID) 
!    }
!
!    nnzrow = L%nnzrow(i) 
!    ja = L%ja(i) 
!    ba = L%ba(i) 
!    for( j = 0  j < nnzrow  j++ ) {
!      icol = ja(j) 
!      sz = B_DIM(bsz,icol) 
!      data = ba(j) 
!      DGEMV( "n",  dim,  sz, alpha, data, dim, x+bsz(icol), inc,
!	     beta, x+nBs, inc )  
!    }
!  }
!  ! Block -- U solve */
!  for( i = n-1  i >= 0  i-- ) {
!    dim = B_DIM(bsz,i) 
!    nnzrow = U%nnzrow(i) 
!    nBs = bsz(i) 
!    ja = U%ja(i) 
!    ba = U%ba(i) 
!    for( j = 0  j < nnzrow  j++ ) {
!      icol = ja(j) 
!      sz = B_DIM(bsz,icol) 
!      data = ba(j) 
!      DGEMV( "n", dim, sz, alpha, data, dim, x+bsz(icol), inc,
!	     beta, x+nBs, inc ) 
!    }
!    data = D(i) 
!    if (OPT == 1) 
!      luinv( dim, data, x+nBs, lu%bf ) 
!    else
!      DGEMV( "n", dim, dim, alpha2, data, dim, x+nBs, inc,
!	     beta2, lu%bf, inc ) 
!	
!    for( bi = 0  bi < dim  bi++ ) {
!      x(nBs+bi) = lu%bf(bi) 
!    }
!  }
!
!  return 0 
!!}



subroutine lumsolC(y,x, lu )
  !----------------------------------------------------------------------
  ! *    performs a forward followed by a backward solve
  ! *    for LU matrix as produced by iluc
  ! *    y  = right-hand-side
  ! *    x  = solution on return
  ! *    lu = LU matrix as produced by iluc.
  ! *--------------------------------------------------------------------*/
  real(dp),pointer,dimension(:,:)::x,y
  type(ilu_type),pointer::lu

  type(cs_type),pointer::L,U
  integer ::n, i, j, nnzrow, nnzL
  type(ja_type),pointer::ia,ja
  type(submat_type),pointer::ma
  real(dp),pointer,dimension(:,:)::D

  n = lu%n
  D => lu%D
  L => lu%L 
  U => lu%U 
  
  do i = 1,lu%n
     do iisub=1,lu%nsubrow    
        x(iisub,i) = y(iisub,i) 
     end do
  end do
  !-------------------- L solve */
  do i = 1,lu%n
     nnzL = L%nnzrow(i) 
     ia => L%pj(i) 
     ma => L%pa(i)
     do iisub = 1,lu%nsubrow 
        do j = 1,nnzL
           do jjsub=1,lu%nsubcol
              x(jjsub,ia%cols(j)) =  x(jjsub,ia%cols(j)) - ma%submat(iisub,jjsub,j) * x(iisub,i) 
           end do
        end do
     end do
  end do
  !-------------------- U solve */
  !do i = lu%n-1  i >= 0  i-- ) {
  do i = lu%n,1,-1
     nnzrow = U%nnzrow(i) 
     ja => U%pj(i) 
     ma => U%pa(i)
     do iisub =1,lu%nsubrow 
        do j = 1,nnzrow
           do jjsub=1,lu%nsubcol
              x(iisub,i) = x(iisub,i) - ma%submat(iisub,jjsub,j) * x(jjsub,ja%cols(j))
           end do
        end do
        x(iisub,i) = x(iisub,i) * D(iisub,i) 
     end do
  end do

  nullify(D)
  nullify(L)
  nullify(U)
  nullify(ia)
  nullify(ja)
  nullify(ma)

  return
end subroutine lumsolC
!|---------------------------------------------------------------------*/ 



subroutine rpermC(mat, perm)
  !----------------------------------------------------------------------
  !  |
  !  | This subroutine permutes the rows of a matrix in SparRow format. 
  !  | rperm  computes B = P A  where P is a permutation matrix.  
  !  | The permutation P is defined through the array perm: for each j, 
  !  | perm(j) represents the destination row number of row number j. 
  !  |
  !  |-----------------------------------------------------------------------
  !  | on entry:
  !  |----------
  !  | (amat) = a matrix stored in SparRow format.
  !  |
  !  |
  !  | on return:
  !  | ----------
  !  | (amat) = P A stored in SparRow format.
  !  |
  !  | integer value returned:
  !  |             0   -% successful return.
  !  |             1   -% memory allocation error.
  !  |---------------------------------------------------------------------*/
  type(cs_type),pointer:: mat
  integer,dimension(:)::perm
  type(ja_ptr),target,dimension(:),allocatable::addj
  type(submat_ptr),target,dimension(:),allocatable::addm
  integer,dimension(:),allocatable::nnz
  integer ::i,istat,permi
  type(ja_type),dimension(:),allocatable,target::row_pj
  type(submat_type),dimension(:),allocatable,target::row_pa

  allocate(addj(mat%n),STAT=istat)
  allocate(addm(mat%n),STAT=istat)
  allocate(nnz(mat%n),STAT=istat)
  allocate(row_pj(mat%n),STAT=istat)
  allocate(row_pa(mat%n),STAT=istat)

write(6,*) "calling rpermC"

  !ceb this function present a problem in that it is not clear how to use fortran pointers to swap
  !    rows for the object mat which is a pointer but whose variables are neither pointer nor target
  !since we cant just point mat%pj(i) to aja (pointer can't point to pointer) 
  ! i guess will just have to copy element by element, which means that we'll have 
  ! to allocate buffer memory 
  ! unfortunately we are not just swapping rows but rearranging them so we cant just store one
  !   row at a time

  !addj,addm points to old rows
  do i=1,mat%n
     !permi=perm(i) 
     addj(i)%ptr => mat%pj(i) 
     addm(i)%ptr => mat%pa(i) 
     nnz( i)     =  mat%nnzrow(i)
  end do

  !assign swapped rows to mat
  do i=1,mat%n
     permi=perm(i) 
     !check to see if current row needs to be swapped
     if (i .ne. permi) then

        !allocate memory to store row that is being replaced
        allocate(row_pj(i)%cols(mat%nnzrow(i)),STAT=istat)
        allocate(row_pa(i)%submat(mat%nsubrow,mat%nsubcol,mat%nnzrow(i)),STAT=istat)

        !copy former mat(i) data to new memory location:  row(i)=matrow(permi)
        do j=1,nnz(i)
           row_pj(i)%cols(j)=mat%pj(i)%cols(j)
           do iisub=1,mat%nsubrow
              do jjsub=1,mat%nsubcol
                 row_pa(i)%submat(iisub,jjsub,j)=mat%pa(i)%submat(iisub,jjsub,j)
              end do
           end do
        end do

        !point add(perm(i)) to recently allocated  memory location so we dont stomp old data in m
        addj(i)%ptr => row_pj(i) !still points to matrow(i) data
        addm(i)%ptr => row_pa(i)  


        ! need to resize row of mat to accept new data: wipe matrow(i)
        if(mat%nnzrow(i) .ne. nnz(permi)) then
           if(allocated(mat%pj(i)%cols)) deallocate(mat%pj(i)%cols)
           if(allocated(mat%pa(i)%submat)) deallocate(mat%pa(i)%submat)
           allocate(mat%pj(i)%cols(nnz(permi)))
           allocate(mat%pa(i)%submat(mat%nsubrow,mat%nsubcol,nnz(permi)))
        end if       



!write(6,*)"swapping row ",i," with row ",perm(i)

        !assign: matrow(i)=matrow(iperm(i))
        !copy data from old row location to new row location
        do j=1,nnz(permi)
           mat%pj(i)%cols(j) = addj(permi)%ptr%cols(j) !mat%pj(i) => addj(i)%ptr =mat(permi)
           do iisub=1,mat%nsubrow
              do jjsub=1,mat%nsubcol
                 mat%pa(i)%submat(iisub,jjsub,j) = addm(permi)%ptr%submat(iisub,jjsub,j)!mat%pa(i) => addm(i)%ptr 
              end do
           end do
        end do
        mat%nnzrow(i) = nnz(permi)   



 
        !check to see which row buffers can be deallocated
        !do j=1,mat%n
        !   if (perm(j) .le. i) then
        !      if(allocated(row_pj(perm(j))%cols)) then
        !         !write(6,*)"eliminating temp row ",perm(j)
        !         deallocate(row_pj(perm(j))%cols)
        !      end if
        !      if(allocated(row_pa(perm(j))%submat)) deallocate(row_pa(perm(j))%submat)
        !   end if
        !end do

     end if
  end do

  !clean pointers
  do i=1,mat%n
     nullify(addj(i)%ptr)
     nullify(addm(i)%ptr)
  end do
  deallocate(row_pj)
  deallocate(row_pa)
  deallocate(addj) 
  deallocate(addm) 
  deallocate(nnz) 
  return 
end subroutine rpermC
!|---------------------------------------------------------------------*/ 



subroutine cpermC(mat, perm, flag) 
  !----------------------------------------------------------------------
  !  |
  !  | This subroutine permutes the columns of a matrix in SparRow format.
  !  | cperm computes B = A P, where P is a permutation matrix.
  !  | that maps column j into column perm(j), i.e., on return 
  !  | The permutation P is defined through the array perm: for each j, 
  !  | perm(j) represents the destination column number of column number j. 
  !  |
  !  |-----------------------------------------------------------------------
  !  | on entry:
  !  |----------
  !  | (mat) = a matrix stored in SparRow format.
  !  |
  !  |
  !  | on return:
  !  | ----------
  !  | (mat) = A P stored in SparRow format.
  !  |
  !  | integer value returned:
  !  |             0   -% successful return.
  !  |             1   -% memory allocation error.
  !  |---------------------------------------------------------------------*/
  type(cs_type),pointer :: mat
  integer,intent(in),dimension(:) :: perm  
  integer,intent(in),optional::flag
  !internal vars
  integer :: istat,size,i,j,k,ii,jj,pos,permj,itmp,newj_nnz,found,test,currnnzrow,newjk,n,nsubrow,nsubcol
  integer,dimension(:),allocatable::newj,rvperm
  integer,pointer,dimension(:)::ja
  real(dp),pointer,dimension(:,:,:)::ma
  !real(dp),dimension(:,:,:),allocatable::new_a
  real(dp),dimension(:,:,:),pointer::new_a
  allocate(rvperm(mat%n),STAT=istat)
  allocate(newj(mat%n),STAT=istat)
  allocate(new_a(mat%nsubrow,mat%nsubcol,mat%n),STAT=istat)
  n=mat%n
  nsubrow=mat%nsubrow
  nsubcol=mat%nsubcol
  !actually only needs to be size nnzrow(i) 
  ! but would have to be resized of each row. mat%n is always >= nnzrow
  write(6,*) "calling cpermC"
  !if(allocated(perm)) then
  do i=1,mat%n
     rvperm(perm(i))=i
  end do
  if (.not.present(flag)) then ! do simple column reordering
     do i=1,n
        do j=1,mat%nnzrow(i)
           newj(j)=perm(mat%pj(i)%cols(j))
           mat%pj(i)%cols(j)=newj(j)!ceb
        end do
     end do
  else ! reorder column indices and matrix elements
     do i=1,n !loop over all rows
        currnnzrow=mat%nnzrow(i)
        ja=>mat%pj(i)%cols
        ma=>mat%pa(i)%submat
        !for this row, create new reordered ja list, adding only nonzero elements
        newj_nnz=0
        do k=1,n
           do j=1,currnnzrow
              !if(ja(j)==perm(k) .and. newj_nnz<currnnzrow) then
              if(ja(j).eq.perm(k)) then
                 newj_nnz=newj_nnz+1
                 newj(newj_nnz)=perm(k)
                 !write(6,*)"i:",i," newja(",newj_nnz,")=",newj(newj_nnz)
                 !temporarily store col entries in new order 
                 do ii=1,nsubrow
                    do jj=1,nsubcol
                       new_a(ii,jj,newj_nnz)=ma(ii,jj,j)
                    end do
                 end do

              end if
           end do
        end do
        !copy reordered submatrix elements back into new location
        do j=1,currnnzrow
           newj(j)=rvperm(newj(j))
           !if(ja(j) .ne. newj(j)) then
              ja(j)=newj(j)
              !write(6,*)"row ",i," final: new ja(",j,")=",mat%pj(i)%cols(j)
              do ii=1,nsubrow
                 do jj=1,nsubcol
                    ma(ii,jj,j)=new_a(ii,jj,j)
                 end do
              end do
           !end if
        end do
     end do
  end if
  nullify(ja)
  nullify(ma)
  deallocate(newj) 
  return 
end subroutine cpermC
!|---------------------------------------------------------------------*/ 



subroutine cpermC2(mat, perm, flag) 
  !----------------------------------------------------------------------
  !  |
  !  | This subroutine permutes the columns of a matrix in SparRow format.
  !  | cperm computes B = A P, where P is a permutation matrix.
  !  | that maps column j into column perm(j), i.e., on return 
  !  | The permutation P is defined through the array perm: for each j, 
  !  | perm(j) represents the destination column number of column number j. 
  !  |
  !  |-----------------------------------------------------------------------
  !  | on entry:
  !  |----------
  !  | (mat) = a matrix stored in SparRow format.
  !  |
  !  |
  !  | on return:
  !  | ----------
  !  | (mat) = A P stored in SparRow format.
  !  |
  !  | integer value returned:
  !  |             0   -% successful return.
  !  |             1   -% memory allocation error.
  !  |---------------------------------------------------------------------*/
  type(cs_type),pointer :: mat
  integer,intent(in),dimension(:) :: perm  
  integer,intent(in),optional::flag
  !internal vars
  integer :: istat,size,i,j,k,ii,jj,pos,permj,itmp,newj_nnz,found,test
  integer,dimension(:),allocatable::newj,rvperm
  integer,pointer,dimension(:)::ja
  real(dp),dimension(:,:,:),allocatable::new_a

  allocate(rvperm(mat%n),STAT=istat)
  allocate(newj(mat%n),STAT=istat)
  allocate(new_a(mat%nsubrow,mat%nsubcol,mat%n),STAT=istat)

  !actually only needs to be size nnzrow(i) 
  ! but would have to be resized of each row. mat%n is always >= nnzrow
  write(6,*) "calling cpermC"
  !if(allocated(perm)) then

       do i=1,mat%n
          rvperm(perm(i))=i
       end do

  if (.not.present(flag)) then ! do simple column reordering

     do i=1,mat%n
        do j=1,mat%nnzrow(i)
           newj(j)=perm(mat%pj(i)%cols(j))
        end do
        do j=1,mat%nnzrow(i)
           mat%pj(i)%cols(j)=newj(j)
        end do
     end do

  else ! reorder column indices and matrix elements

     do i=1,mat%n
        
        ja=>mat%pj(i)%cols
   
        do j=1,mat%n
           newj(j)=perm(j)
!write(6,*)"i:",i," newja(",j,")=",newj(j)
        end do

        newj_nnz=0
        ja_indx=0
        k=0
        do  !k=1,mat%n
           if(k.ge.mat%n .or. newj_nnz.ge.mat%nnzrow(i)) exit
           k=k+1

           found=0
           do j=1,mat%nnzrow(i)
!write(6,*)"1:  i:",i," ja(",j,")=",ja(j)," newj(",k,")=",newj(k)

              if(ja(j).eq. newj(k)) then
!write(6,*)"1:  i:",i," new ja(",j,")=",newj(k)
                 newj_nnz=newj_nnz+1 
                 found=1
                 ja_indx=j
                 !add to new_a
                 do ii=1,mat%nsubrow
                    do jj=1,mat%nsubcol
                       new_a(ii,jj,newj_nnz)=mat%pa(i)%submat(ii,jj,ja_indx)
                    end do
                 end do
                 exit !exit ja search and go to next newj element
              end if

           end do
          
           if (found.eq.0) then
!write(6,*)"shifting "
              !element was not found in ja so we eliminate it from newj
              do j=newj_nnz+1,mat%n-1
!write(6,*)"shifting newj(",j,")=",newj(j+1)
                 newj(j)=newj(j+1)
!write(6,*)"2:  i:",i," newj(",j,")=",newj(k)
              end do
              k=k-1
           end if

        end do



        do j=1,mat%nnzrow(i)
           newj(j)=rvperm(newj(j))
           mat%pj(i)%cols(j)=newj(j)
           !write(6,*)"row ",i," final: new ja(",j,")=",mat%pj(i)%cols(j)

           do ii=1,mat%nsubrow
              do jj=1,mat%nsubcol
                 mat%pa(i)%submat(ii,jj,j)=new_a(ii,jj,j)
              end do
           end do

        end do



        !ceb sort
        !itmp=1
        !call qsort_col(mat%pa(i)%submat,mat%pj(i)%cols,mat%nsubrow,itmp,mat%nnzrow(i))!ceb
     end do
  end if

  nullify(ja)
  deallocate(newj) 
  return 
end subroutine cpermC2
!|---------------------------------------------------------------------*/ 



!void qsort2C(int *ja, double *ma, int left, int right, int abval){
recursive subroutine qsort_col(ma,ja,nsub,left,right)
!/*----------------------------------------------------------------------
!|
!| qqsort: sort ma[left]...ma[right] into increasing order
!| from Kernighan & Ritchie
!|
!| ja holds the column indices
!| abval = 1: consider absolute values
!|         0: values
!|
!|---------------------------------------------------------------------*/
  integer,intent(inout),dimension(:)::ja
  real(dp),intent(inout),dimension(:,:,:)::ma
  integer,intent(inout) :: left,right
  integer,intent(in)::nsub
  !internal vars
  integer :: i, last,itmp
  if (left >= right)  return
  itmp=(left+right)/2
  call swapj(ja, left, itmp)
  last = left
  do i=left+1,right
     if (ja(i) < ja(left)) then
        last=last+1
        call swapj(ja, last, i)
        !call swap_submat(ma,nsub,last,i)
     end if
  end do
  call swapj(ja, left, last)
  !call swap_submat(ma,nsub,left,last)
  itmp=last-1
  call qsort_col(ma,ja,nsub,left,itmp)
  itmp=last+1
  call qsort_col(ma,ja,nsub,itmp,right)
  return
end subroutine qsort_col



subroutine dpermC(mat, perm) 
  !----------------------------------------------------------------------
  !  |
  !  | This subroutine permutes the rows and columns of a matrix in 
  !  | SparRow format.  dperm computes B = P^T A P, where P is a permutation 
  !  | matrix.
  !  |
  !  |-----------------------------------------------------------------------
  !  | on entry:
  !  |----------
  !  | (amat) = a matrix stored in SparRow format.
  !  |
  !  |
  !  | on return:
  !  | ----------
  !  | (amat) = P^T A P stored in SparRow format.
  !  |
  !  | integer value returned:
  !  |             0   -% successful return.
  !  |             1   -% memory allocation error.
  !  |---------------------------------------------------------------------*/
  type(cs_type),pointer::mat
  integer,dimension(:)::perm
  !if (
  call rpermC(mat, perm)!) return 1 
  !if (
  call cpermC(mat, perm)!) return 1 
  return !0 
end subroutine dpermC
!|---------------------------------------------------------------------*/ 



!subroutine CSparTran( csptr amat, csptr bmat, CompressType *compress )
subroutine CSparTran(amat,bmat)
  !----------------------------------------------------------------------
  !  | Finds the compressed transpose of a matrix stored in SparRow format.
  !  | Patterns only.
  !  |-----------------------------------------------------------------------
  !  | on entry:
  !  |----------
  !  | (amat)     = a matrix stored in SparRow format.
  !  | (compress) = quotient graph of matrix amat
  !  |
  !  | on return:
  !  | ----------
  !  | (bmat)     = the compressed transpose of (mata) stored in SparRow
  !  |              format.
  !  |
  !  | integer value returned:
  !  |             0   -% successful return.
  !  |             1   -% memory allocation error.
  !  |---------------------------------------------------------------------*/
  type(cs_type),pointer::amat,bmat
  integer i, j, nnzrow, pos 
  type(ja_type),pointer::aja
  integer,pointer,dimension(:)::ind

  ind => bmat%nnzrow 
  do i=1,amat%n
     ind(i) = 0 
  end do
  !-------------------- compute lengths  */
  do i=1,amat%n 
    !if( compress(i).grp .ne. -1 ) continue 
    aja => amat%pj(i) 
    nnzrow = amat%nnzrow(i) 
    do j=1, nnzrow
      pos = aja%cols(j) 
      !if( compress(pos).grp == -1 ) {
	ind(pos)=ind(pos)+1 
      !}
    end do
  end do
  !--------------------  allocate space  */
  do i=1,amat%n
     if( ind(i) .eq. 0 ) then
        deallocate(bmat%pj(i)%cols) 
        continue !exit do loop?
     end if
     !bmat%pj(i) = (integer *)Malloc( ind(i)*sizeof(int), "CSparTran" ) 
     allocate(bmat%pj(i)%cols(ind(i)))
     ind(i) = 0  ! indicate next available position of each row */
  end do
  !--------------------  now do the actual copying  */
  do i=1,amat%n
    !if( compress(i).grp .ne. -1 ) continue 
    aja => amat%pj(i) 
    nnzrow = amat%nnzrow(i) 
    do j = 1,nnzrow
      pos = aja%cols(j) 
      !if( compress(pos).grp == -1 ) {
	bmat%pj(pos)%cols(ind(pos)) = i 
	ind(pos)=ind(pos)+1 
      !}
    end do
  end do
  nullify(aja)
  nullify(ind)
  return 
end subroutine CSparTran
!|---------------------------------------------------------------------*/ 







subroutine SparTran(amat, bmat, job, flag)
  !/*----------------------------------------------------------------------
  !| Finds the transpose of a matrix stored in SparRow format.
  !|
  !|-----------------------------------------------------------------------
  !| on entry:
  !|----------
  !| (amat) = a matrix stored in SparRow format.
  !|
  !| job    = integer to indicate whether to fill the values (job.eq.1)
  !|          of the matrix (bmat) or only the pattern.
  !|
  !| flag   = integer to indicate whether the matrix has been filled
  !|          0 - no filled
  !|          1 - filled
  !|
  !| on return:
  !| ----------
  !| (bmat) = the transpose of (mata) stored in SparRow format.
  !|
  !| integer value returned:
  !|             0   --> successful return.
  !|             1   --> memory allocation error.
  !|---------------------------------------------------------------------*/
  type(cs_type),pointer ::amat,bmat
  integer:: job,flag
  !internal vars
  integer::  i,j,pos,istat,nsubrow,nsubcol,iisub,jjsub
  !integer,dimension(:),allocatable::iab
  integer,dimension(:),allocatable::ind

  !write(6,*)"spartran size=",amat%n
  nsubrow=amat%nsubrow
  nsubcol=amat%nsubcol
  !allocate(iab(amat%n+1),STAT=istat)
  allocate(ind(amat%n),STAT=istat)

  do i=1,amat%n
     ind(i) = 0
  end do
  !iab(1)=1
  
  if(flag.eq.0) then
     !/*--------------------  compute lengths  */
     do i=1,amat%n
        do j=1,amat%nnzrow(i)
           pos = amat%pj(i)%cols(j)
!write(6,*) "spartran checkpt 1:   ind(",pos,")=",ind(pos),">amat%nnzrow(",i,"):",amat%nnzrow(i)
           ind(pos)=ind(pos)+1
!if(ind(pos)>amat%nnzrow(i)) 
!write(6,*) "spartran checkpt 1.5: iab(",pos,")=",iab(pos),">amat%nnzrow(",i,"):",amat%nnzrow(i)
        end do
     end do

     !do i=1,amat%n
     !   iab(i+1)=iab(i+1)+iab(i)
     !end do

     !/*--------------------  allocate space  */
     do i=1,bmat%n
        !write(6,*)"spartran checkpt 1.6: alloc bmat%pj(",i,")%cols(",ind(i),")"
        bmat%nnzrow(i) = ind(i)!iab(i+1)-iab(i) 
        !write(6,*)"bmat%nnzrow(",i,")",bmat%nnzrow(i)!ceb
        allocate(bmat%pj(i)%cols(bmat%nnzrow(i)),STAT=istat) !this causes memory access errors in ind()
        if (job .eq. 1) allocate(bmat%pa(i)%submat(nsubrow,nsubcol,bmat%nnzrow(i)),STAT=istat)
        !ind(i) = 1!reset ind
     end do
  end if


  do i=1,amat%n
     ind(i)=1
  end do

!write(6,*)"spartran checkpt 2"

  !/*--------------------  now do the actual copying  */
  do i=1,amat%n
!write(6,*)"amat%nnzrow(",i,")=",amat%nnzrow(i)," bmat%nnzrow(",i,")=",bmat%nnzrow(i)
     do j=1,amat%nnzrow(i)
        pos = amat%pj(i)%cols(j)

!write(6,*)"spartran checkpt 3: j=",j," ind(pos:",pos,")=",ind(pos)," bmat%nnzrow(",i,")",bmat%nnzrow(i)
        bmat%pj(pos)%cols(ind(pos)) = i
        if(i==pos) bmat%pd(pos)=ind(pos)!ceb
        if (job .eq. 1) then !copy submat
           do iisub=1,nsubrow
              do jjsub=1,nsubcol
                 !bmat%pa(pos)%submat(iisub,jjsub,ind(pos)) = amat%pa(i)%submat(iisub,jjsub,j)
                 bmat%pa(pos)%submat(iisub,jjsub,ind(pos)) = amat%pa(i)%submat(jjsub,iisub,j)
              end do
           end do
        end if
        ind(pos)=ind(pos)+1
        !if(pos > amat%n .or. pos.eq.0) write(6,*)" index ind(",pos,")=",ind(pos)," i=",i," j=",j," amat%nnzrow(i)",amat%nnzrow(i)
     end do
  end do

!write(6,*)"spartran checkpt 7"
!deallocate(iab,STAT=istat)
!write(6,*)"spartran checkpt 8"
deallocate(ind,STAT=istat)
!write(6,*)"spartran checkpt 9"
  return
end subroutine SparTran
!/*-------------- end of SparTran ---------------------------------------
!|---------------------------------------------------------------------*







!/*----------------------------------------------------------------------
!| Finds the transpose of a matrix stored in SparRow format.
!|-----------------------------------------------------------------------
!serial routie only (n must be equal to max column index (no phantom nodes))
subroutine SparTran_mod(ia,ja,iaT,jaT,n,nnz)
  !
  integer,intent(in):: n,nnz
  integer,intent(in)::ia(n+1),ja(nnz)
  integer,intent(inout)::iaT(n+1),jaT(nnz)
  integer i,j,pos,jstart,jend
  integer :: ind(n)
  !/*--------------------  compute lengths  */

  !in FORTRAN these arrays must be initialized
  do i=1,n
     iaT(i)=0
  end do
  
  !write(6,*) "SparTran chkpt 1"
  iaT(1)=1

  do i=1,n 
     jstart=ia(i)
     jend=ia(i+1)-1
     do j=jstart,jend
        iaT(ja(j)+1) = iaT(ja(j)+1)+1
     end do
  end do
  
  do i=1,n
     iaT(i+1)=iaT(i+1)+iaT(i)
  end do
  
  !write(6,*) "SparTran chkpt 3"
  
  !ceb need to store jcol with appropriate j index  for jaT
  do i=1,n
     !write(6,*)"iaT(",i,")=",iaT(i)
     ind(i)=iaT(i)
  end do
  !write(6,*)"iaT(",n+1,")=",iaT(n+1)
  !write(6,*) "SparTran chkpt 4"
  
  !/*--------------------  now do the actual copying  */
  do i=1,n 
     jstart=ia(i)
     jend=ia(i+1)-1
     do j=jstart,jend
        pos = ja(j)
        !write(6,*) "i=",i," ind(",pos,")=",ind(pos)
        jaT(ind(pos)) = i
        ind(pos) = ind(pos)+1
     end do
  end do
  
  !  do i=1,n 
  !     jstart=iaT(i)
  !     jend=iaT(i+1)-1
  !     do j=jstart,jend
  !          write(6,*) "i=",i," jaT(",j,")=",jaT(j)
  !     end do
  !  end do
  !write(6,*) "SparTran chkpt 5"
end subroutine SparTran_mod
!/*---------------------------------------------------------------------





!/*---------------------------------------------------------------------
!|     defines weights based on diagonal dominance ratios
!|--------------------------------------------------------------------*/
!subroutine weightsC(mat,w,n,nsub)       
subroutine weightsC_mod(A,ia,ja,n,nsub,nnz,w)       
  !
  integer,intent(in):: n,nsub,nnz
  integer,intent(in):: ia(n+1), ja(nnz)
  real(dp),intent(in)::A(nsub,nsub,nnz)
  real(dp),intent(inout)::w(n)
  !internal vars
  real(dp):: tdia, wmax, tnorm,submat_norm
  real(dp) ::tmat(nsub,nsub)
  integer:: irow, k,kstart,kend
  wmax=0.0

  do irow=1,n
     kstart=ia(irow)
     kend=ia(irow+1)-1
     !kz = (kstart-kend+1)!mat->nnzrow(irow)!row nnz array
     !kr = A!mat->pa(irow)!matrix row array
     !kj = ja!mat->pj(irow)!matrix col index array
     tnorm = 0.0
     tdia = 0.0
     do k=kstart,kend
        call get_submat(A,tmat,k,nnz,nsub)
        submat_norm= L1Norm(tmat,nsub)
        if (ja(k) == irow) tdia = submat_norm !if (kj(k) == irow) tdia = abs(kr(k))
        tnorm = tnorm + submat_norm !tnorm = tnorm + abs(kr(k))
     end do
     if (tnorm > 0.0) tnorm =  tdia / tnorm
     w(irow) = tnorm
     if (tnorm > wmax) wmax = tnorm
  end do
  do irow=1,n
     w(irow) = w(irow)/wmax
     !write(6,*)"row ",irow," weight ",w(irow)
  end do
end subroutine weightsC_mod
!|---------------------------------------------------------------------*/ 






subroutine weightsC(mat, w)       
  !/*---------------------------------------------------------------------
  !|     defines weights based on diagonal dominance ratios
  !|--------------------------------------------------------------------*/
  type(cs_type),pointer :: mat
  real(dp),intent(inout),dimension(:)::w !ceb should this be 2d?
  !internal vars
  integer:: irow, k, n,kz,nsubrow
  real(dp)::tdia,wmax,tnorm,kr_norm
  type(ja_type),pointer :: kj
  type(submat_type),pointer::kr

  n=mat%n
  nsubrow=mat%nsubrow
  wmax=0.0

  do irow=1,n
     kz = mat%nnzrow(irow)
     kr => mat%pa(irow)
     kj => mat%pj(irow)

     tnorm = 0.0
     tdia = 0.0
     do k=1,kz
        kr_norm=L1Norm2(kr%submat,k,nsubrow)
        if (kj%cols(k) .eq. irow) tdia = kr_norm 
        tnorm= tnorm + kr_norm
     end do
      if (tnorm > 0.0) tnorm =  tdia / tnorm
      w(irow) = tnorm
      if (tnorm > wmax) wmax = tnorm

   end do
   do irow=1,n
      w(irow) = w(irow)/wmax
      !write(6,*)"weight(",irow,")=",w(irow)
   end do

   nullify(kr)
   nullify(kj)
   return
end subroutine weightsC
!|--------------------------------------------------------------------*/




subroutine lumatvec(mat, x, y)
  !---------------------------------------------------------------------
  !  | This function does the matrix vector product y = A x.
  !  |----------------------------------------------------------------------
  !  | on entry:
  !  | mata  = the matrix (in LUSparMat form)
  !  | x     = a vector
  !  |
  !  | on return
  !  | y     = the product A * x
  !  |--------------------------------------------------------------------*/
  type(ilu_type),pointer::mat
  real(dp),pointer,dimension(:,:)::x,y
  !   local variables    */
  integer :: i,iisub, j,jjsub, nnzcol, nnzrow 
  type(ja_type),pointer::ia,ja 
  type(submat_type),pointer::ma
  type(cs_type),pointer::L,U 

  L => mat%L
  U => mat%U 

  do i = 1,mat%n
     do iisub=1,mat%nsubrow
        y(iisub,i) = 0.0
     end do
  end do
  do i = 1, mat%n
    nnzcol = L%nnzrow(i) 
    ia => L%pj(i) 
    ma => L%pa(i) 
    do iisub=1,mat%nsubrow
       do j = 1,nnzcol
          do jjsub=1,mat%nsubcol
             y(iisub,ia%cols(j)) = y(iisub,ia%cols(j)) + ma%submat(iisub,jjsub,j) * x(iisub,i)
          end do
       end do
    end do
  end do
  do i = 1,mat%n
     nnzrow = U%nnzrow(i) 
     ja => U%pj(i) 
     ma => U%pa(i) 
     do iisub=1,mat%nsubrow
        do j = 1, nnzrow
           do jjsub=1,mat%nsubcol
              y(iisub,i) = y(iisub,i) + ma%submat(iisub,jjsub,j) * x(jjsub,ja%cols(j))
           end do
        end do
     end do
  end do

  nullify(ma)
  nullify(ia)
  nullify(ja)
  nullify(L)
  nullify(L)
  nullify(U)
  return
end subroutine lumatvec
!|---------------------------------------------------------------------*/ 



!real(dp) function vbnorm2(sz, a)!ceb what should dimensions of a be? vector or submat
!  integer ::sz
!  real(dp),dimension(:) ::a
!  !-------------------- return average norm among a(sz) */
!  integer :: tmp
!  tmp = 1 
!  vbnorm2=GNRM2(sz, a, tmp) / sz
!  return  
!end function vbnorm2



!integer condestLU( iluptr lu, real(dp) *y, real(dp) *x, FILE *fp )
subroutine condestLU(lu, y, x)
  type(ilu_type),pointer::lu
  real(dp),pointer,dimension(:,:)::y,x
  integer :: i,iisub
  real(dp) :: norm 

  norm = 0.0 

  do i = 1,lu%n 
     do iisub=1,lu%nsubrow
        y(iisub,i) = 1.0 
     end do
  end do
  call lusolC(y, x, lu) 
  do i = 1,lu%n
     do iisub=1,lu%nsubrow
        norm = max( norm, abs(x(iisub,i)) )
     end do
  end do

  !norm is estimated condition number

  !fprintf( fp, "ILU inf-norm lower bound : %16.2f\n", norm ) 
  !if( norm > 1e30 )  return -1
  return! 0 
end subroutine condestLU
!|---------------------------------------------------------------------*/ 



!integer condestArms(arms armspre, real(dp) *y, FILE *fp )
subroutine condestArms(armspre, y)
  type(arms_data_type),pointer::armspre
  real(dp),pointer,dimension(:,:)::y
  !-------------------- simple estimate of cond. number of precon */
  integer ::i,iisub 
  real(dp) ::norm
  norm=0.0
  do i = 1,armspre%n 
     do iisub=1,armspre%nsubrow 
        y(iisub,i) = 1.0 
     end do
  end do
  call armsol2(y, armspre)   
  do i = 1,armspre%n
     do iisub=1,armspre%nsubrow
        norm = max(norm,abs(y(iisub,i))) 
     end do
  end do

  !norm is estimated condition number

  !fprintf( fp, "ARMS inf-norm lower bound = : %16.2f\n", norm ) 
  !if( norm > 1e30 ) return -1 
  return! 0 
end subroutine condestArms
!|---------------------------------------------------------------------*/ 



subroutine matvecCSR(mat, x, y)
  real(dp),pointer,dimension(:,:)::x,y
  type(SMat_type),pointer::mat
  !-------------------- matvec for csr format using the SMatptr struct*/
  call matvec(mat%CSR, x, y )   
  return
end subroutine matvecCSR
!|---------------------------------------------------------------------*/ 


!subroutine matvecVBR(mat, x, y)
!  real(dp),pointer,dimension(:,:)::x,y
!  type(SMat_type),pointer::mat
!  !-------------------- matvec for vbr format using the SMat struct*/
!  call vbmatvec(mat%VBCSR, x, y )   
!  return
!end subroutine matvecVBR


subroutine matvecLDU(mat, x, y)
  real(dp),pointer,dimension(:,:)::x,y
  type(SMat_type),pointer::mat
  !-------------------- matvec for ldu format using the Smatrix struct*/
  call lumatvec(mat%LDU, x, y )   
  return
end subroutine matvecLDU
!|---------------------------------------------------------------------*/ 

!-------------------- preconditioning operations */

subroutine  preconILU(x, y, mat)
  real(dp),pointer,dimension(:,:)::x,y
  type(SPre_type),pointer::mat
  !-------------------- precon for csr format using the SPre struct*/
  call lusolC(x, y, mat%ILU)   
  return
end subroutine preconILU
!|---------------------------------------------------------------------*/ 


!subroutine  preconVBR(x, y, mat)
!  real(dp),pointer,dimension(:,:)::x,y
!  type(SPre_type),pointer::mat
!  !-------------------- precon for ldu format using the SPre struct*/
!  call vblusolC(x, y, mat%VBILU)   
!  return
!end subroutine preconVBR


subroutine preconLDU(x, y, mat)
  real(dp),pointer,dimension(:,:)::x,y
  type(SPre_type),pointer::mat!{
  !-------------------- precon for vbr format using the SPre struct*/
  call lumsolC(x, y, mat%ILU)   
  return
end subroutine preconLDU
!|---------------------------------------------------------------------*/ 


subroutine preconARMS(x,y,mat)
  real(dp),pointer,dimension(:,:)::x,y
  type(SPre_type),pointer::mat
  !-------------------- precon for ldu format using the SPre struct*/
  integer i,iisub  
  !memcpy(y, x, mat%ARMS%n*sizeof(real(dp)))   
  do i=1,mat%ARMS%n
     do iisub=1,mat%ARMS%nsubrow
        y(iisub,i)=x(iisub,i)
     end do
  end do
  call armsol2(y, mat%ARMS)   
  return
end subroutine preconARMS
!|---------------------------------------------------------------------*/ 


subroutine  Lsolp(start, mata, b, y)
  !---------------------------------------------------------------------
  !  |
  !  | This routine does the forward solve L y = b. where L is upper or
  !  | bottom part of local low triangular matrix
  !  |
  !  | Can be done in place.
  !  |
  !  | Zhongze Li, Aug. 17th, 2001
  !  |
  !  |----------------------------------------------------------------------
  !  | on entry:
  !  | start = the index of the first component
  !  | n     = one ore than the index owned by the processor 
  !  | b     = a vector
  !  | mata  = the matrix (in SparRow form)
  !  |
  !  | on return
  !  | y     = the product L^{-1} * b
  !  |
  !  |--------------------------------------------------------------------*/
  !   local variables    */
  type(cs_type),pointer::mata
  real(dp),pointer,dimension(:,:)::b,y
  integer ::start, i, k,iisub,kksub
  type(submat_type),pointer::kr 
  type(ja_type),pointer ::ki
  
  do i=1,mata%n
     if ( mata%nnzrow(i) .gt. 0 ) then
        kr => mata%pa(i) 
        ki => mata%pj(i) 
        do iisub=1,mata%nsubrow
           y(iisub,i) = b(iisub,i) 
           do k=1,mata%nnzrow(i)
              if(ki%cols(k) < start) then
                 do kksub=1,mata%nsubcol           
                    y(iisub,i) = y(iisub,i) - kr%submat(iisub,kksub,k)*y(kksub,ki%cols(k)) 
                 end do
              end if
           end do
        end do
     end if
  end do
  nullify(kr)
  nullify(ki)
  return
end subroutine Lsolp
!|---------------------------------------------------------------------*/ 


!ceb may have to modify the following routine to properly handle the block diagonal
!void Usolp(integer start, csptr mata, real(dp) *y, real(dp) *x)
subroutine Usolp(start, mata, y, x)
  !{
  !  !---------------------------------------------------------------------
  !    |
  !    | This routine does the backward solve U x = y, where U is upper or
  !    | bottom part of local upper triangular matrix
  !    |
  !    | Can be done in place.
  !    |
  !    | Zhongze Li, Aug. 17th, 2001
  !    |
  !    |----------------------------------------------------------------------
  !    | on entry:
  !    | start = the index of the first component
  !    | n     = one ore than the index owned by the processor 
  !    | y     = a vector
  !    | mata  = the matrix (in SparRow form)
  !    |
  !    | on return
  !    | x     = the product U^{-1} * y
  !    |
  !    |---------------------------------------------------------------------*/
  integer :: start
  type(cs_type),pointer:: mata
  real(dp),pointer,dimension(:,:)::x,y
  !     local variables    */
  integer :: i, k,iisub,kksub 
  type(ja_type),pointer::ki
  type(submat_type),pointer::kr
  real(dp)::t(mata%nsubrow)

  do i=start,1,-1
     kr => mata%pa(i)
     ki => mata%pj(i)
     do iisub=1,mata%nsubrow
        x(iisub,i) = y(iisub,i)
        do k=1,mata%nnzrow(i)
           do kksub=1,mata%nsubcol
              !t(iisub) = t(iisub)- kr%submat(iisub,kksub,k) * x(kksub,ki%cols(k))
              x(iisub,i) = x(iisub,i) - kr%submat(iisub,jjsub,k) * x(kksub,ki%cols(k)) !ceb
           end do
        end do
        !x(iisub,i) = t(iisub)*kr%submat(iisub,iisub,1)
     end do
     !end do
     
     !ceb added
     !//x[i]=D[i]*x[i]
     do iisub=1,mata%nsubrow
        t(iisub)=0.0
        do kksub=1,mata%nsubcol		
           t(iisub) = t(iisub) + kr%submat(iisub,kksub,1)*x(kksub,i)
        end do
     end do
     do iisub=1,mata%nsubrow
        x(iisub,i)=t(iisub)
     end do
     !ceb
  end do

  nullify(kr)
  nullify(ki)
  return
end subroutine Usolp
!----------------------------------------------------------------------*/




!void dscale(int n, double *dd, double *x, double * y)
subroutine dscale(n, nsub, dd, x, y)
  !/* Computes  y == DD * x                               */
  !/* scales the vector x by the diagonal dd - output in y */
  integer :: n,nsub
  real(dp),pointer,dimension(:,:)::dd
  real(dp),pointer,dimension(:,:)::x,y
  integer :: k,kksub
!write(6,*)"dscale checkpt 1"
  do k=1,n
     do kksub=1,nsub
!write(6,*)"dscale: y(",kksub,k,")=dd:",dd(kksub,k)," * x:",x(kksub,k),"=",dd(kksub,k)*x(kksub,k)
        y(kksub,k) = dd(kksub,k)*x(kksub,k)
     end do
  end do
  return
end subroutine dscale
!----------------------------------------------------------------------*/


!subroutine alloc_nnz(array,size)
!  integer,dimension(:),pointer::array
!  integer :: size,istat
!  allocate(array(size),STAT=istat)
!end subroutine alloc_nnz

!subroutine alloc_ja(array,size)
!  type(ja_type),dimension(:),pointer::array
!  integer :: size,istat
!  allocate(array(size),STAT=istat)
!end subroutine alloc_ja

!subroutine alloc_ma(array,size)
!  type(submat_type),dimension(:),pointer::array
!  integer :: size,istat
!  allocate(array(size),STAT=istat)
!end subroutine alloc_ma
!=======================================================================================================

!subroutine alloc_cols(row,size)
!  type(ja_type),pointer :: row
!  integer :: size,istat
!  allocate(row%cols(size),STAT=istat)
!  return
!end subroutine alloc_cols
!=======================================================================================================

!subroutine alloc_submat(row,size,rsub,csub)
!  type(submat_type),pointer :: row
!  integer :: size,rsub,csub,istat
!  allocate(row%submat(rsub,csub,size),STAT=istat)
!end subroutine alloc_submat
!=======================================================================================================


subroutine print_dataP4(amat) 
  !/*----------------------------------------------------------------------
  !| initialize PerMat4 struct given the F, E, blocks.  
  !|----------------------------------------------------------------------
  type(p4_type),intent(in)::amat
  !write(6,*)
  write(6,*)"P4 data"
  write(6,*) "mylev=",amat%mylev
  write(6,*) "n=",amat%n
  !write(6,*) "nsubrow=",amat%nsubrow
  !write(6,*) "nsubcol=",amat%nsubcol
  write(6,*) "nB=",amat%nB
  !write(6,*) "symperm=",amat%symperm
  write(6,*)

  return
end subroutine print_dataP4
!/*---------------------------------------------------------------------
!|     end of setupP4 
!|--------------------------------------------------------------------*/


end module matsol_lib
!=======================================================================================================
!=======================================================================================================
