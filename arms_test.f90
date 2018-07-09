
!===================================================================
subroutine ARMS(A,x,b,dim,neqn,nnz_mat,ia,ja,iau,droptol,fill,schur_start,nnz_pc,ia_pc,ja_pc,iau_pc,lev)
!{
  integer, intent(in) ::fill,dim,neqn,nnz_mat
  integer, intent(inout) ::schur_start,lev
  real(dp), intent(in):: droptol !dropping threshold
  integer, intent(in), dimension(dim+1) :: ia
  integer, intent(in), dimension(nnz_mat) :: ja
  integer, intent(in), dimension(dim) :: iau
  real(dp), intent(in), dimension(neqn,neqn,nnz_mat) :: A
  real(dp), intent(inout), dimension(neqn,dim) :: x
  real(dp), intent(in), dimension(neqn,dim) :: b
  integer, intent(inout) :: nnz_pc
  integer, intent(inout), allocatable, dimension(:) :: ia_pc
  integer, intent(inout), allocatable, dimension(:) :: ja_pc
  integer, intent(inout), allocatable, dimension(:) :: iau_pc
  real(dp), intent(inout), allocatable, dimension(:,:,:):: ILU
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

  integer :: iord(nnodes)
  integer :: riord(nnodes)

  integer,  allocatable,dimension(:)     :: nnz_schur
  integer,  allocatable,dimension(:)     :: jw
  real(dp), allocatable,dimension(:,:,:) :: w
  integer :: levmax
  levmax=2

  !increment recursion level
  lev=lev+1
  !1. reorder current level A matrix and separate into coarse [B] and fine [C] sets
  ! may need to alter the indepSetReorder function to allow starting and ending indices
  ! different than 1..nnodes
  call indepSetReorder(A,ia,ja,nnodes,nsub,nnz,bsize,iord,riord,ind_nnodes,indtol,nbnd) 
  ![A_lev]=|B F| 
  !        |E C| 
  !reordered A is kept throughout as it is defined and allocated external to this function

  !2. factor coarse [B]=[L][U] matrix via ILU
  ! |B F|=|[L]       [0]||[U] [L^-1][F]| 
  ! |E C| |[E][U^-1] [I]||[0] [A_lev+1]| 
  ! This ilut takes matrix A_lev and factors the coarse set [B] into LU and also 
  ! creates the schur complement matrix A_lev+1
  ! also creates W=[L^-1][F]
  !              G=[E][U^-1]
  !ceb need to add reordering map to blkilut or can we just send in ILU as A reordered
!ceb may want to modify this to output L and U instead of combined ILU
! or maybe even L, U, W, G, S
  call BLKILUT(dim,neqn,nnz_mat,ia,ja,iau,A,ILU,droptol,fill,ind_nnodes+1,nnz_pc,ia_pc,ja_pc,iau_pc)
  !call BLKILUT_SCHUR(dim,neqn,nnz_mat,ia,ja,iau,A,iord,riord,L,U,droptol,fill,ind_nnodes+1,nnz_pc,ia_pc,ja_pc,iau_pc)



  !3. solve [L]b' = b,  b=rhs for independent set  
  !i.e. b'=[L^-1]b
  ! alter this to use L from BLKILUT_SCHUR
  !//compute bp[i]= b[i]-L[p,j]*bp[j]
  !ceb indices may need to be shifted as we recurse to higher levels
  do i=1,ind_nnodes
     p=iord(i)
     do ii=1,neqn
        bp(ii,p)=b(ii,p)!//x=rhs
        do j=1,L_rownnz(i)!//Lower Matrix
           jnode=iord(L_ja(i)%cols(j))
              do jj=1,neqn
                 bp(ii,p)=bp(ii,p)-L(ii)%submat(ii,jj,j)*bp(jj,jnode)
              end do
        end do
     end do
  end do



  !4. compute h' = h-[E][U^-1]b' = h-[G]b', h=rhs for schur complement set 
! modify to compute as h-[E][U^-1]b'

  !compute [U^-1]bp
  !U[i,j]*x[j]=y[i](off diagonals only)!ceb fix this in representation of U,L
  do i=1,ind_nnodes
     p=iord(i)
     do ii=1,neqn
        bp(ii,p)=b(ii,p)!//rhs
        do j=1,U_rownnz(i)!//Upper Matrix
           jnode=iord(U_ja(i)%cols(j))
              do jj=1,neqn
                 bp(ii,p)=bp(ii,p)-U(ii)%submat(ii,jj,j)*bp(jj,jnode)
              end do
        end do
     end do
  end do

  do i=ind_nnodes+1,nnodes
     p=i
     p=iord(i)
     m1=ia(p)
     m2=iau(p)-1
     do ii=1,neqn
        hp(ii,p)=h(ii,p)
        do j=m1,m2!//
           jnode=ja(j)
           jnode=iord(ja(j))
           if(jnode .le. ind_nnodes) then
              !compute h-[E][U^-1]b'
              do jj=1,neqn
                 hp(ii,p)= hp(ii,p)-A(ii,jj,j)*Ubp(jj,jnode)
              end do
           end if
        end do
     end do
  end do


  ! Now we need the schur complement to pass to the next level
  ! A_lev+1 = S = [C]-[E][B^-1][F]

  ! z is solution for dependent part of solution vector x
  ! ideally we could send a separate A to either of these functions 
  ! without having to include previous ordering info
  !i.e. if A_lev+1 were a separate matrix with its own indexing although
  ! it is a submatrix of A_lev+0
  if (lev==levmax) then
  !5.   solve [A_lev+1]z=hp via GMRES+preconditioner
  !  call GMRES(A_lev+1,z,hp)!(need version that can work with the schur complement (split) matrices)
     ! or we can figure out how to send it the A_+1 matrix and vector as a normal matrix
  else
  !6.  call ARMS(A_lev+1,z,hp) (recursive)
  end if

  !ceb No longer need [A_lev+1] or [G] or [L]
  !    we need only [U],[W] or [U],[L],[F]
  !    so free unecessary memory

  !7. compute g = b'-[L^-1][F]z = b'-[W]z

  do i=1,ind_nnodes
     p=iord(i) !use p in place of i
     !compute [F]z
     m1=iau(p)+1
     m2=ia(p+1)-1
     do j=m1,m2
        jcol=ja(j)
        jnode=iord(jcol)
        !only dependent set columns make up submatrix F
        if(jcol.gt.ind_nnodes) then
           do ii=1,neqn
              Fz(ii,p)=0
              do jj=1,neqn
                 Fz(ii,p)=Fz(ii,p)+A(ii,jj,j)*z(jj,jnode)
              end do
           end do
        end if
     end do
  end do

  !we could reduce storage of Fz by recomputing Fz for each jnode
  !in the following loop
  do i=1,ind_nnodes
     p=iord(i) !use p in place of i
     ! compute g = b'-[L^-1][F]z
     do ii=1,neqn
        g(ii,p)=bp(ii,p)
        do j=1,L_rownnz(i) !//
           jnode=iord(L_ja(i)%cols(j))
           do jj=1,neqn
              g(ii,p)=g(ii,p)-L(i)%submat(ii,jj,j)*Fz(jj,jnode)
           end do
        end do
     end do
  end do



  ! y is solution for independent part of solution vector x
  !8. backsub y=[U^-1]g
  !//compute g[p]= y[p]-U[p,j]*g[j]
  do i=nnodes_indep,1,-1
     p=iord(i) !use p in place of i  diag=iau(p)
     do ii=1,neqn
        do j=1,U_rownnz(i)!//Upper Matrix
           jnode=iord(U_ja(i)%cols(j)) 
           do jj=1,neqn
              x(ii,p) = x(ii,p)-U(i)%submat(ii,jj,j)*g(jj,jnode)
           end do
        end do
     end do
     !ceb [D] is currently a combination of L and U
     !//x[p]=D[p]*x[p]
     do ii=1,neqn
        t(ii)=0.0
        do jj=1,neqn		
           t(ii) = t(ii)+U(i)%submat(ii,jj,1)*x(jj,p)
        end do
     end do
     do kk=1,neqn
        x(kk,p)=t(kk)
     end do
  end do

!}
end subroutine ARMS
!===================================================================

