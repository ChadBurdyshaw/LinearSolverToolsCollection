module blas
  use kinddefs, only : dp
  implicit none


interface tmemcpy
   function tmemcpy_r(a,b,size) result(retval)
     use kinddefs, only : dp
     real(dp)::a(size), b(size)
     integer::size,retval
   end function tmemcpy_r
   function tmemcpy_i(a,b,size) result(retval)
     integer::a(size), b(size)
     integer::size,retval
   end function tmemcpy_i
end interface

interface tblank
   function tblank_r(a,size) result(retval)
     use kinddefs, only : dp
     real(dp)::a(size)
     integer::size,retval
   end function tblank_r
   function tblank_i(a,size) result(retval)
     integer::a(size)
     integer::size,retval
   end function tblank_i
end interface


contains

  !=============================================================================
  real(dp) function snrm2 ( n, sx, incx)
      integer :: n,incx,next
      integer :: nn,i,j
      !real(dp),dimension(:)::sx 
      !real(dp)::sx(1) 
      real(dp)::sx(:) 
      real(dp) :: cutlo, cuthi, hitest, sum, xmax, zero, one,snrm

      zero =0.0e0
      one= 1.0e0
!
!     euclidean norm of the n-vector stored in sx() with storage
!     increment incx .
!     if    n .le. 0 return with result = 0.
!     if n .ge. 1 then incx must be .ge. 1
!
!           c.l.lawson, 1978 jan 08
!
!     four phase method     using two built-in constants that are
!     hopefully applicable to all machines.
!         cutlo = maximum of  sqrt(u/eps)  over all known machines.
!         cuthi = minimum of  sqrt(v)      over all known machines.
!     where
!         eps = smallest no. such that eps + 1. .gt. 1.
!         u   = smallest positive no.   (underflow limit)
!         v   = largest  no.            (overflow  limit)
!
!     brief outline of algorithm..
!
!     phase 1    scans zero components.
!     move to phase 2 when a component is nonzero and .le. cutlo
!     move to phase 3 when a component is .gt. cutlo
!     move to phase 4 when a component is .ge. cuthi/m
!     where m = n for x() real and m = 2*n for complex.
!
!     values for cutlo and cuthi..
!     from the environmental parameters listed in the imsl converter
!     document the limiting values are as follows..
!     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
!                   univac and dec at 2**(-103)
!                   thus cutlo = 2**(-51) = 4.44089e-16
!     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
!                   thus cuthi = 2**(63.5) = 1.30438e19
!     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
!                   thus cutlo = 2**(-33.5) = 8.23181d-11
!     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
!     data cutlo, cuthi / 8.232d-11,  1.304d19 /
!     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      cutlo=4.441e-16
      cuthi=1.304e19
!
      if(n .gt. 0) go to 10
         snrm2  = zero
         go to 300 
!
   10 assign 30 to next
      sum = zero
      nn = n * incx
!                                                 begin main loop
      i = 1
   20    go to next,(30, 50, 70, 110)
   30 if( abs(sx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
!
!                        phase 1.  sum is zero
!
   50 if( sx(i) .eq. zero) go to 200
      if( abs(sx(i)) .gt. cutlo) go to 85
!
!                                prepare for phase 2.
      assign 70 to next
      go to 105
!
!                                prepare for phase 4.
!
  100 i = j
      assign 110 to next
      sum = (sum / sx(i)) / sx(i)
  105 xmax = abs(sx(i))
      go to 115
!
!                   phase 2.  sum is small.
!                             scale to avoid destructive underflow.
!
   70 if( abs(sx(i)) .gt. cutlo ) go to 75
!
!                     common code for phases 2 and 4.
!                     in phase 4 sum is large.  scale to avoid overflow.
!
  110 if( abs(sx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / sx(i))**2
         xmax = abs(sx(i))
         go to 200
!
  115 sum = sum + (sx(i)/xmax)**2
      go to 200
!
!
!                  prepare for phase 3.
!
   75 sum = (sum * xmax) * xmax
!
!
!     for real or d.p. set hitest = cuthi/n
!     for complex      set hitest = cuthi/(2*n)
!
   85 hitest = cuthi/float( n )
!
!                   phase 3.  sum is mid-range.  no scaling.
!
      do 95 j =i,nn,incx
      if(abs(sx(j)) .ge. hitest) go to 100
   95    sum = sum + sx(j)**2
      snrm2 = sqrt( sum )
      go to 300
!
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
!
!              end of main loop.
!
!              compute square root and adjust for scaling.
!
      snrm2 = xmax * sqrt(sum)
  300 continue
      return
   
end function snrm2





!=============================================================================
subroutine sscal(n,sa,sx,incx)
  !
  !     scales a vector by a constant.
  !     uses unrolled loops for increment equal to 1.
  !     jack dongarra, linpack, 3/11/78.
  !
  !real(dp) sa,sx(1)
  real(dp) :: sa
  real(dp),dimension(:)::sx
  integer :: i,incx,m,mp1,n,nincx

!write(6,*) "sscal n=",n," sa=",sa," "
  !
  if(n.le.0)return
!  if(incx.eq.1)go to 20
  !
  !        code for increment not equal to 1
  !
  nincx = n*incx
  do  i = 1,nincx,incx
     sx(i) = sa*sx(i)
     !10   continue
  end do
  return
  !
  !        code for increment equal to 1
  !
  !
  !        clean-up loop
  !
!20 m = mod(n,5)
!  if( m .eq. 0 ) go to 40
  do i = 1,m
     sx(i) = sa*sx(i)
     !30      continue
  end do
!  if( n .lt. 5 ) return
!40 mp1 = m + 1
!  do i = mp1,n,5
!     sx(i) = sa*sx(i)
!     sx(i + 1) = sa*sx(i + 1)
!     sx(i + 2) = sa*sx(i + 2)
!     sx(i + 3) = sa*sx(i + 3)
!     sx(i + 4) = sa*sx(i + 4)
!  end do
  !50 continue
  return
end subroutine sscal





!=============================================================================
subroutine scopy(n,sx,incx,sy,incy)
  !
  !     copies a vector, x, to a vector, y.
  !     uses unrolled loops for increments equal to 1.
  !     jack dongarra, linpack, 3/11/78.
  !
  !real(dp) sx(1),sy(1)
  real(dp),intent(in)   :: sx(:)
  real(dp),intent(inout):: sy(:)
  integer :: i,incx,incy,ix,iy,m,mp1,n
  !
  if(n.le.0)return
  if(incx.eq.1.and.incy.eq.1)go to 20
  
  ! code for unequal increments or equal increments
  ! not equal to 1
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do  i = 1,n
     sy(iy) = sx(ix)
     ix = ix + incx
     iy = iy + incy
  end do
  return

  !code for both increments equal to 1
20 continue!ceb
  do  i = 1,n
     sy(i) = sx(i)
  end do

  return
end subroutine scopy





!=============================================================================
subroutine saxpy(n,sa,sx,incx,sy,incy)
  !
  !     constant times a vector plus a vector.
  !     uses unrolled loop for increments equal to one.
  !     jack dongarra, linpack, 3/11/78.
  !
  real(dp),intent(in) ::sa
  real(dp),intent(in)::sx(:)
  real(dp),intent(inout)::sy(:)
  integer ::i,incx,incy,ix,iy,m,mp1,n
  !
  if(n.le.0)return
  if (sa .eq. 0.0) return
  if(incx.eq.1.and.incy.eq.1)go to 20
  !
  !        code for unequal increments or equal increments
  !          not equal to 1
  !
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i = 1,n
     sy(iy) = sy(iy) + sa*sx(ix)
     ix = ix + incx
     iy = iy + incy
  end do
  return
  !
  !        code for both increments equal to 1
  !
  !
  !        clean-up loop
  !
20 m = mod(n,4)
  if( m .eq. 0 ) go to 40
  do i = 1,m
     sy(i) = sy(i) + sa*sx(i)
  end do
  if( n .lt. 4 ) return
40 mp1 = m + 1
  do i = mp1,n,4
     sy(i) = sy(i) + sa*sx(i)
     sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
     sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
     sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
  end do
  return
end subroutine saxpy






!=============================================================================
real(dp) function sdot(n,sx,incx,sy,incy)
  !
  !     forms the dot product of two vectors.
  !     uses unrolled loops for increments equal to one.
  !     jack dongarra, linpack, 3/11/78.
  !
  !real(dp),dimension(:)::sx,sy
  !real(dp)::sx(1),sy(1)
  real(dp),intent(in)::sx(:),sy(:)
  real(dp) ::stemp
  integer :: i,incx,incy,ix,iy,m,mp1,n
  !
  stemp = 0.0e0
  sdot = 0.0e0
  if(n.le.0)return
  if(incx.eq.1.and.incy.eq.1)go to 20
  !
  !        code for unequal increments or equal increments
  !          not equal to 1
  !
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i = 1,n
     stemp = stemp + sx(ix)*sy(iy)
     ix = ix + incx
     iy = iy + incy
  end do
  sdot = stemp
  return
  !
  !        code for both increments equal to 1
  !
  !
  !        clean-up loop
  !
20 m = mod(n,5)
  if( m .eq. 0 ) go to 40
  do i = 1,m
     stemp = stemp + sx(i)*sy(i)
  end do
  if( n .lt. 5 ) go to 60
40 mp1 = m + 1
  do i = mp1,n,5
     stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) + &
          sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
  end do
60 sdot = stemp
  return
end function sdot



!=============================================================================
real(dp) function sdot2(n,sx,incx,sy,incy)
  !
  !     forms the dot product of two vectors.
  !     uses unrolled loops for increments equal to one.
  !     jack dongarra, linpack, 3/11/78.
  !
  !real(dp),dimension(:)::sx,sy
  !real(dp)::sx(1),sy(1)
  real(dp)::sx(:),sy(:,:)
  real(dp) ::stemp
  integer :: i,j,incx,incy,ix,iy,m,mp1,n
  !
  stemp = 0.0e0
  sdot2 = 0.0e0
  if(n.le.0)return
!  if(incx.eq.1.and.incy.eq.1)go to 20
  !
  !        code for unequal increments or equal increments
  !          not equal to 1
  !
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i = 1,n
     do j=1,incy
     stemp = stemp + sx(ix)*sy(i,j)
  end do
     ix = ix + incx
     !iy = iy + incy
  end do
  sdot2 = stemp
  return

  !
  !        code for both increments equal to 1
  !
  !
  !        clean-up loop
  !
!20 m = mod(n,5)
!  if( m .eq. 0 ) go to 40
!  do i = 1,m
!     stemp = stemp + sx(i)*sy(i)
!  end do
!  if( n .lt. 5 ) go to 60
!40 mp1 = m + 1
!  do i = mp1,n,5
!     stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) + &
!          sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
!  end do
!60 sdot = stemp
!  return
end function sdot2








!=============================================================================
subroutine srotg(sa,sb,c,s)
  !
  !     construct givens plane rotation.
  !     jack dongarra, linpack, 3/11/78.
  !
  real(dp) sa,sb,c,s,roe,scale,r,z,one
  !
  one=1.0
  roe = sb
  if( abs(sa) .gt. abs(sb) ) roe = sa
  scale = abs(sa) + abs(sb)
  
  if( scale .ne. 0.0 ) then
     r = scale*sqrt((sa/scale)**2 + (sb/scale)**2)
     r = sign(one,roe)*r
     c = sa/r
     s = sb/r
  else
     c = 1.0
     s = 0.0
     r = 0.0
  end if
20 z = 1.0
  if( abs(sa) .gt. abs(sb) ) z = s
  if( abs(sb) .ge. abs(sa) .and. c .ne. 0.0 ) z = 1.0/c
  sa = r
  sb = z
  return
end subroutine srotg

!      function second()
!      real tarray(2)
!      second = etime(tarray)
!      return
!      end





!================================================
!//computes frobenius norm of square elemental matrix
real(dp) function FNorm(A,neqn)
  integer :: i,j,neqn
  real(dp) :: A(neqn,neqn)
  real(dp) sum
  sum=0
  do i=1,neqn
     do j=1,neqn
        sum = sum + A(j,i)*A(j,i)
     end do
  end do
  FNorm=sqrt(sum/(neqn*neqn))
  return 
end function FNorm

!================================================
!//computes frobenius norm of square elemental matrix
real(dp) function L1Norm(A,neqn)
  integer :: i,j,neqn
  real(dp) :: A(neqn,neqn)
  real(dp) sum
  sum=0
  do i=1,neqn
     do j=1,neqn
        sum = sum + abs(A(j,i))
     end do
  end do
  L1Norm= sum/(neqn*neqn)
  return 
end function L1Norm

!================================================
!//computes frobenius norm of square elemental matrix
real(dp) function FNorm2(A,i,neqn,nnz)
  integer :: i,ii,jj,neqn,nnz
  real(dp) :: A(neqn,neqn,nnz)
  real(dp) sum
  sum=0
  do ii=1,neqn
     do jj=1,neqn
        sum = sum + A(jj,ii,i)*A(jj,ii,i)
     end do
  end do
  FNorm2=sqrt(sum/(neqn*neqn))
  return 
end function FNorm2

!================================================
!//computes frobenius norm of square elemental matrix
!real(dp) function L1Norm2(A,i,neqn,nnz)
real(dp) function L1Norm2(A,i,neqn)
  integer :: i,ii,jj,neqn
  real(dp),dimension(:,:,:) :: A
  real(dp) :: sum
  sum=0
  do ii=1,neqn
     do jj=1,neqn
        sum = sum + abs(A(jj,ii,i))
     end do
  end do
  L1Norm2= sum/(neqn*neqn)
  return 
end function L1Norm2



!/* this function mutiplies a block matrix A of dimension dim in compressed row storage by a vector x of length BLOCKSIZE*dim and stores the solution in the vector Ax with dimension [BLOCKSIZE*dim] */
!======================================================
subroutine matvect(ia,ja, A, x, Ax, dim, neqn, nnz,blank)
  integer :: dim,neqn,nnz,blank
  integer ,intent(in), dimension(dim+1) :: ia
  integer ,intent(in), dimension(nnz) :: ja
  real(dp),intent(in), dimension(neqn,neqn,nnz) :: A
  real(dp),intent(in), dimension(neqn,dim) :: x
  real(dp),intent(inout), dimension(neqn,dim)  :: Ax
  integer :: j, k, jstart, jend, ii, jj, icol

  if (blank.eq.1) call blank_2D_mat(Ax,neqn,dim)

  do k=1,dim 
     jstart = ia(k)
     jend = ia(k+1)
     do j=jstart,jend-1	 
        icol=ja(j)
        do ii=1,neqn
           do jj=1,neqn
              Ax(ii,k) = Ax(ii,k) + A(ii,jj,j)*x(jj,icol)
           end do
        end do
     end do
  end do
  return
end subroutine matvect





!================================================
subroutine blank_1D_mat(A,n1)
  integer :: n1
  integer ii,jj
  real(dp) :: A(n1)
  do ii=1,n1
        A(ii)=0.0
  end do
  return
end subroutine blank_1D_mat

!================================================
subroutine blank_2D_mat(A,n1,n2)
  integer :: n1,n2
  integer ii,jj
  real(dp) :: A(n1,n2)
  do ii=1,n1
     do jj=1,n2
        A(ii,jj)=0.0
     end do
  end do
  return
end subroutine blank_2D_mat

!================================================
subroutine blank_3D_mat(A,n1,n2,n3)
  integer :: n1,n2,n3
  integer ii,jj,kk
  real(dp) :: A(n1,n2,n3)
  do kk=1,n3
     do jj=1,n2
        do ii=1,n1
           A(ii,jj,kk)=0.0
        end do
     end do
  end do
  return
end subroutine blank_3D_mat
!================================================
subroutine blank_3D_submat(A,i,n1,n2)
  integer :: i,n1,n2
  integer ii,jj
  real(dp) ,dimension(:,:,:):: A!(n1,n2,n3)
  do jj=1,n2
     do ii=1,n1
        A(ii,jj,i)=0.0
     end do
  end do
  return
end subroutine blank_3D_submat


!================================================
!/* This routine solves X=A^-1 using LU factorization, where X and B are matrices */
subroutine InvertSubMatrx(A,X,neqn)
  integer ii,jj,kk,neqn
  real(dp) :: LUtmp(neqn,neqn)
  real(dp) :: Z(neqn,neqn)
  real(dp) :: Xtemp(neqn,neqn)
  real(dp) :: B(neqn,neqn)
  real(dp),intent(in),dimension(neqn,neqn) :: A
  real(dp) :: X(neqn,neqn)
  real(dp) :: Xtmp(neqn,neqn)

! set B as identity
  do ii=1,neqn
     do jj=1,neqn
        B(ii,jj)=0.0
        if (ii.eq.jj) B(ii,jj) = 1.0
     end do
  end do
  
  call LU_sub(LUtmp,A,neqn)
  !// solving L factors
  do ii=1,neqn
     do jj=1,neqn	  
        Z(ii,jj) = B(ii,jj)
        do kk=1,ii-1	      
           Z(ii,jj) = Z(ii,jj) - LUtmp(ii,kk)*Z(kk,jj)
        end do
     end do
  end do
  !// solving U factors
  do ii=neqn,1,-1
     if(abs(LUtmp(ii,ii))<1.0e-10) then
        return
     end if
     do jj=1,neqn	
        do kk=neqn,ii+1,-1
           Z(ii,jj) = Z(ii,jj) - LUtmp(ii,kk)*Xtmp(kk,jj)
        end do
        Xtmp(ii,jj) = Z(ii,jj)/LUtmp(ii,ii)
     end do
  end do
  !// writing solution onto solution matrix
  call copy_submat(Xtmp,X,neqn)
  
  return
end subroutine InvertSubMatrx


!================================================
!//converts A diagonal submatrix to LU form and scales rows by diagonal entry
subroutine LU_sub(LU,A,neqn)
  integer :: neqn
  real(dp) :: LU(neqn,neqn)
  real(dp) :: A(neqn,neqn)
  call copy_submat(A,LU,neqn)
  call LUgeneral_rm(LU,neqn)
  return
end subroutine LU_sub




!================================================
!// LU decomposition of a matrix in array format, row major order.
subroutine LUgeneral_rm(D,n)
  integer i,j,k,n
  real(dp) sum
  real(dp) :: D(n,n)

  ! Segment 1!
  ! Note that the lower triangular matrix has D(.,1) in the first column so nothing gets done
  ! Complete the first row of the upper triangular component  
  do j = 2,n
     D(1,j) = D(1,j)/D(1,1)
  end do
  ! Now loop over each "segment" starting with row 2
  do i=2,n
     ! First do the column corresponding to the segment
     do k = i,n
        sum = 0.
        do j = 1,i-1
           sum = sum + D(k,j)*D(j,i)
        end do
        D(k,i) = D(k,i) - sum
     end do    
     ! Now do the upper triangular portion (the row)
     do k = i+1,n
        sum = 0.
        do j = 1,i-1
           sum = sum + D(i,j)*D(j,k)
        end do
        D(i,k) = (D(i,k) - sum)/D(i,i)
     end do
  end do

  return
end subroutine LUgeneral_rm




!================================================
subroutine get_submat(A,Asub,i,nnz,nsub)
  integer :: i,nnz,nsub
  integer ii,jj
  real(dp) :: A(nsub,nsub,nnz),Asub(nsub,nsub)
  do ii=1,nsub
     do jj=1,nsub
        Asub(ii,jj)=A(ii,jj,i)
     end do
  end do
  return
end subroutine get_submat

!================================================
subroutine set_submat(A,Asub,i,nnz,nsub)
  integer :: i,nnz,nsub
  integer ii,jj
  real(dp) :: A(nsub,nsub,nnz),Asub(nsub,nsub)
  do ii=1,nsub
     do jj=1,nsub
        A(ii,jj,i)=Asub(ii,jj)
     end do
  end do
  return
end subroutine set_submat

!================================================
subroutine copy_submat(Aold,Anew,nsub)
  integer :: nsub
  integer ii,jj
  real(dp) :: Aold(nsub,nsub),Anew(nsub,nsub)
  do ii=1,nsub
     do jj=1,nsub
        Anew(ii,jj)=Aold(ii,jj)
     end do
  end do
  return
end subroutine copy_submat

!================================================
subroutine mult_submat_const(A,B,C,neqn)
  integer :: neqn
  integer ii,jj
  real(dp) :: A(neqn,neqn),B(neqn,neqn),C
  do ii=1,neqn
     do jj=1,neqn
        A(ii,jj)=C*B(ii,jj)
     end do
  end do
  return
end subroutine mult_submat_const

!================================================
!subroutine add_submat_const(A,B,C,neqn)
!  integer :: neqn
!  integer ii,jj
!  real(dp) :: A(neqn,neqn),B(neqn,neqn),C
!  do ii=1,neqn
!     do jj=1,neqn
!        A(ii,jj)=C+B(ii,jj)
!     end do
!  end do
!  return
!end subroutine add_submat_const


!================================================
!subroutine mult_submat_submat(A,B,neqn)
!  integer :: neqn
!  integer ii,jj,kk
!  real(dp) :: M(neqn,neqn),A(neqn,neqn),B(neqn,neqn)
!  
!  do ii=1,neqn
!     do jj=1,neqn
!        do kk=1,neqn
!           M(ii,jj)=M(ii,jj)+A(ii,kk)*B(kk,jj)
!        end do
!     end do
!  end do
!  !copy M into A
!  do ii=1,neqn
!     do jj=1,neqn
!        A(ii,jj)=M(ii,jj)
!     end do
!  end do
!  return
!end subroutine mult_submat_submat

!================================================
!subroutine mult_submat_submat(M,A,B,neqn)
!  integer :: neqn
!  integer ii,jj,kk
!  real(dp) :: M(neqn,neqn),A(neqn,neqn),B(neqn,neqn)
!  !blank M
!  do ii=1,neqn
!     do jj=1,neqn
!        M(ii,jj)=0
!     end do
!  end do
!  
!  do ii=1,neqn
!     do jj=1,neqn
!        do kk=1,neqn
!           M(ii,jj)=M(ii,jj)+A(ii,kk)*B(kk,jj)
!        end do
!     end do
!  end do
!
!  return
!end subroutine mult_submat_submat

!================================================
!subroutine add_submat_submat(A,B,neqn)
!  integer :: neqn
!  integer ii,jj
!  real(dp) :: A(neqn,neqn),B(neqn,neqn)
!  
!  do ii=1,neqn
!     do jj=1,neqn
!        A(ii,jj)=A(ii,jj)+B(ii,jj)
!     end do
!  end do
!
!  return
!end subroutine add_submat_submat




!===============================================================
! here a is an array of neqn*neqn submatrices
subroutine qsplitC(a,ind,n,ncut,neqn,ofstind)
!----------------------------------------------------------------------
!|     does a quick-sort split of a real array.
!|     on input a[1 : n] is a real array
!|     on output is permuted such that its elements satisfy:
!|
!|     abs(a[i]) >= abs(a[ncut]) for i < ncut and
!|     abs(a[i]) <= abs(a[ncut]) for i > ncut
!|
!|     ind[1 : n] is an integer array permuted in the same way as a.
!|---------------------------------------------------------------------
  integer,intent(in) :: neqn,ncut,n
  integer,intent(in) :: ofstind
  real(dp),intent(inout) ::a(neqn,neqn,n)
  integer,intent(inout) :: ind(n)
!internal vars
  integer :: iisub,jjsub,first
  real(dp) ::tmat(neqn,neqn)
  real(dp) :: abskey
  real(dp) :: tmp
  integer :: j, itmp, mid, last

  first=1
  last = n

!write(6,*) "qsplitC: checkpt0"
!write(6,*) "qsplitC:first=",first," ncut=",ncut," last=",last

  if ( ncut<first .or. ncut>last) then !ncut not in domain
     !write(6,*) "qsplitC: element ",ncut," not in domain"
     return
  end if
  
  !/* outer loop -- while mid != ncut */
  do !label1:
     mid = first
     if(mid.eq.ncut) exit !exit outer do loop
!write(6,*) "qsplitC: checkpt1"

     abskey = L1Norm2(a,mid+ofstind,neqn)
!write(6,*) "qsplitC: checkpt2"

     do j=first+1,last 
        if (L1Norm2(a,j+ofstind,neqn) > abskey) then
!write(6,*) "qsplitC: 0 swapping col ",mid," with col",j!ceb
           mid=mid+1
           itmp = ind(mid+ofstind)
           ind(mid+ofstind) = ind(j+ofstind)
           ind(j+ofstind) = itmp
           do iisub =1,neqn
              do jjsub =1,neqn
                 tmp = a(iisub,jjsub,mid+ofstind)
                 a(iisub,jjsub,mid+ofstind) = a(iisub,jjsub,j+ofstind)
                 a(iisub,jjsub,j+ofstind) = tmp
              end do
           end do
        end if
     end do

!write(6,*) "qsplitC: checkpt3"

     !ceb why here
     if (mid.ne.first) then!ceb
!write(6,*) "qsplitC: 1 swapping col ",mid," with col",first
        !/* interchange */
        do iisub=1,neqn
           do jjsub =1,neqn
              tmp=a(iisub,jjsub,mid+ofstind) !tmp = a[mid]
              a(iisub,jjsub,mid+ofstind)=a(iisub,jjsub,first+ofstind) !a[mid] = a[first];
              a(iisub,jjsub,first+ofstind)=tmp !a[first]  = tmp;
           end do
        end do
        
        itmp = ind(mid+ofstind)
        ind(mid+ofstind) = ind(first+ofstind)
        ind(first+ofstind) = itmp
     end if

!write(6,*) "qsplitC: checkpt4"

     !/* test for while loop */
     if (mid == ncut) return   
     if (mid > ncut) then 
        last = mid-1
     else
        first = mid+1
     end if
  end do ! goto label1
end subroutine qsplitC
!==================================================



!===============================================================
! here a is an array of neqn*neqn submatrices
subroutine qsplit(a,ind,n,ncut,neqn)
!----------------------------------------------------------------------
!|     does a quick-sort split of a real array.
!|     on input a[1 : n] is a real array
!|     on output is permuted such that its elements satisfy:
!|
!|     abs(a[i]) >= abs(a[ncut]) for i < ncut and
!|     abs(a[i]) <= abs(a[ncut]) for i > ncut
!|
!|     ind[1 : n] is an integer array permuted in the same way as a.
!|---------------------------------------------------------------------
  integer,intent(in) :: neqn,ncut,n
  !integer,intent(inout) :: first
  real(dp),intent(inout) ::a(neqn,neqn,n)
  integer,intent(inout) :: ind(n)
!internal vars
  integer :: first
  integer :: iisub,jjsub
  real(dp) ::tmat(neqn,neqn)
  real(dp) :: abskey
  real(dp) :: tmp
  integer :: j, itmp, mid, last

!write(6,*) "qsplit: checkpt0"

  !if(first <= 0) 
  first=1
  last = n
!write(6,*) "qsplit: 0: first=",first," ncut=",ncut," last=",last

  if ( ncut<first .or. ncut>last) then !ncut not in domain
     write(6,*) "qsplitC: element ",ncut," not in domain"
     return
  end if
  
  !/* outer loop -- while mid != ncut */
  do !label1:
     mid = first
     if(mid.eq.ncut) exit !exit outer do loop
!write(6,*)"qsplit: calling L1Norm2"
     abskey = L1Norm2(a,mid,neqn)

     do j=first+1,last 
        if (L1Norm2(a,j,neqn) > abskey) then
           !write(6,*) "qsplit: 0 swapping col ",mid," with col",j!ceb
           mid=mid+1
           itmp = ind(mid)
           ind(mid) = ind(j)
           ind(j) = itmp
           do iisub =1,neqn
              do jjsub =1,neqn
                 tmp = a(iisub,jjsub,mid)
                 a(iisub,jjsub,mid) = a(iisub,jjsub,j)
                 a(iisub,jjsub,j) = tmp
              end do
           end do
        end if
     end do

     !ceb why here
     if (mid.ne.first) then!ceb
        !write(6,*) "qsplit: 1 swapping col ",mid," with col",first
        !/* interchange */
        do iisub=1,neqn
           do jjsub =1,neqn
              tmp=a(iisub,jjsub,mid) !tmp = a[mid]
              a(iisub,jjsub,mid)=a(iisub,jjsub,first) !a[mid] = a[first];
              a(iisub,jjsub,first)=tmp !a[first]  = tmp;
           end do
        end do
        
        itmp = ind(mid)
        ind(mid) = ind(first)
        ind(first) = itmp
     end if

     !/* test for while loop */
     if (mid == ncut) return   


     if (mid > ncut) then 
        last = mid-1
     else
        first = mid+1
     end if
  end do ! goto label1
  
!write(6,*) "qsplit: 2: first=",first," last=",last



end subroutine qsplit
!==================================================

   


recursive subroutine qqsort(ja,ma,left,right)
!/*----------------------------------------------------------------------
!|
!| qqsort: sort ja[left]...ja[right] into increasing order
!| from Kernighan & Ritchie
!|
!| ma holds the real values
!|
!|---------------------------------------------------------------------*/
  real(dp),intent(inout),dimension(:)::ma
  integer,intent(inout),dimension(:)::ja
  integer,intent(inout) :: left,right
  !internal vars
  integer :: i, last,itmp 
  if (left >= right)  return
  itmp=(left+right)/2
  call swapj(ja, left, itmp)
  call swapm(ma, left, itmp)
  last = left
  do i=left+1,right
    if (ja(i) < ja(left)) then
       last=last+1
       call swapj(ja, last, i)
       call swapm(ma, last, i)
    end if
  end do
  call swapj(ja, left, last)
  call swapm(ma, left, last)
  itmp=last-1
  call qqsort(ja, ma, left, itmp)
  itmp=last+1
  call qqsort(ja, ma, itmp, right)
  return
end subroutine qqsort


!void qsort2C(int *ja, double *ma, int left, int right, int abval){
recursive subroutine qsort_intarray(ja,left,right)
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
  integer,intent(inout) :: left,right
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
     end if
  end do
  call swapj(ja, left, last)
  itmp=last-1
  call qsort_intarray(ja,left,itmp)
  itmp=last+1
  call qsort_intarray(ja,itmp,right)
  return
end subroutine qsort_intarray






!void qsortR2I(double *wa, int *cor1, int *cor2, int left, int right){
recursive subroutine qsortR2I(wa, cor1, cor2, left, right)
!/*----------------------------------------------------------------------
!|
!| qqsort: sort wa[left]...wa[right] into decreasing order
!| from Kernighan & Ritchie
!|
!|---------------------------------------------------------------------*/
  real(dp),intent(inout),dimension(:)::wa
  integer,intent(inout),dimension(:)::cor1,cor2
  integer,intent(in):: left,right
  !internal vars
  integer:: i, last,itmp
  if (left >= right) return
  call swapm(wa, left, (left+right)/2)!swap submat
  call swapj(cor1, left, (left+right)/2)
  call swapj(cor2, left, (left+right)/2)
  last = left
  do i=left+1,right
     if (wa(i) > wa(left)) then
        last=last+1
        call swapm(wa, last, i)!swap submat
        call swapj(cor1, last, i)
        call swapj(cor2, last, i)
     end if
  end do
  call swapm(wa, left, last)!swap submat
  call swapj(cor1, left, last)
  call swapj(cor2, left, last)
  itmp=last-1
  call qsortR2I(wa, cor1, cor2, left, itmp)
  itmp=last+1
  call qsortR2I(wa, cor1, cor2, itmp, right)
  return
end subroutine qsortR2I


!|---------------------------------------------------------------------*/
subroutine swap_submat(w,left,right,nsub)
  integer:: nsub,left,right
  real(dp),intent(inout),dimension(:,:,:)::w
  integer iisub,jjsub
  real(dp),dimension(nsub,nsub)::s
  do iisub=1,nsub
     do jjsub=1,nsub
        s(iisub,jjsub)=w(iisub,jjsub,left)
        w(iisub,jjsub,left)=w(iisub,jjsub,right)
        w(iisub,jjsub,right)=s(iisub,jjsub)
     end do
  end do
  return
end subroutine swap_submat
!|---------------------------------------------------------------------*/


!void qsortR2I(double *wa, int *cor1, int *cor2, int left, int right){
recursive subroutine qsortR2I_submat(wa, cor1, cor2, left, right, nsubrow)
!/*----------------------------------------------------------------------
!|
!| qqsort: sort wa[left]...wa[right] into decreasing order
!| from Kernighan & Ritchie
!|
!|---------------------------------------------------------------------*/
  real(dp),intent(inout),dimension(:,:,:)::wa
  integer,intent(inout),dimension(:)::cor1,cor2
  integer,intent(in):: left,right,nsubrow
  !internal vars
  integer:: i, last,itmp
  if (left >= right) return
  call swap_submat(wa,last,(left+right)/2,nsubrow)
  call swapj(cor1, left, (left+right)/2)
  call swapj(cor2, left, (left+right)/2)
  last = left
  do i=left+1,right

     if (L1Norm2(wa,i,nsubrow) > L1Norm2(wa,left,nsubrow)) then
        last=last+1
        call swap_submat(wa,last,i,nsubrow)
        call swapj(cor1, last, i)
        call swapj(cor2, last, i)
     end if
  end do
  call swap_submat(wa,left,last,nsubrow)
  call swapj(cor1, left, last)
  call swapj(cor2, left, last)
  itmp=last-1
  call qsortR2I_submat(wa, cor1, cor2, left, itmp,nsubrow)
  itmp=last+1
  call qsortR2I_submat(wa, cor1, cor2, itmp, right,nsubrow)
  return
end subroutine qsortR2I_submat




!void qsort2C(int *ja, double *ma, int left, int right, int abval){
recursive subroutine qsort2C(ja,ma,left,right,abval)
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
  real(dp),intent(inout),dimension(:)::ma
  integer,intent(inout),dimension(:)::ja
  integer,intent(inout) :: left,right,abval
  !internal vars
  integer :: i, last,itmp
  if (left >= right)  return
  if (abval.gt.0) then
     itmp=(left+right)/2
     call swapj(ja, left, itmp)
     call swapm(ma, left, itmp)
     last = left
     do i=left+1,right
        if (abs(ma(i)) < abs(ma(left))) then
           last=last+1
           call swapj(ja, last, i)
           call swapm(ma, last, i)
        end if
     end do
     call swapj(ja, left, last)
     call swapm(ma, left, last)
     itmp=last-1
     call qsort2C(ja, ma, left, itmp, abval)
     itmp=last+1
     call qsort2C(ja, ma, itmp, right, abval)
  else 
     itmp=(left+right)/2
     call swapj(ja, left, itmp)
     call swapm(ma, left, itmp)
     last = left
     do i=left+1,right 
        if (ma(i) < ma(left)) then
           last=last+1
           call swapj(ja, last, i)
           call swapm(ma, last, i)
        end if
     end do
     call swapj(ja, left, last)
     call swapm(ma, left, last)
     itmp=last-1
     call qsort2C(ja, ma, left, itmp, abval)
     itmp=last+1
     call qsort2C(ja, ma, itmp, right, abval)
  end if
  return
end subroutine qsort2C






  
recursive subroutine qsortC(ja, ma, left, right, abval)
!/*----------------------------------------------------------------------
!|
!| qqsort: sort ma[left]...ma[right] into decreasing order
!| from Kernighan & Ritchie
!|
!| ja holds the column indices
!| abval = 1: consider absolute values
!|         0: values
!|
!|---------------------------------------------------------------------*/
  real(dp),intent(inout),dimension(:)::ma
  integer,intent(inout),dimension(:)::ja
  integer,intent(inout) :: left,right,abval
  !internal vars
  integer :: i, last,itmp
  if (left >= right)  return
  if (abval.gt.0) then
     itmp=(left+right)/2
     call swapj(ja, left, itmp)
     call swapm(ma, left, itmp)
     last = left
     do i=left+1,right
        if (abs(ma(i)) > abs(ma(left))) then
           last=last+1
           call swapj(ja, last, i)
           call swapm(ma, last, i)
        end if
     end do
     call swapj(ja, left, last)
     call swapm(ma, left, last)
     itmp=last-1
     call qsortC(ja, ma, left, itmp, abval)
     itmp=last+1
     call qsortC(ja, ma, itmp, right, abval)
  else 
     itmp=(left+right)/2
     call swapj(ja, left, itmp)
     call swapm(ma, left, itmp)
     last = left
     do i=left+1,right
        if (ma(i) > ma(left)) then
           last=last+1
           call swapj(ja, last, i)
           call swapm(ma, last, i)
        end if
     end do
    call swapj(ja, left, last)
    call swapm(ma, left, last)
    itmp=last-1
    call qsortC(ja, ma, left, itmp, abval)
    itmp=last+1
    call qsortC(ja, ma, itmp, right, abval)
 end if
 return
end subroutine qsortC





subroutine swapj(v,i,j)
  integer,intent(inout),dimension(:) :: v
  integer :: i,j
  integer :: temp
  temp = v(i)
  v(i) = v(j)
  v(j) = temp
end subroutine swapj

subroutine swapm(m,i,j) 
  real(dp),intent(inout),dimension(:) :: m
  integer,intent(in) :: i,j
  real(dp) :: temp
  temp = m(i)
  m(i) = m(j)
  m(j) = temp
end subroutine swapm

subroutine swapmat(m,i,j,nsub) 
  real(dp),intent(inout),dimension(:,:,:) :: m
  integer,intent(in) :: i,j,nsub
  real(dp),dimension(nsub,nsub) :: temp
  integer :: isub,jsub
  do isub=1,nsub
     do jsub=1,nsub
        temp(isub,jsub) = m(isub,jsub,i)
        m(isub,jsub,i)  = m(isub,jsub,j)
        m(isub,jsub,j)  = temp(isub,jsub)
     end do
  end do
end subroutine swapmat

!==================================================================================================
end module blas


function tmemcpy_r(a,b,size) result(retval)
  use kinddefs, only : dp
  real(dp)::a(size), b(size)
  integer::size,type
  integer::i,retval
  do i=1,size
     a(i)=b(i)
  end do
  retval=0
  return
end function tmemcpy_r

function tmemcpy_i(a,b,size) result(retval)
  integer::a(size), b(size)
  integer::size,type
  integer::i,retval
  do i=1,size
     a(i)=b(i)
  end do
  retval=0
  return
end function tmemcpy_i

function tblank_r(a,size) result(retval)
  use kinddefs, only : dp
  real(dp)::a(size)
  integer::size,type
  integer::i,retval
  do i=1,size
     a(i)=0
  end do
  retval=0
  return
end function tblank_r

function tblank_i(a,size) result(retval)
  integer::a(size)
  integer::size,type
  integer::i,retval
  do i=1,size
     a(i)=0
  end do
  retval=0
  return
end function tblank_i
