program dog
  parameter(maxnnz=88000,maxnodes=5500,maxnsub=6)
  dimension A(maxnsub*maxnsub*maxnnz),qnode(maxnodes*maxnsub)
  integer ia(maxnodes),iau(maxnodes),ja(maxnnz)
  integer ianew(maxnodes),iaunew(maxnodes),janew(maxnnz)
  integer list(maxnodes),tag(maxnodes)
  
  
  character(len=130)::matrix_file
  matrix_file="matrix.txt"
  !trim(matrix_file)
  open(UNIT=5,FILE=matrix_file,form='formatted',STATUS='old')
  !
  ! Read the input
  !
  read(5,*) !Dimension
  read(5,*)nnodes
  read(5,*) !block size
  read(5,*)nsub
  read(5,*)!number of non zero entries in matrix
  read(5,*)nnz
  
  close(5)
  
  iadim = nnodes + 1
  
  write(6,*)"About to call reader"
  !
  ! Notice that I am using the tag array twice
  !
  call READER(nnodes,nsub,nnz,A,ia,iau,ja,list,tag,ianew,iaunew,janew,tag,qnode)
  
  stop
end program dog



!============================ BANDWIDTH ================================
!
! Finds the maximum bandwidth
!
!===================================================================
subroutine BANDWIDTH(ia,iau,ja,nnodes,nnz)
  integer ia(nnodes + 1),iau(nnodes),ja(nnz)
  integer lowerBand,upperBand
  integer bw
  !
  bw = 0
  do  i = 1,nnodes
     jstart = ia(i)
     jend   = ia(i+1) - 1
     lowerBand  = i - ja(jstart)
     upperBand = ja(jend) - i
     bw = max(lowerBand,upperBand,bw)
  end do
  write(6,*)"Bandwidth = ",bw
  return
end subroutine BANDWIDTH


!============================ READER ================================
!
! Reads the matrix in Chad's format and converts to mine
! The only difference is that he reads in nnodes+1 iau's
! This routine also gives the option to transpose 
! or divide by the diagonal
!
!===================================================================
subroutine READER(nnodes,nsub,nnz,A,ia,iau,ja,list,tag,ianew,iaunew,janew,oldtonew,qnode)
  dimension A(nsub,nsub,nnz),qnode(nsub,nnodes)
  integer ia(nnodes + 1),iau(nnodes),ja(nnz)
  integer ianew(nnodes + 1),iaunew(nnodes),janew(nnz)
  integer list(nnodes),tag(nnodes)
  integer oldrow,newrow,oldcol,newcol
  integer k,m,ii,jj
  real s(nsub)
  
  character(len=130)::matrix_file
  matrix_file="matrix.txt"
  !trim(matrix_file)
  open(UNIT=5,FILE=matrix_file,form='formatted',STATUS='old')
  
  
  read(5,*) !Dimension
  read(5,*)nnodes
  read(5,*) !block size
  read(5,*)nsub
  read(5,*)!number of non zero entries in matrix
  read(5,*)nnz
  
  read(5,*)
  do  i = 1,nnodes+1
     read(5,*) ia(i)
  end do
  read(5,*)
  do  i = 1,nnodes
     read(5,*) iau(i)
  end do
  !     read(5,*) ! Chad has nnodes+1 iau
  read(5,*)
  do  i = 1,nnz
     read(5,*) ja(i)
  end do
  read(5,*)
  do  i = 1,nnz
     !read in column major
     do j=1,nsub
        read(5,*) (A(j,k,i),k=1,nsub)
     end do
  end do
  close(5)
  
  
  
  ! Write out the matrix and the rhs so the answer is all 1's
  !
  ! First transpose A
  !
  !     call transposeA(nnodes,nnz,ia,ja,iau,A)
  !     call scaleAD(nnodes,nnz,ia,ja,iau,A)
  !
  ! Generate RHS. Use qnode for scratch
  !
  iread = 0
  if(iread.eq.0)then
     do  i = 1,nnodes
        do j=1,nsub
           qnode(j,i) = 0.
        end do
     end do
     do  i = 1,nnodes
        jstart = ia(i) ! Beginning of row
        jend   = ia(i+1) - 1 ! Ending of row
        do k=1,nsub
           s(k)=0.
        end do
        do  j=jstart,jend
           icol   = ja(j) ! We don't care about the column because we have all 1's as the solution
           do k=1,nsub
              do m=1,nsub
                 s(k) =  s(k) + A(k,m,j)*1.  
              end do
           end do
        end do
        do k=1,nsub
           qnode(k,i) = qnode(k,i) + s(k)
        end do
     end do
  else
     read(5,*)
     do i = 1,nnodes
        do k=1,nsub
           read(5,*)qnode(k,i) 
        end do
     end do
  end if
  
  
  ! If you want to read in a right hand side, do it here
  !       read(5,*)
  !       do i = 1,nnodes
  !         read(5,*)(qnode(j,i),j=1,4)
  !         A(1,1,iau(i)) = A(1,1,iau(i)) + 1.
  !         A(2,2,iau(i)) = A(2,2,iau(i)) + 1.
  !         A(3,3,iau(i)) = A(3,3,iau(i)) + 1.
  !         A(4,4,iau(i)) = A(4,4,iau(i)) + 1.
  !       end do
  ! 
  ! We are finished generating the right-hand side
  ! Now, modify A to give us more diagonal dominance
  !
  !     call setDDA(A,nnodes,nnz,ia,ja,iau)
  !
  ! Re-order the matrix using Cuthill-McKee
  call BANDWIDTH(ia,iau,ja,nnodes,nnz)
  call reOrder(nnodes,ndim,nnz,ia,ja,iau,A,list,tag,ianew,iaunew,janew,oldtonew)
  call BANDWIDTH(ianew,iaunew,janew,nnodes,nnz)
  
  
  
  matrix_file="matrix.bin"
  !trim(matrix_file)
  open(UNIT=299,FILE=matrix_file,form='unformatted',STATUS='replace')
  
  !
  ! Write out the matrix
  !
  !write(299)"Dimension"
  write(99,*)"Dimension"
  write(99,*)nnodes
  write(299)nnodes

  !write(299)"block"
  write(99,*)"block"
  write(99,*)nsub
  write(299)nsub

  !write(299)" nnz"
  write(99,*)" nnz"
  write(99,*)nnz
  write(299)nnz

  !write(299)"IA dimension=",nnodes+1
  write(99,*)"IA dimension=",nnodes+1
  do  i = 1,nnodes+1
     write(99,*)ianew(i)
     write(299)ianew(i)
     !write(99,*)ia(i)
     !write(299)ia(i)
  end do
  !write(299)"IAU dimension=",nnodes
  write(99,*)"IAU dimension=",nnodes
  do  i = 1,nnodes
     write(99,*)iaunew(i)
     write(299)iaunew(i)
     !write(99,*)iau(i)
     !write(299)iau(i)
  end do
  !write(299)"JA nnz=",nnz
  write(99,*)"JA nnz=",nnz
  do  i = 1,nnz
     write(99,*)janew(i)
     write(299)janew(i)
     !write(99,*)ja(i)
     !write(299)ja(i)
  end do
  !write(299)"A dimension=",nsub*nnz," (nsub*nnz)"
  write(99,*)"A dimension=",nsub*nnz," (nsub*nnz)"
  do  i = 1,nnodes
     newrow = i
     oldrow = list(i)
     jstart = ianew(newrow)
     jend   = ianew(newrow +1) - 1
     kstart = ia(oldrow)
     kend   = ia(oldrow +1) - 1
     do  j = jstart,jend
        newcol = list(janew(j)) ! This is the old number for the new column
        do  k = kstart,kend
           oldcol = ja(k)
           if(newcol.eq.oldcol)then
              do ii=1,nsub
                 write(99,*)(A(ii,jj,k),jj=1,nsub) 
                 write(299)(A(ii,jj,k),jj=1,nsub) 
              end do
           end if
        end do
     end do
  end do
  
  !
  ! Write out the right-hand side
  !
  !write(299)"Right hand side dimension=",nnodes
  write(99,*)"Right hand side dimension=",nnodes
  do i = 1,nnodes
     index = list(i)
     write(99,*)(qnode(j,index),j=1,nsub)
     write(299)(qnode(j,index),j=1,nsub)
  end do
  
  close(299)
  return
end subroutine READER


!============================ SETDDA ================================
!
! Compute some information related to the diagonal dominance of A
! We are looking to monitor D/sum(off-diagonals)
!
!===================================================================
subroutine setDDA(A,nnodes,ndim,nnz,ia,ja,iau)
  dimension A(ndim,ndim,nnz)
  integer ia(1),ja(1),iau(1)
  integer i,j,k,ii,jj
  real d(ndim)
  real ratio(ndim)
  real ratmin(ndim)
  real ratave(ndim)
  real sum(ndim)
  real tol
  ! This is code to examine the diagonal dominance
  !
  do i=1,ndim
     ratmin(i)=1.e9
     ratave(i)=0.0
  end do
  !     write(6,*)"Enter tol"
  !     read(5,*)tol
  do  i = 1,nnodes
     ii = iau(i)
     do j=1,ndim
        d(j) = abs(A(j,j,ii))
     end do
     !write(67,202)i,A(1,1,ii),A(2,2,ii),A(3,3,ii),A(4,4,ii)
     !
     ! Sum off-diagonals from diagonal block
     !
     do j=1,ndim
        do k=1,ndim
           if(k .ne. j) sum(j) = sum(j) + abs(A(j,k,ii))
        end do
     end do
     !
     ! Now get off-diagonals from non-diagonal blocks
     !
     jstart = ia(i) ! Beginning of row
     jend   = ia(i+1) - 1 ! Ending of row
     do  j = jstart,jend
        icol = ja(j)
        if(icol.ne.i)then
           do k=1,ndim
              do jj=1,ndim
                 sum(k) = sum(k) + abs(A(k,jj,j))
              end do
           end do
        end if
     end do
     
     tol = 0.00
     do k=1,ndim
        ratio(k)=d(k)/sum(k)
        ratave(k)=ratave(k)+ratio(k)
        if(ratio(k).le.ratmin(k))ratmin(k)=ratio(k)
        if(ratio(k).le.tol)then
           if(A(k,k,ii).lt.0.0)A(k,k,ii) = -tol*sum(k)
           if(A(k,k,ii).gt.0.0)A(k,k,ii) =  tol*sum(k)
        end if
     end do
     
     !        write(66,202)i,ratio(1),ratio(2),ratio(3),ratio(4)
202  format(1h ,i6,1x,4(e14.7,1x))
  end do
  
  write(6,*)"Diagonal dominance indicators"
  do k=1,ndim
     ratave(k) = ratave(k)/nnodes
     write(6,*)"Ratio ",k," = ",ratmin(k)," average = ",ratave(k)
  end do
  
  return
end subroutine setDDA


!============================ SCALEAD ==============================
!
! Scale A by the diagonal
!
!===================================================================
subroutine scaleAD(nnodes,ndim,nnz,ia,ja,iau,A)
  !parameter( neqn= 4)
  dimension A(ndim,ndim,nnz)
  dimension rowEntry(ndim,ndim),columnEntry(ndim,ndim)
  integer ia(nnodes + 1),ja(nnz),iau(nnodes)
  integer column,i,j,k,ii,jj
  real d(ndim)
  !
  do  i = 1,nnodes
     do ii=1,ndim
        d(ii) = A(ii,ii,iau(i))
     end do
     
     !     dmin = max(d1,d2,d3,d4)
     !     d1 = dmin
     !     d2 = dmin
     !     d3 = dmin
     !     d4 = dmin
     jstart = ia(i)
     jend   = ia(i+1) - 1
     do  j = jstart,jend
        do ii=1,ndim
           do jj=1,ndim
              A(ii,jj,j) = A(ii,jj,j)/d(ii)
           end do
        end do
     end do
  end do
  return
end subroutine scaleAD


!============================ REORDER ==============================
!
! Re-order A using Cuthill-McKee
!
!===================================================================
subroutine reOrder(nnodes,ndim,nnz,ia,ja,iau,A,list,tag,ianew,iaunew,janew,oldtonew)
  dimension A(ndim,ndim,nnz)
  integer ia(nnodes+1),iau(nnodes),ja(nnz)
  integer ianew(nnodes+1),iaunew(nnodes),janew(nnz)
  integer list(nnodes),tag(nnodes)
  integer oldnode,oldcol,newcol
  integer oldtonew(nnodes)
  !
  do  i = 1,nnodes
     list(i) = i
     iaunew(i) = iau(i)
     ianew(i)  = ia(i)
  end do
  list(nnodes+1) = nnodes + 1
  ianew(nnodes+1) = ia(nnodes+1)
  do  i = 1,nnz
     janew(i) = ja(i)
  end do
  ireorder = 0 
  if(ireorder.ne.1)return
  !
  do  i = 1,nnodes
     tag(i) = 0
  end do
  !
  ! Now go grab all the entries that connect to this node in the number list
  !
  !     write(6,*)"Entering reOrder"
  icount = 1
  istart = 1
  list(1) = 1
  tag(1)  = 1
6000 continue
  jstart = ia(list(istart))
  jend   = ia(list(istart)+ 1) - 1
  do  j = jstart,jend
     inode = ja(j)
     if(tag(inode).eq.0)then
        icount = icount + 1
        list(icount) = inode
        tag(inode) = 1
        if(icount.eq.nnodes)goto 7000
     end if
  end do
  istart = istart + 1
  goto 6000
7000 continue
  list(nnodes+1) = nnodes + 1
  !     write(6,*)"Exiting reOrder"
  !
  ! At this point, list(i) is a list of the old node
  ! numbers stored in the new order
  ! For example, if list() = {1,3,4,6,5,2}
  ! then list(1) = old node 1
  !      list(2) = old node 3
  !      list(3) = old node 4
  !      list(4) = old node 6
  !      list(5) = old node 5
  !      list(6) = old node 2
  !
  ! Now we can construct new ia, iau, and ja
  !
  ianew(1) = 1
  do  i = 1,nnodes
     oldnode = list(i)
     nonzeros = ia(oldnode + 1) - ia(oldnode)
     ianew(i+1) = ianew(i) + nonzeros
  end do
  !
  ! Fill oldtonew
  !
  do  i = 1,nnodes
     oldtonew(list(i)) = i
  end do
  !
  ! Now get ja
  !
  jcount = 0
  do  i = 1,nnodes
     oldnode = list(i)
     jstart = ia(oldnode)
     jend   = ia(oldnode + 1) - 1
     do  j = jstart,jend
        oldcol = ja(j)
        newcol = oldtonew(oldcol)
        jcount = jcount + 1
        janew(jcount) = newcol
     end do
     istart = ianew(i)
     iend   = ianew(i+1) - 1
     call SORTER(istart,iend,janew,iaunew,i)
  end do
  !
  ! Now get iau
  !
  !     do 8040 i = 1,nnodes
  !       jstart = ianew(i)
  !       jend   = ianew(i+1) - 1
  !       do 8050 j = jstart,jend
  !         if(janew(j).eq.i)iaunew(i) = j
  !8050   continue
  !8040 continue
  iaunew(nnodes + 1) = -1
  !
  !     write(6,*)"New ia and iau"
  !     do 9000 i = 1,nnodes+1
  !       write(6,*)ianew(i),iaunew(i)
  !9000 continue
  !
  !     write(6,*)"New ja"
  !     do 9001 i = 1,nnodes
  !       jstart = ianew(i)
  !       jend = ianew(i+1) - 1
  !       write(6,*)(janew(j),j=jstart,jend)
  !9001 continue
  !     stop
  !
  return
end subroutine reOrder


!===================================================================
!
! Sort each of our bins
!
!===================================================================
subroutine SORTER(istart,iend,ja,iau,inode)
  integer ja(1),iau(1)
  !
  do  i = istart,iend
     min = ja(i)
     minsave = ja(i)
     jsave = i
     do  j = i+1,iend
        if(ja(j).lt.min)then
           min = ja(j)
           jsave = j
        end if
     end do
     ja(i) = min
     ja(jsave) = minsave
     if(ja(i).eq.inode)iau(inode) = i
  end do
  !
  return
end subroutine SORTER


!============================ TRANSPOSEA ===========================
!
! Transpose A
!
!===================================================================
subroutine transposeA(nnodes,ndim,nnz,ia,ja,iau,A)
  dimension A(ndim,ndim,nnz)
  dimension rowEntry(ndim,ndim),columnEntry(ndim,ndim)
  integer ia(nnodes + 1),ja(nnz),iau(nnodes)
  integer column
  !
  do  i = 1,nnodes
     jstart = iau(i)      ! Start on diagonal
     jend   = ia(i+1) - 1 ! Go to the end of the row
     do  j = jstart,jend
        column = ja(j)
        do  k = 1,ndim
           do  L = 1,ndim
              rowEntry(k,L) = A(k,L,j)
           end do
        end do
        kstart = ia(column)
        kend   = ia(column+1) - 1
        do  k = kstart,kend
           kindex = k
           if(ja(k).eq.i)goto 1050
        end do
1050    continue
        do  k = 1,ndim
           do  L = 1,ndim
              columnEntry(k,L) = A(L,k,kindex)
              !                 A(k,L,j) = columnEntry(k,L)
              !                 A(k,L,kindex) = rowEntry(L,k)
           end do
        end do
        do  k = 1,ndim
           do  L = 1,ndim
              A(k,L,j) = columnEntry(k,L)
              A(k,L,kindex) = rowEntry(L,k)
           end do
        end do
     end do
  end do
  return
end subroutine transposeA
