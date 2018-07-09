      program dog
      parameter(maxnnz=2000000,maxnodes=95000)
      dimension A(4*4*maxnnz),qnode(maxnodes*4)
      integer ia(maxnodes),iau(maxnodes),ja(maxnnz)
      integer ianew(maxnodes),iaunew(maxnodes),janew(maxnnz)
      integer list(maxnodes),tag(maxnodes)
c
c Read the input
c
      read(5,*)
      read(5,*)nnodes
      read(5,*)
      read(5,*)nnz
      iadim = nnodes + 1
c
c Notice that I am using the tag array twice
c
      call READER(nnodes,nnz,A,ia,iau,ja,list,tag,
     1            ianew,iaunew,janew,tag)

      stop
      end
c============================ BANDWIDTH ================================
c
c Finds the maximum bandwidth
c
c===================================================================
      subroutine BANDWIDTH(ia,iau,ja,nnodes,nnz)
      integer ia(nnodes + 1),iau(nnodes),ja(nnz)
      integer lowerBand,upperBand
      integer bw
c
      bw = 0
      do 1000 i = 1,nnodes
        jstart = ia(i)
        jend   = ia(i+1) - 1
        lowerBand  = i - ja(jstart)
        upperBand = ja(jend) - i
        bw = max(lowerBand,upperBand,bw)
 1000 continue
      write(6,*)"Bandwidth = ",bw
      return
      end
c============================ READER ================================
c
c Reads the matrix in Chad's format and converts to mine
c The only difference is that he reads in nnodes+1 iau's
c This routine also gives the option to transpose 
c or divide by the diagonal
c
c===================================================================
      subroutine READER(nnodes,nnz,A,ia,iau,ja,list,tag,
     1            ianew,iaunew,janew,oldtonew)
      dimension A(4,4,nnz),qnode(4,nnodes)
      integer ia(nnodes + 1),iau(nnodes),ja(nnz)
      integer ianew(nnodes + 1),iaunew(nnodes),janew(nnz)
      integer list(nnodes),tag(nnodes)
      integer oldrow,newrow,oldcol,newcol

      read(5,*)
      do 1000 i = 1,nnodes+1
        read(5,*)ia(i)
 1000 continue
      read(5,*)
      do 1010 i = 1,nnodes
        read(5,*)iau(i)
 1010 continue
c     read(5,*) ! Chad has nnodes+1 iau
      read(5,*)
      do 1020 i = 1,nnz
        read(5,*)ja(i)
 1020 continue
      read(5,*)
      do 1030 i = 1,nnz
        read(5,*)A(1,1,i),A(1,2,i),A(1,3,i),A(1,4,i)
        read(5,*)A(2,1,i),A(2,2,i),A(2,3,i),A(2,4,i)
        read(5,*)A(3,1,i),A(3,2,i),A(3,3,i),A(3,4,i)
        read(5,*)A(4,1,i),A(4,2,i),A(4,3,i),A(4,4,i)
 1030 continue
c
c Write out the matrix and the rhs so the answer is all 1's
c
c First transpose A
c
c     call transposeA(nnodes,nnz,ia,ja,iau,A)
c     call scaleAD(nnodes,nnz,ia,ja,iau,A)
c
c Generate RHS. Use qnode for scratch
c
      do 5030 i = 1,nnodes
        qnode(1,i) = 0.
        qnode(2,i) = 0.
        qnode(3,i) = 0.
        qnode(4,i) = 0.
 5030 continue
      do 5040 i = 1,nnodes
        jstart = ia(i) ! Beginning of row
        jend   = ia(i+1) - 1 ! Ending of row
        s1 = 0.
        s2 = 0.
        s3 = 0.
        s4 = 0.
        do 5050 j=jstart,jend
        icol   = ja(j) ! We don't care about the column because we have all 1's as the solution
        s1 =  s1 + A(1,1,j)*1. + A(1,2,j)*1. + A(1,3,j)*1. + A(1,4,j)*1.
        s2 =  s2 + A(2,1,j)*1. + A(2,2,j)*1. + A(2,3,j)*1. + A(2,4,j)*1.
        s3 =  s3 + A(3,1,j)*1. + A(3,2,j)*1. + A(3,3,j)*1. + A(3,4,j)*1.
        s4 =  s4 + A(4,1,j)*1. + A(4,2,j)*1. + A(4,3,j)*1. + A(4,4,j)*1.
 5050   continue
        qnode(1,i) = qnode(1,i) + s1
        qnode(2,i) = qnode(2,i) + s2
        qnode(3,i) = qnode(3,i) + s3
        qnode(4,i) = qnode(4,i) + s4
 5040 continue
c 
c We are finished generating the right-hand side
c Now, modify A to give us more diagonal dominance
c
c     call setDDA(A,nnodes,nnz,ia,ja,iau)
c
c Re-order the matrix using Cuthill-McKee
      call BANDWIDTH(ia,iau,ja,nnodes,nnz)
      call reOrder(nnodes,nnz,ia,ja,iau,A,list,tag,
     1             ianew,iaunew,janew,oldtonew)
      call BANDWIDTH(ianew,iaunew,janew,nnodes,nnz)
c
c Write out the matrix
c
      write(99,*)"Dimension"
      write(99,*)nnodes
      write(299)nnodes
      write(99,*)" nnz"
      write(99,*)nnz
      write(299)nnz
      write(99,*)"IA dimension=",nnodes+1
      do 5000 i = 1,nnodes+1
        write(99,*)ianew(i)
        write(299)ianew(i)
 5000 continue
      write(99,*)"IAU dimension=",nnodes
      do 5010 i = 1,nnodes
      write(99,*)iaunew(i)
      write(299)iaunew(i)
 5010 continue
      write(99,*)"JA nnz=",nnz
      do 5080 i = 1,nnz
         write(99,*)janew(i)
         write(299)janew(i)
 5080 continue
      write(99,*)"A dimension=",4*nnz," (4*nnz)"
      do 6000 i = 1,nnodes
        newrow = i
        oldrow = list(i)
        jstart = ianew(newrow)
        jend   = ianew(newrow +1) - 1
        kstart = ia(oldrow)
        kend   = ia(oldrow +1) - 1
        do 6010 j = jstart,jend
          newcol = list(janew(j)) ! This is the old number for the new column
          do 6020 k = kstart,kend
            oldcol = ja(k)
            if(newcol.eq.oldcol)then
              write(99,*)A(1,1,k),A(1,2,k),A(1,3,k),A(1,4,k)
              write(99,*)A(2,1,k),A(2,2,k),A(2,3,k),A(2,4,k)
              write(99,*)A(3,1,k),A(3,2,k),A(3,3,k),A(3,4,k)
              write(99,*)A(4,1,k),A(4,2,k),A(4,3,k),A(4,4,k)
              write(299)A(1,1,k),A(1,2,k),A(1,3,k),A(1,4,k)
              write(299)A(2,1,k),A(2,2,k),A(2,3,k),A(2,4,k)
              write(299)A(3,1,k),A(3,2,k),A(3,3,k),A(3,4,k)
              write(299)A(4,1,k),A(4,2,k),A(4,3,k),A(4,4,k)
c             write(99,*)"----------------------------"
c             kk = k + 1
            end if
c           kstart = kk
 6020     continue
 6010   continue
 6000 continue
c
c     do 5020 i = 1,nnz
c     write(99,*)A(1,1,i),A(1,2,i),A(1,3,i),A(1,4,i)
c     write(99,*)A(2,1,i),A(2,2,i),A(2,3,i),A(2,4,i)
c     write(99,*)A(3,1,i),A(3,2,i),A(3,3,i),A(3,4,i)
c     write(99,*)A(4,1,i),A(4,2,i),A(4,3,i),A(4,4,i)
c5020 continue
c
c Write out the right-hand side
c
      write(99,*)"Right hand side dimension=",nnodes
      do 5060 i = 1,nnodes
        index = list(i)
        write(99,*)(qnode(j,index),j=1,4)
        write(299)(qnode(j,index),j=1,4)
 5060 continue
      
      return
      end
c============================ SETDDA ================================
c
c Compute some information related to the diagonal dominance of A
c We are looking to monitor D/sum(off-diagonals)
c
c===================================================================
      subroutine setDDA(A,nnodes,nnz,ia,ja,iau)
      dimension A(4,4,nnz)
      integer ia(1),ja(1),iau(1)

c
c This is code to examine the diagonal dominance
c
      rat1min = 1.e9
      rat2min = 1.e9
      rat3min = 1.e9
      rat4min = 1.e9
      rat1ave = 0.0
      rat2ave = 0.0
      rat3ave = 0.0
      rat4ave = 0.0
c     write(6,*)"Enter tol"
c     read(5,*)tol
      do 7000 i = 1,nnodes
        ii = iau(i)
        d1 = abs(A(1,1,ii))
        d2 = abs(A(2,2,ii))
        d3 = abs(A(3,3,ii))
        d4 = abs(A(4,4,ii))
        write(67,202)i,A(1,1,ii),A(2,2,ii),A(3,3,ii),A(4,4,ii)
c
c Sum off-diagonals from diagonal block
c
        sum1 = abs(A(1,2,ii)) + abs(A(1,3,ii)) + abs(A(1,4,ii))
        sum2 = abs(A(2,1,ii)) + abs(A(2,3,ii)) + abs(A(2,4,ii))
        sum3 = abs(A(3,1,ii)) + abs(A(3,2,ii)) + abs(A(3,4,ii))
        sum4 = abs(A(4,1,ii)) + abs(A(4,2,ii)) + abs(A(4,3,ii))
c
c Now get off-diagonals from non-diagonal blocks
c
        jstart = ia(i) ! Beginning of row
        jend   = ia(i+1) - 1 ! Ending of row
        do 7001 j = jstart,jend
          icol = ja(j)
          if(icol.ne.i)then
             sum1 = sum1 + abs(A(1,1,j)) + abs(A(1,2,j)) 
     1                   + abs(A(1,3,j)) + abs(A(1,4,j))
             sum2 = sum2 + abs(A(2,1,j)) + abs(A(2,2,j)) 
     1                   + abs(A(2,3,j)) + abs(A(2,4,j))
             sum3 = sum3 + abs(A(3,1,j)) + abs(A(3,2,j)) 
     1                   + abs(A(3,3,j)) + abs(A(3,4,j))
             sum4 = sum4 + abs(A(4,1,j)) + abs(A(4,2,j)) 
     1                   + abs(A(4,3,j)) + abs(A(4,4,j))
          end if
 7001   continue
        ratio1 = d1/sum1
        ratio2 = d2/sum2
        ratio3 = d3/sum3
        ratio4 = d4/sum4
        rat1ave = rat1ave + ratio1
        rat2ave = rat2ave + ratio2
        rat3ave = rat3ave + ratio3
        rat4ave = rat4ave + ratio4
        if(ratio1.le.rat1min)rat1min = ratio1
        if(ratio2.le.rat2min)rat2min = ratio2
        if(ratio3.le.rat3min)rat3min = ratio3
        if(ratio4.le.rat4min)rat4min = ratio4

        tol = 0.00
        if(ratio1.le.tol)then
c          A(1,1,ii) = tol*sum1
           if(A(1,1,ii).lt.0.0)A(1,1,ii) = -tol*sum1
           if(A(1,1,ii).gt.0.0)A(1,1,ii) =  tol*sum1
        end if
        if(ratio2.le.tol)then
c          A(2,2,ii) = tol*sum2
           if(A(2,2,ii).lt.0.0)A(2,2,ii) = -tol*sum2
           if(A(2,2,ii).gt.0.0)A(2,2,ii) =  tol*sum2
        end if
        if(ratio3.le.tol)then
c          A(3,3,ii) = tol*sum3
           if(A(3,3,ii).lt.0.0)A(3,3,ii) = -tol*sum3
           if(A(3,3,ii).gt.0.0)A(3,3,ii) =  tol*sum3
        end if
        if(ratio4.le.tol)then
c          A(4,4,ii) = tol*sum4
           if(A(4,4,ii).lt.0.0)A(4,4,ii) = -tol*sum4
           if(A(4,4,ii).gt.0.0)A(4,4,ii) =  tol*sum4
        end if

        write(66,202)i,ratio1,ratio2,ratio3,ratio4
  202   format(1h ,i6,1x,4(e14.7,1x))
 7000 continue
      rat1ave = rat1ave/nnodes
      rat2ave = rat2ave/nnodes
      rat3ave = rat3ave/nnodes
      rat4ave = rat4ave/nnodes
      write(6,*)"Diagonal dominance indicators"
      write(6,*)"Ratio 1 = ",rat1min," average = ",rat1ave
      write(6,*)"Ratio 2 = ",rat2min," average = ",rat2ave
      write(6,*)"Ratio 3 = ",rat3min," average = ",rat3ave
      write(6,*)"Ratio 4 = ",rat4min," average = ",rat4ave
c     stop
      return
      end
c============================ SCALEAD ==============================
c
c Scale A by the diagonal
c
c===================================================================
      subroutine scaleAD(nnodes,nnz,ia,ja,iau,A)
      parameter( neqn= 4)
      dimension A(neqn,neqn,nnz)
      dimension rowEntry(neqn,neqn),columnEntry(neqn,neqn)
      integer ia(nnodes + 1),ja(nnz),iau(nnodes)
      integer column
c
      do 1000 i = 1,nnodes
      d1 = A(1,1,iau(i))
      d2 = A(2,2,iau(i))
      d3 = A(3,3,iau(i))
      d4 = A(4,4,iau(i))

c     dmin = max(d1,d2,d3,d4)
c     d1 = dmin
c     d2 = dmin
c     d3 = dmin
c     d4 = dmin
      jstart = ia(i)
      jend   = ia(i+1) - 1
      do 1010 j = jstart,jend
         A(1,1,j) = A(1,1,j)/d1
         A(1,2,j) = A(1,2,j)/d1
         A(1,3,j) = A(1,3,j)/d1
         A(1,4,j) = A(1,4,j)/d1
         A(2,1,j) = A(2,1,j)/d2
         A(2,2,j) = A(2,2,j)/d2
         A(2,3,j) = A(2,3,j)/d2
         A(2,4,j) = A(2,4,j)/d2
         A(3,1,j) = A(3,1,j)/d3
         A(3,2,j) = A(3,2,j)/d3
         A(3,3,j) = A(3,3,j)/d3
         A(3,4,j) = A(3,4,j)/d3
         A(4,1,j) = A(4,1,j)/d4
         A(4,2,j) = A(4,2,j)/d4
         A(4,3,j) = A(4,3,j)/d4
         A(4,4,j) = A(4,4,j)/d4
 1010 continue
 1000 continue
      return
      end

c============================ REORDER ==============================
c
c Re-order A using Cuthill-McKee
c
c===================================================================
      subroutine reOrder(nnodes,nnz,ia,ja,iau,A,list,tag,
     1                   ianew,iaunew,janew,oldtonew)
      dimension A(4,4,nnz)
      integer ia(nnodes+1),iau(nnodes),ja(nnz)
      integer ianew(nnodes+1),iaunew(nnodes),janew(nnz)
      integer list(nnodes),tag(nnodes)
      integer oldnode,oldcol,newcol
      integer oldtonew(nnodes)
c
      do 1100 i = 1,nnodes
        list(i) = i
        iaunew(i) = iau(i)
        ianew(i)  = ia(i)
 1100 continue
      list(nnodes+1) = nnodes + 1
      ianew(nnodes+1) = ia(nnodes+1)
      do 1101 i = 1,nnz
        janew(i) = ja(i)
 1101 continue
      ireorder = 1 
      if(ireorder.ne.1)return
c
      do 1000 i = 1,nnodes
        tag(i) = 0
 1000 continue
c
c Now go grab all the entries that connect to this node in the number list
c
c     write(6,*)"Entering reOrder"
      icount = 1
      istart = 1
      list(1) = 1
      tag(1)  = 1
 6000 continue
      jstart = ia(list(istart))
      jend   = ia(list(istart)+ 1) - 1
      do 1010 j = jstart,jend
        inode = ja(j)
        if(tag(inode).eq.0)then
          icount = icount + 1
          list(icount) = inode
          tag(inode) = 1
          if(icount.eq.nnodes)goto 7000
        end if
 1010 continue
      istart = istart + 1
      goto 6000
 7000 continue
      list(nnodes+1) = nnodes + 1
c     write(6,*)"Exiting reOrder"
c
c At this point, list(i) is a list of the old node
c numbers stored in the new order
c For example, if list() = {1,3,4,6,5,2}
c then list(1) = old node 1
c      list(2) = old node 3
c      list(3) = old node 4
c      list(4) = old node 6
c      list(5) = old node 5
c      list(6) = old node 2
c
c Now we can construct new ia, iau, and ja
c
      ianew(1) = 1
      do 8000 i = 1,nnodes
       oldnode = list(i)
       nonzeros = ia(oldnode + 1) - ia(oldnode)
       ianew(i+1) = ianew(i) + nonzeros
 8000 continue
c
c Fill oldtonew
c
      do 8030 i = 1,nnodes
        oldtonew(list(i)) = i
 8030 continue
c
c Now get ja
c
      jcount = 0
      do 8010 i = 1,nnodes
        oldnode = list(i)
        jstart = ia(oldnode)
        jend   = ia(oldnode + 1) - 1
        do 8020 j = jstart,jend
          oldcol = ja(j)
          newcol = oldtonew(oldcol)
          jcount = jcount + 1
          janew(jcount) = newcol
 8020   continue
        istart = ianew(i)
        iend   = ianew(i+1) - 1
        call SORTER(istart,iend,janew,iaunew,i)
 8010 continue
c
c Now get iau
c
c     do 8040 i = 1,nnodes
c       jstart = ianew(i)
c       jend   = ianew(i+1) - 1
c       do 8050 j = jstart,jend
c         if(janew(j).eq.i)iaunew(i) = j
c8050   continue
c8040 continue
         iaunew(nnodes + 1) = -1
c
c     write(6,*)"New ia and iau"
c     do 9000 i = 1,nnodes+1
c       write(6,*)ianew(i),iaunew(i)
c9000 continue
c
c     write(6,*)"New ja"
c     do 9001 i = 1,nnodes
c       jstart = ianew(i)
c       jend = ianew(i+1) - 1
c       write(6,*)(janew(j),j=jstart,jend)
c9001 continue
c     stop
c
      return
      end
c===================================================================
c
c Sort each of our bins
c
c===================================================================
      subroutine SORTER(istart,iend,ja,iau,inode)
      integer ja(1),iau(1)
c
      do 1000 i = istart,iend
        min = ja(i)
        minsave = ja(i)
        jsave = i
        do 1010 j = i+1,iend
          if(ja(j).lt.min)then
            min = ja(j)
            jsave = j
          end if
 1010   continue
        ja(i) = min
        ja(jsave) = minsave
        if(ja(i).eq.inode)iau(inode) = i
 1000 continue
c
      return
      end

c============================ TRANSPOSEA ===========================
c
c Transpose A
c
c===================================================================
      subroutine transposeA(nnodes,nnz,ia,ja,iau,A)
      parameter( neqn= 4)
      dimension A(neqn,neqn,nnz)
      dimension rowEntry(neqn,neqn),columnEntry(neqn,neqn)
      integer ia(nnodes + 1),ja(nnz),iau(nnodes)
      integer column
c
      do 1000 i = 1,nnodes
         jstart = iau(i)      ! Start on diagonal
         jend   = ia(i+1) - 1 ! Go to the end of the row
         do 1010 j = jstart,jend
            column = ja(j)
            do 1020 k = 1,neqn
               do 1030 L = 1,neqn
                  rowEntry(k,L) = A(k,L,j)
 1030          continue
 1020       continue
            kstart = ia(column)
            kend   = ia(column+1) - 1
            do 1040 k = kstart,kend
               kindex = k
               if(ja(k).eq.i)goto 1050
 1040       continue
 1050       continue
            do 1060 k = 1,neqn
               do 1070 L = 1,neqn
                  columnEntry(k,L) = A(L,k,kindex)
c                 A(k,L,j) = columnEntry(k,L)
c                 A(k,L,kindex) = rowEntry(L,k)
 1070          continue
 1060       continue
            do 2060 k = 1,neqn
               do 2070 L = 1,neqn
                  A(k,L,j) = columnEntry(k,L)
                  A(k,L,kindex) = rowEntry(L,k)
 2070          continue
 2060       continue
 1010    continue
 1000 continue
      return
      end