module arms2
  use kinddefs
  !use linalg
  use blas
  use armslib
  use ilutp
  use pilu_lib
  implicit none
contains  

  !/*----------------------------------------------------------------------
  ! * Parallel Multi-Level Block ILUT PRECONDITIONER
  ! * arms2    :  parallel arms2
  ! *--------------------------------------------------------------------*/
!int parms_arms_vcsr(parms_Mat self, parms_FactParam param, void *mat, parms_Operator *op)
subroutine arms_factorization(self, param, mat, op)
  !/*---------------------------------------------------------------------
  !| MULTI-LEVEL BLOCK ILUT PRECONDITIONER.
  !| ealier  version:  June 23, 1999  BJS -- 
  !| version2: Dec. 07th, 2000, YS  [reorganized ]
  !| version 3 (latest) Aug. 2005.  [reorganized + includes ddpq]
  !+---------------------------------------------------------------------- 
  !| ON ENTRY:
  !| ========= 
  !| ( Amat ) = original matrix A stored in C-style Compressed Sparse
  !|            Row (CSR) format -- 
  !|            see LIB/heads.h for the formal definition of this format.
  !|
  !| ipar[1:18]  = integer array to store parameters for 
  !|       arms construction (arms2) 
  !|
  !|       ipar[1]:=nlev.  number of levels (reduction processes). 
  !|                       see also "on return" below. 
  !| 
  !|       ipar[2]:= level-reordering option to be used.  
  !|                 if ipar[2]==0 ARMS uses a block independent set ordering
  !|                  with a sort of reverse cutill Mc Kee ordering to build 
  !|                  the blocks. This yields a symmetric ordering. 
  !|                  in this case the reordering routine called is indsetC
  !|                 if ipar[2] == 1, then a nonsymmetric ordering is used.
  !|                  In this case, the B matrix is constructed to be as
  !|                  diagonally dominant as possible and as sparse as possble.
  !|                  in this case the reordering routine called is ddPQ.
  !|                 
  !|       ipar[3]:=bsize. for indset  Dimension of the blocks. 
  !|                  bsize is only a target block size. The size of 
  !|                  each block can vary and is >= bsize. 
  !|                  for ddPQ: this is only the smallest size of the 
  !|                  last level. arms will stop when either the number 
  !|                  of levels reaches nlev (ipar[1]) or the size of the
  !|                  next level (C block) is less than bsize.
  !|
  !|       ipar[4]:=iout   if (iout > 0) statistics on the run are 
  !|                       printed to FILE *ft
  !|
  !|	ipar[5]:= Krylov subspace dimension for last level 
  !|		    ipar[5] == 0 means only backward/forward solve
  !|		    is performed on last level.                  
  !|	ipar[6]:=  maximum # iterations on last level
  !|
  !|       ipar[7-10] NOT used [reserved for later use] - set to zero.
  !| 
  !| The following set method options for arms2. Their default values can
  !| all be set to zero if desired. 
  !|
  !|       ipar[11-14] == meth[1:4] = method flags for interlevel blocks
  !|       ipar[15-18] == meth[1:4] = method flags for last level block - 
  !|       with the following meaning in both cases:
  !|            meth[1] nonsummetric permutations of  1: yes. affects rperm
  !|                    USED FOR LAST SCHUR COMPLEMENT 
  !|            meth[2] permutations of columns 0:no 1: yes. So far this is
  !|                    USED ONLY FOR LAST BLOCK [ILUTP instead of ILUT]. 
  !|                    (so ipar[12] does no matter - enter zero). If 
  !|                    ipar[16] is one then ILUTP will be used instead 
  !|                    of ILUT. Permutation data stored in: perm2. 
  !|            meth[3] diag. row scaling. 0:no 1:yes. Data: D1
  !|            meth[4] diag. column scaling. 0:no 1:yes. Data: D2
  !|       all transformations related to parametres in meth[*] (permutation, 
  !|       scaling,..) are applied to the matrix before processing it 
  !| 
  !| ft       =  file for printing statistics on run
  !|
  !| droptol  = Threshold parameters for dropping elements in ILU 
  !|            factorization.
  !|            droptol[1:5] = related to the multilevel  block factorization
  !|            droptol[6:7] = related to ILU factorization of last block.
  !|            This flexibility is more than is really needed. one can use
  !|            a single parameter for all. it is preferable to use one value
  !|            for droptol[1:5] and another (smaller) for droptol[6:7]
  !|            droptol[1] = threshold for dropping  in L [B]. See piluNEW.c:
  !|            droptol[2] = threshold for dropping  in U [B].
  !|            droptol[3] = threshold for dropping  in L^{-1} F 
  !|            droptol[4] = threshold for dropping  in E U^{-1} 
  !|            droptol[5] = threshold for dropping  in Schur complement
  !|            droptol[6] = threshold for dropping  in L in last block
  !|              [see ilutpC.c]
  !|            droptol[7] = threshold for dropping  in U in last block
  !|              [see ilutpC.c]
  !|             This provides a rich selection - though in practice only 4
  !|             parameters are needed [which can be set to be the same 
  !              actually] -- indeed it makes sense to take
  !|             droptol[1] = droptol[2],  droptol[3] = droptol[4], 
  !|             and droptol[5] = droptol[6]
  !|
  !| lfil     = lfil[1:7] is an array containing the fill-in parameters.
  !|            similar explanations as above, namely:
  !|            lfil[1] = amount of fill-in kept  in L [B]. 
  !|            lfil[2] = amount of fill-in kept  in U [B].
  !|            lfil[3] = amount of fill-in kept  in E L\inv 
  !|            lfil[4] = amount of fill-in kept  in U \inv F
  !|            lfil[5] = amount of fill-in kept  in S    .
  !|            lfil[6] = amount of fill-in kept  in L_S  .
  !|            lfil[7] = amount of fill-in kept  in U_S 
  !|             
  !| tolind   = tolerance parameter used by the indset function. 
  !|            a row is not accepted into the independent set if 
  !|            the *relative* diagonal tolerance is below tolind.
  !|            see indset function for details. Good values are 
  !|            between 0.05 and 0.5 -- larger values tend to be better
  !|            for harder problems.
  !| 
  !| ON RETURN:
  !|=============
  !|
  !| (PreMat)  = arms data structure which consists of two parts:
  !|             levmat and ilsch. 
  !|
  !| ++(levmat)= permuted and sorted matrices for each level of the block 
  !|             factorization stored in PerMat4 struct. Specifically
  !|             each level in levmat contains the 4 matrices in:
  !|
  !|
  !|            |\         |       |
  !|            |  \   U   |       |
  !|            |    \     |   F   |
  !|            |  L   \   |       |
  !|            |        \ |       |
  !|            |----------+-------|
  !|            |          |       |
  !|            |    E     |       |
  !|            |          |       |
  !|            
  !|            plus a few other things. See LIB/heads.h for details.
  !|
  !| ++(ilsch) = IluSpar struct. If the block of the last level is:
  !|
  !|                        |  B    F |
  !|                  A_l = |         | 
  !|                        |  E    C |
  !|
  !|             then IluSpar contains the block C and an ILU
  !|             factorization (matrices L and U) for the last 
  !|             Schur complement [S ~ C - E inv(B) F ]
  !|             (modulo dropping) see LIB/heads.h for details.
  !|
  !| ipar[0]   = number of levels found (may differ from input value) 
  !|
  !+---------------------------------------------------------------------*/
  type(parms_Mat)::self
  type(parms_Operator),pointer:: op
  type(parms_FactParam),target:: param !points to various parameter values
  type(cs_type),intent(in),pointer :: mat
  !/*-------------------- local variables  (initialized)   */
  type(parms_Map):: is

  type(p4_type),pointer :: levc, levmat !all should be null
  type(p4_type),pointer:: levn!,levp !all should be null

  type(cs_type),pointer :: schur, B, F, E, C !all should be null
  type(ilut_type),pointer :: ilsch !should be null 

  real(dp),pointer,dimension(:,:)::dd1,dd2
  real(dp),pointer,dimension(:):: droptol 
  real(dp) ::tolind
  integer :: nlev, bsize, iout, ierr, total_nnz, total_nnz_sch,total_nnz_BLU,total_nnz_E,total_nnz_F
  integer,dimension(:),pointer :: lfil !pointers?
  integer,dimension(4) :: methL, methS
  !/*--------------------  local variables  (not initialized)   */
  integer :: nA, nB, nC, j, jjsub,n, ilev, symperm, schur_start, nbnd, i,iisub,istat,itmp,pos,k
  integer :: nsubrow,nsubcol,nBlast
  !FILE *ft
  !/*--------------------    work arrays:    */
  integer,pointer,dimension(:):: iwork_mem, uwork_mem 
  integer,pointer,dimension(:):: iwork, uwork 

  !/*   timer arrays:  */ 
  !/*   double *symtime, *unstime, *factime, *tottime;*/
  nullify(dd1)
  nullify(dd2)
  nullify(levc)
  nullify(levmat)
  nullify(schur)
  nullify(B)
  nullify(F)
  nullify(E)
  nullify(C)
  nullify(ilsch)
  nullify(droptol)
  nullify(lfil)
  nullify(iwork_mem)
  nullify(uwork_mem)
  nullify(iwork)
  nullify(uwork)
  ierr=0

  !/*---------------------------BEGIN ARMS-------------------------------*/
  !/*   schur matrix starts being original A */ 
  !write(6,*)"arms_fact checkpt 1"  
  !/*-------------------- begin                                         */
  n = param%n	
  nsubrow=mat%nsubrow
  nsubcol=mat%nsubcol					
  nbnd = param%schur_start
  schur_start = param%schur_start
  !is = self%is !parms_Map
  if (schur_start .eq. -1) nbnd = self%is%schur_start
  !ceb this allocation and pointer assignment is fairly confusing
  !  there must be a less muddled way to do this
  !obtain parameters for current A matrix
  !if (param%isalloc.eq.0) then !?
     !call parms_OperatorCreate(newOpt)
     !PARMS_MEMCPY(newOpt%ops, parms_arms_sol_vptr,1)!copies parameters?
     !create new data object
     !allocate(data) 
     call createP4(op%data%levmat)
     levmat => op%data%levmat !should be null here
     call createILUT(op%data%ilus)
     ilsch =>  op%data%ilus !should be null here
     op%data%nsubrow=nsubrow
     op%data%nsubcol=nsubcol
     !write(6,*)"arms_fact checkpt 4"  
     !op => newOpt !pointer to procedure?
  !else 
  !   data => op%data
  !end if
  
!write(6,*)"0: schur_start=",schur_start," nbnd=",nbnd

  !/* compute the number of nonzero entries in mat */
  op%data%nnz_mat = 0
  do i=1, mat%n 
     op%data%nnz_mat = op%data%nnz_mat + mat%nnzrow(i)
  end do
  
  !/* retrieve data from input objects */
  droptol => param%droptol
  tolind  = param%tolind

  nlev = param%ipar(1)!ipar%nlev
!write(6,*) "arms_fact nlev=",nlev
  bsize= param%ipar(3)!ipar%bsize
!write(6,*) "arms_fact bsize=",bsize
  iout = param%ipar(4)!ipar%iout

  !write(6,*) "fil[0,1,2,4,5,6]=",param%lfil(1),param%lfil(2),param%lfil(3),param%lfil(4),param%lfil(5),param%lfil(6),param%lfil(7)
  !write(6,*)"droptol[1,2,3,4,5,6,7]=",droptol(1),droptol(2),droptol(3),droptol(4),droptol(5),droptol(6),droptol(7)
 
  ! do i=1,18
  !   write(6,*)"ipar(",i,")=",param%ipar(i)
  !end do

  !if (iout > 0 ) then
  !   ft = stdout
  !else
  !   ft = NULL
  !end if
  ierr = 0

  !/*---------------------------------------------------------------------
  !| The matrix (a,ja,ia) plays role of Schur compl. from the 0th level.
  !+--------------------------------------------------------------------*/
  nC = mat%n
  nA = mat%n
  n  = mat%n
  if (bsize .ge. n) bsize = n-1 !block size must be less than total amat size?
  levmat%n = n
  levmat%nB = 0!?
  levmat%mylev=1

  schur => mat!pointer to Amat
  levc => op%data%levmat

  !/*--------------------------------------- */ 
  levc%n = 0
  !copy parameter list
  do iisub=1,4
     methL(iisub)=param%ipar(10+iisub)!ipar(11-14):method flags for interlevel blocks
     methS(iisub)=param%ipar(14+iisub)!ipar(15-19):method flags for last level block
     !write(6,*)"methL(",iisub,")=",methL(iisub)
     !write(6,*)"methS(",iisub,")=",methS(iisub)
  end do

  !/*---------------------------------------------------------------------
  !  | The preconditioner construction is divided into two parts:
  !  |   1st part: construct and store multi-level L and U factors;
  !  |   2nd part: construct the ILUT factorization for the coarsest level
  !  +--------------------------------------------------------------------*/
  !  if ( (iout > 0)  && (nlev > 0) ) then !{
  !     fprintf(ft,"  \n");
  !     fprintf(ft,"Level   Total Unknowns    B-block   Coarse set\n");
  !     fprintf(ft,"=====   ==============    =======   ==========\n");
  !     !}
  !  end if
  



  !/*---------------------------------------------------------------------
  !  | main loop to construct multi-level LU preconditioner. Loop is on the
  !  | level ilev. nA is the dimension of matrix A_l at each level.
  !  +--------------------------------------------------------------------*/
  do ilev=1,nlev !

write(6,*)"arms2: ilev=",ilev
     !/*-------------------- new nA is old nC -- from previous level */
     nA = nC
!write(6,*)"arms_fact checkpt 9: nA=",nA," bsize=",bsize  
     if ( nA .le. bsize )  goto 1000  !exit ilev loop

!write(6,*)"arms_fact checkpt 10"  

     !/*-------------------- allocate work space                        */ 
     allocate(iwork_mem(nA),STAT=istat)!allocate work index array
     iwork=>iwork_mem
     !symperm = 0;    !/* 0nly needed in cleanP4 */
     if (param%ipar(2) .eq. 1) then!level_reordering
        !write(6,*)"arms_fact checkpt 10.2"  
        allocate(uwork_mem(nA),STAT=istat)!allocate work element array
        uwork=>uwork_mem
     else
        symperm = 1
        uwork => iwork
     end if
     
!ceb if we scale before pqperm the scaling array is not matched
!ceb maybe we can do row scaling before permutation 

!     !/*-------------------- SCALING*/
!     !ceb scales original input matrix (scale before reordering via pq or indset)
        if (methL(3).eq.1) then !{
!write(6,*)"arms_fact checkpt 11: nA=",nA," calling roscalC(schur)"  
           if(associated(dd1))deallocate(dd1)
           allocate(dd1(nsubrow,nA),STAT=istat)
           call roscalC(schur, dd1,1)
        end if

        if (methL(4).eq.1) then !{
           !write(6,*)"arms_fact checkpt 12: nA=",nA," calling coscalC(schur)"    
           if(associated(dd2))deallocate(dd2)
           allocate(dd2(nsubrow,nA),STAT=istat)
           call coscalC(schur, dd2,1)
           !call coscalC(schur, dd2,0)
        end if

!write(6,*)"arms_fact checkpt 13"  

     !/*--------------------independent-sets-permutation-------------------
     !  |  do reordering -- The matrix and its transpose are used.
     !  +--------------------------------------------------------------------*/
     !ceb get reordering lists for input matrix
     if (param%ipar(2) .eq. 1) then !use nonsymmetric ordering diagonally dominant PQ permutation
!write(6,*),"arms_fact checkpt 14.1 nbnd=",nbnd," tolind",tolind
        call PQperm(schur, bsize, uwork, iwork, nB, tolind, nbnd)
     else !use block independent symmetric ordering ,similar to cuthill-mckee 
        call indsetC(schur, bsize, iwork, nB, tolind, nbnd) !reorder current matrix A=|[B][F]| 
!write(6,*)"arms_fact checkpt 15 nbnd=",nbnd," nB=",nB," bsize=",bsize  
!nB=schur%n
!do i=1,nB
!iwork(i)=i
!uwork(i)=i
!end do 
    end if


     !/*---------------------------------------------------------------------
     ! permutes input matrix according to renumbering arrays iwork and uwork
     !/*---------------------------------------------------------------------
     call rpermC(schur,uwork)!permute schur rows
     call cpermC(schur,iwork,1)!permute schur columns
write(6,*),"arms_fact checkpt 16.0"  

     !/*---------------------------------------------------------------------
     !  | nB is the total number of nodes (rows) in the independent set.
     !  | nC : nA - nB = the size of the reduced system.
     !  +--------------------------------------------------------------------*/
     nC = nA - nB
     nbnd = nbnd- nB

     write(6,*),"arms_fact checkpt 16.1 nB=",nB," nC=",nC," nbnd=",nbnd

!do j=1,schur%n
!   write(6,*)"main: rperm|cperm(",j,")=",uwork(j),iwork(j)
!end do 

     !/*   if the size of B or C is zero , exit the main loop  */
     if ( nB .eq. 0 .or. nC .eq. 0 )  then
        write(6,*) "Exiting reduction loop nB=",nB," nC=",nC
          goto 1000 !exit ilev loop
     end if

!write(6,*),"arms_fact checkpt 16.2"  
     !/*---------------------------------------------------------------------
     !  | The matrix for the current level is in (schur).
     !  | The permutations arrays are in iwork and uwork (row).
     !  | The routines rpermC, cpermC permute the matrix in place.
     !  *-----------------------------------------------------------------------*/
     !/*   DEBUG : SHOULD THIS GO BEFORE GOTO LABEL1000 ?? */

!  do i=1,schur%n
!     do j=1,schur%nnzrow(i)
!        write(20,*) i,schur%pj(i)%cols(j)
!        write(21,*) i,schur%pj(i)%cols(j)
!  !      !write(21,*) iwork(i),iwork(schur%pj(i)%cols(j))
!  !      write(22,*) ((schur%pa(i)%submat(iisub,jjsub,j),iisub=1,schur%nsubrow),jjsub=1,schur%nsubcol)
!     end do
!  end do

!do j=1,schur%n
!   write(6,*)"main: rperm|cperm(",j,")=",uwork(j),iwork(j)
!end do

!     !/*---------------------------------------------------------------------
!     ! permutes input matrix according to renumbering arrays iwork and uwork
!     !/*---------------------------------------------------------------------
     !call rpermC(schur,uwork)!permute schur rows
     !call cpermC(schur,iwork,1)!permute schur columns



!  do i=1,schur%n
!     do j=1,schur%nnzrow(i)
!!  !      write(94,*) i, schur%pj(i)%cols(j)
!        write(879,*) i, schur%pj(i)%cols(j),((schur%pa(i)%submat(iisub,jjsub,j),iisub=1,schur%nsubrow),jjsub=1,schur%nsubcol)
!     end do
!  end do


     !/*-----------------------------------------------------------------------
     !  | If this is the first level, the permuted matrix is stored in 
     !  | (levc) = (levmat).  Otherwise, the next level is created in (levc).
     !  +--------------------------------------------------------------------*/
     if (ilev > 1) then
write(6,*)"arms_fact checkpt 22.0 ilev:",ilev  
        ! /*-   delete C matrix of any level except last one (no longer needed) */
        call cleanCS(C) 
        !pointer shuffle
        nullify(levn)
        call createP4(levn)!ceb need to create levn?
        !next levc%prev is current levc
        levn%prev=>levc
        !current levc%next=>levn
        levc%next => levn !levn is new target
        !new levc=>levn (i.e. old levc%next)  
        levc => levn      !levn is target
     end if

     !/*-------------------- p4ptr struct for current schur complement */
!write(6,*)"arms_fact checkpt 21.2 nB=",nB," nC=",nC
     call csSplit4(schur, nB, nC, B, F, E, C)!matrices BFEC created
     call setupP4(levc, nB, nC, F, E, ilev)

!write(6,*)"arms_fact checkpt 21.2 levc%nB=",levc%nB
 !write(6,*)"mat levc:"
!call print_dataP4(levc)

     !/*--------------------     copy a few pointers       ---- */   
     !Attach permutation arrays to levc matrix   
     levc%perm  => iwork
     levc%rperm => uwork 

     !ceb we want to keep the memory that iwork,uwork point to so we will nullify these 
     !pointers. Not sure if allocating uwork_mem again will overwrite that space or 
     !open up a new memory location.
     nullify(iwork)
     nullify(uwork)
     nullify(iwork_mem)
     nullify(uwork_mem)

     levc%symperm = symperm
     levc%D1=>dd1
     nullify(dd1)!keep memory by taking this pointer away

     !do i=1,levc%n
     !   do iisub=1,levc%nsubrow
     !      write(6,*)"arms2: levc%D1(",iisub,i,")=",levc%D1(iisub,i)
     !   end do
     !end do

     levc%D2=>dd2 
     nullify(dd2)!keep memory by taking this pointer away

!write(6,*)"arms_fact checkpt 24.6 levc%n=",levc%n  
!     do j=1,levc%n
!        write(6,*)"arms_fact checkpt 24  ilev:",ilev," arms: levc%rperm(",j,")=",levc%rperm(j)
!    end do


     !/*---------------------------------------------------------------------
     !  | a copy of the matrix (schur) has been permuted. Now perform the 
     !  | block factorization: 
     !  |
     !  | | B   F |       | L       0 |     | U  L^-1 F |
     !  | |       |   =   |           |  X  |           | = L x U
     !  | | E   C |       | E U^-1  I |     | 0    A1   |
     !  |   
     !  | The factors E U^-1 and L^-1 F are discarded after the factorization.
     !  |
     !  +--------------------------------------------------------------------*/ 
     !if (iout > 0) fprintf(ft,"%3d %13d %13d %10d\n", ilev+1,nA,nB,nC);
     !/*---------------------------------------------------------------------
     !  | PILUT constructs one level of the block ILU fact.  The permuted matrix
     !  | is in (levc).  The L and U factors will be stored in the p4mat struct.
     !  | destroy current Schur  complement - no longer needed  - and set-up new
     !  | one for next level...
     !  +--------------------------------------------------------------------*/

     nullify(schur) !nullify rather than deallocate so we dont destroy original mat

     itmp=1
     call setupCS(schur,nC,nsubrow,nsubcol,itmp)!create new schur
!write(6,*)"arms_fact checkpt 27 calling pilu"  
     !/*----------------------------------------------------------------------
     !  | calling PILU to construct this level block factorization
     !  | ! core dump in extreme case of empty matrices.
     !  +----------------------------------------------------------------------*/
     call pilu(levc, B, C, droptol, param%lfil, schur)!factorize levc into B=LU->levc, C, F,E,schur

!write(6,*)"arms_fact checkpt 28: schur%n=",schur%n  
     call cleanCS(B)!delete B, we have it in LU form in levc

  end do!end loop over ilev
  
  !/*---------------------------------------------------------------------
  !  |   done with the reduction. Record the number of levels in ipar[0] 
  !  +--------------------------------------------------------------------*/
1000 continue







  !/* printf (" nnz_Schur %d \n",cs_nnz (schur)); */
  !write(6,*)"arms_fact checkpt 30:   ilev=",ilev 
  if (ilev.ge.nlev)ilev=nlev
  write(6,*)"arms_fact checkpt 30.1: ilev=",ilev 

  nullify(levc%next)
  param%ipar(1) = ilev
  op%data%nlev = ilev
  op%data%n = n
  nC = schur%n





write(6,*)"arms_fact checkpt 31 nC=schur%n=",nC  

     call setupILUT(ilsch,nC,nsubrow,nsubcol)!create matrix ilsch
  

   
!write(6,*)"arms_fact checkpt 32"  
     
     !/*--------------------------------------------------------------------*/
     !/* define C-matrix (member of ilsch) to be last C matrix */ 
     !if (ilev > 0) ilsch%C=>C 
     if (ilev > 1) ilsch%C=>C 
  
     !ceb skip scaling and permutations for case that nC==nA
     if (nC.ne.nA) then
        !/*-------------------- for ilut fact of schur matrix */
        if(nC < n) then
           !/*  SCALING  */
           nullify(ilsch%D1)! = NULL;
           if (methS(3).eq.1) then
              !           write(6,*)"arms_fact checkpt 33.1 nC=",nC," calling roscalC(schur)"  
              allocate(ilsch%D1(ilsch%nsubrow,nC),STAT=istat)
              call roscalC(schur, ilsch%D1, 1) 
           end if
           !  
           nullify(ilsch%D2)!  = NULL;
           if (methS(4).eq.1) then
              !           write(6,*)"arms_fact checkpt 33.3 nC=",nC," calling coscalC(schur)"  
              allocate(ilsch%D2(ilsch%nsubrow,nC),STAT=istat)
              call coscalC(schur, ilsch%D2, 1)
              !call coscalC(schur, ilsch%D2, 0)!inf norm
           end if
           
        else !schur is complete matrix which may have been scaled already
           ilsch%D1=>dd1
           ilsch%D2=>dd2
        end if

        
        !/*---------------------------------------------------------------------
        !  |     get ILUT factorization for the last reduced system.
        !  +--------------------------------------------------------------------*/
        if(associated(iwork_mem)) deallocate(iwork_mem)! = NULL;
        nullify(iwork)
        if(associated(uwork_mem)) deallocate(uwork_mem)! = NULL;
        nullify(uwork)
        
        if (methS(1).eq.1) then !pq permute schur complement
           !write(6,*) "arms2 checkpt 34.9"
           allocate(iwork_mem(nC),STAT=istat)
           allocate(uwork_mem(nC),STAT=istat)
           !if(.not.associated(uwork,iwork)) then !if uwork does not point to iwork
           !   allocate(uwork(nC),STAT=istat)
           !end if
           iwork=>iwork_mem
           uwork=>uwork_mem
           
           tolind = 0.0 
           !write(6,*) "arms2 checkpt 35.1, nB=",nb," nbnd=",nbnd," tolind=",tolind 
           call PQperm(schur, bsize, uwork, iwork, nB, tolind, nbnd)!get permutations for schur matrix
           
           !do j=1,schur%n
           !   write(6,*)"schur comp: rperm|cperm(",j,")=",uwork(j),iwork(j)
           !end do
           !      write(6,*) "arms2 checkpt 35.2"
           call rpermC(schur,uwork)!permute schur rows
           call cpermC(schur,iwork,1)!permute schur columns
           !  do i=1,schur%n
           !     do j=1,schur%nnzrow(i)
           !        write(780,*) i, schur%pj(i)%cols(j),((schur%pa(i)%submat(iisub,jjsub,j),iisub=1,schur%nsubrow),jjsub=1,schur%nsubcol)
           !     end do
           !  end do
        end if
     else !C=> last B
        ilsch%D1=>dd1
        ilsch%D2=>dd2
     end if

write(6,*) "arms2 checkpt 35.9"
     ilsch%rperm => uwork!should this be a copy or a pointer 
     ilsch%perm  => iwork!should this be a copy or a pointer 

     nullify(uwork_mem)     
     nullify(iwork_mem) 

     nullify(ilsch%perm2)!for column pivoting in ilutp
     !=====================================================
     if (methS(1) .eq. 1) then !call ilut that doesn't do partial pivoting
        !write(6,*)"arms_fact checkpt 36: calling ilutD "  
        !call ilutD(schur, droptol, param%lfil, ilsch)!factor schur into LU form and store in ilsch       

        if(.not.associated(ilsch%perm2)) allocate(ilsch%perm2(nC),STAT=istat)
        do j=1,nC
           ilsch%perm2(j) = j
        end do
        !!write(6,*)"arms_fact checkpt 37: calling ilutpC nC=",nC  
        !!!factors final schur complement into LU threshold form and store in ilsch,
        !!! allows independent column reordering for partial pivoting
        call ilutpC(schur, droptol, param%lfil, PERMTOL, nC, ilsch) 
     else!=================================================
        if(.not.associated(ilsch%perm2)) allocate(ilsch%perm2(nC),STAT=istat)
        do j=1,nC
           ilsch%perm2(j) = j
        end do
        write(6,*)"arms_fact checkpt 37: calling ilutpC nC=",nC  
        !factors final schur complement into LU threshold form and store in ilsch,
        ! allows independent column reordering for partial pivoting
        call ilutpC(schur, droptol, param%lfil, PERMTOL, nC, ilsch) 
     end if!===============================================
     




     
  !/*---------- OPTIMIZATION: NEED TO COMPOUND THE TWO
  !  RIGHT PERMUTATIONS -- CHANGES HERE AND IN 
  !  USCHUR SOLVE ==  compound permutations */     
  
  if (ierr.eq.1) then
     !fprintf(ft," ERROR IN  ILUT -- IERR = %d\n", ierr); 
     return!(1); 
  end if
    
  !/* Last Schur complement no longer needed */
  ! info now stored in ilusch, we can now eliminate A_lev
  call cleanCS(schur)!ceb what if we need for GMRes iterations? (V or W cycle)

  !data%nnz_prec = nnz_arms(data)!ceb  not sure where this is defined
  op%data%ind = n - ilsch%n
  if (ilev.eq.1) then
     op%data%schur_start = (n - ilsch%n)+1
!write(6,*)"1: schur_start=",op%data%schur_start
  else
     !is = self%is !parms_Map
     op%data%schur_start = self%is%schur_start !
!write(6,*)"2 :schur_start=",self%is%schur_start
  end if  
  !copy input parameters to factorization matrix objects
  do iisub=1,18
     op%data%ipar(iisub)=param%ipar(iisub)
  enddo
  op%data%pgfpar(1)=param%pgfpar(1)
  op%data%pgfpar(2)=param%pgfpar(2)

  
  !ceb output matrix structure to file


!write(6,*)"arms_fact checkpt 43 ilev=",ilev
!  write(6,*)"levc:"
!  call print_dataP4(levc)
!  write(6,*)"levc%prev:"
! if(associated(levc%prev)) call print_dataP4(levc%prev) 
!write(6,*)"arms_fact checkpt 43.5:"


 
  !descend to lowest level matrix
  do
     !write(6,*)"43.5 ilev=",ilev
     if (.not.associated(levc%prev)) exit
     levc=>levc%prev
     ilev=ilev-1
  end do
 
 
!write(6,*)"arms2 43.6 ilev=",ilev  

  !write LU's from level data
  !call print_dataP4(levc)

  nBlast=0
  total_nnz_BLU=0
  total_nnz_F=0
  total_nnz_E=0

  do
     write(6,*)"arms2 checkpt 45: ilev=",ilev
     call print_dataP4(levc)
     !     write(6,*)"arms2 checkpt 44: levc%nB=",levc%nB
     
     !write L,U cols for each row in levc
     do i=1,levc%nB
        if(associated(levc%L) .and. allocated(levc%L%nnzrow)) then
           total_nnz_BLU=total_nnz_BLU+levc%L%nnzrow(i)+levc%U%nnzrow(i)
           do j=1,levc%L%nnzrow(i)
              write(55,*) i+nBlast,levc%L%pj(i)%cols(j)+nBlast 
           end do
           do k=1,levc%U%nnzrow(i)
              write(55,*) i+nBlast,levc%U%pj(i)%cols(k)+nBlast
           end do
        end if
     end do

     !write(6,*)"arms2 checkpt 44.5"
     !do j=1,levc%n
     !   write(6,*)"lev:",ilev," arms: levc%rperm(",j,")=",levc%rperm(j)
     !end do

     do i=1,levc%nB
        if(allocated(levc%F%nnzrow)) then
           total_nnz_F=total_nnz_F+levc%F%nnzrow(i)
           do j=1,levc%F%nnzrow(i)
              write(56,*) i+nBlast,         (n-(levc%n-levc%nB))+levc%F%pj(i)%cols(j)
           end do
        end if
     end do

     do i=1,levc%n-levc%nB
        if(allocated(levc%E%nnzrow)) then
           total_nnz_E=total_nnz_E+levc%E%nnzrow(i)
           do j=1,levc%E%nnzrow(i)
             ! write(57,*) i+nBlast+levc%nB,         (n-(levc%nB))+levc%E%pj(i)%cols(j)
              !write(57,*) i+nBlast+levc%nB,         (n-(levc%nB)-(levc%n-levc%nB))+levc%E%pj(i)%cols(j)
              write(57,*) i+nBlast+levc%nB,         (n-levc%n)+levc%E%pj(i)%cols(j)
           end do
        end if
     end do

     nBlast=nBlast+levc%nB



     if (.not.associated(levc%next)) exit
     levc=>levc%next
     ilev=ilev+1
  end do

  
  write(6,*)" total LU nnz=",total_nnz_BLU," nBlast:",nBlast
!     call print_dataP4(levc)

!write(6,*)"arms2 checkpt 45"

  
  ! write schur
  total_nnz_sch=0
  do i=1,ilsch%n
!write(6,*) "ilsch:",i
     total_nnz_sch=total_nnz_sch+ilsch%L%nnzrow(i)+ilsch%U%nnzrow(i)
     !if(ilsch%L%nnzrow(i)>0) then
        do j=1,ilsch%L%nnzrow(i)
!write(6,*) "L:",i,ilsch%L%nnzrow(i), j
           write(58,*) i+nBlast,ilsch%L%pj(i)%cols(j)+nBlast
        end do
     !end if
     do k=1,ilsch%U%nnzrow(i)
 !write(6,*) "U:",i,ilsch%U%nnzrow(i), k
       write(58,*) i+nBlast,ilsch%U%pj(i)%cols(k)+nBlast
     end do
  end do

 
  total_nnz=total_nnz_BLU+total_nnz_F+total_nnz_E+total_nnz_sch

write(6,*)" total LU nnz=",total_nnz_BLU
write(6,*)" total F nnz=",total_nnz_F
write(6,*)" total E nnz=",total_nnz_E
write(6,*)" total schur nnz=",total_nnz_sch
write(6,*)" total LU+last schur nnz=",total_nnz
  

   nullify(levc)
   nullify(ilsch)
   nullify(iwork)
   nullify(uwork)
   nullify(iwork_mem)
   nullify(uwork_mem)
   nullify(schur)
   nullify(levmat)
   if(associated(B)) nullify(B)
   nullify(F)
   nullify(E)
   nullify(C)
   nullify(droptol)
   nullify(lfil)
   nullify(dd1)
   nullify(dd2)
!   write(6,*)"arms_fact checkpt 99"  

   !values returned in op%data%ilus
   !               and op%data%levmat

  return! 0;
  !}!/*-----end-of-ARMS2----------------------------------------------------
  !  +--------------------------------------------------------------------*/
end subroutine arms_factorization
!|--------------------------------------------------------------------*/




















end module arms2



