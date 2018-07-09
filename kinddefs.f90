module kinddefs
  implicit none
  private

  public :: dp, jp, odp, dqp 
  public :: intarray,realarray

  public :: ja_type, ja_ptr, submat_type, submat_ptr 
  public :: matrix_type, matrixptr, cs_type, csptr
  public :: vbs_type, vbsptr, vbilu_type, vbiluptr
  public :: ilu_type, iluptr, ilut_type, ilutptr, p4_type, p4ptr
  public :: arms_data_type, arms_dataptr
  public :: SMat_type, SMatptr, SPre_type, SPreptr
  public :: parms_Mat,parms_Map,parms_FactParam,parms_Operator
  public :: parms_Operator_ops
  public :: DBL_EPSILON,PERMTOL,ofst
  public :: solver_params_type

  public :: system_r4, system_r8, system_r16, r4, r8
  public :: system_i1, system_i2, system_i4, system_i8
! single precision (IEEE 754)
  integer, parameter :: system_r4=selected_real_kind(6, 37)
! double precision (IEEE 754)
  integer, parameter :: system_r8=selected_real_kind(15, 307)
! quad precision (IEEE 754?)
  integer, parameter :: system_r16=selected_real_kind(33, 4931)
! one byte integer (greater than 10e2)
  integer, parameter :: system_i1=selected_int_kind(2)
! two byte integer (greater than 10e3)
  integer, parameter :: system_i2=selected_int_kind(3)
! four byte integer (greater than 10e5)
  integer, parameter :: system_i4=selected_int_kind(5)
! eight byte integer (greater than 10e10)
  integer, parameter :: system_i8=selected_int_kind(10)
  integer, parameter :: dp =system_r8 ! default precision kind
! Types controlling precision of linear solves.
  integer, parameter :: jp =system_r8 ! Diagonal LU Jacobian precision kind
  integer, parameter :: dqp=system_r4 ! Delta-q precision kind
  integer, parameter :: odp=system_r4 ! Off-diagonal Jacobian precision kind
  integer, parameter :: r4=system_r4 ! shorthand for system_r4
  integer, parameter :: r8=system_r8 ! shorthand for system_r8

  real(dp), PARAMETER :: DBL_EPSILON = 2.2204460492503131e-16 !// double epsilon
  real(dp), PARAMETER :: PERMTOL = 0.99 !  0 --> no permutation 0.01 to 0.1 good

  integer, parameter :: ofst = 1 !offset from c style index to fortran
! types facilitating matrix definitions
  !type :: nnz_type
  !end type nnz_type

  !type :: rownnz_type
  !   integer  :: nnz ! number of nonzero block cols in row
  !end type rownnz_type

  type :: solver_params_type
     character(len=130)::matrix_file
     integer :: mformat,method,precon,nsubiters,krylovdirs
     integer :: restarts,block_size,nlev,mreord
     integer :: fill_BL,fill_BU,fill_EU,fill_LF,fill_S,fill_SL,fill_SU
     integer :: rscale_B,cscale_B,rscale_S,cscale_S,pq_S
     real(dp) :: convgtol,droptol
  end type solver_params_type

  type :: intarray
     integer,dimension(:),allocatable :: array
  end type intarray

  type :: realarray
     real(dp),dimension(:),allocatable :: array
  end type realarray

  type :: ja_type
     !integer, target,dimension(:), allocatable :: cols ! list of nonzero block cols in row
     integer, dimension(:), allocatable :: cols ! list of nonzero block cols in row
  end type ja_type
  !|--------------------------------------------------------------------*/
  type :: ja_ptr
     type(ja_type),pointer::ptr
  end type ja_ptr
  !|--------------------------------------------------------------------*/

  type :: submat_type
     !real(dp),target,dimension(:,:,:), allocatable :: submat
     real(dp),dimension(:,:,:), allocatable :: submat
  end type submat_type
  !|--------------------------------------------------------------------*/
  type::submat_ptr
     type(submat_type),pointer::ptr
  end type submat_ptr
  !|--------------------------------------------------------------------*/

  type :: matrix_type
     integer :: nnodes
     integer :: nsubrow
     integer :: nsubcol
     integer,          dimension(:), allocatable :: rownnz
     integer,dimension(:),allocatable :: iau  !/* index (pj) of diagonal in each row */
     type(ja_type),    dimension(:), allocatable :: ja
     type(submat_type),dimension(:), allocatable :: row
  end type matrix_type
  !|--------------------------------------------------------------------*/
  type :: matrixptr
     type(matrix_type),pointer::ptr
  end type matrixptr
  !|--------------------------------------------------------------------*/

  ! Description: struct parms_vcsr.
  type :: cs_type 
     integer :: n        ! the dimension of the matrix (number of rows) 
     integer :: nsubrow
     integer :: nsubcol
     integer,dimension(:),allocatable :: nnzrow  !/* length of each row   */
     integer,dimension(:),allocatable :: pd  !/* index (pj) of diagonal in each row */
     !integer,dimension(:),pointer :: nnzrow  !/* length of each row   */
     !integer,dimension(:),allocatable :: space !//!<  length of space ( a work array)       
     integer,dimension(:,:),allocatable :: space !//!<  length of space ( a work array)      
     type(ja_type),    dimension(:), allocatable :: pj ! Parameter: pj = An indirect pointer to store column indices.
     !type(ja_type)    ,pointer,dimension(:) :: pj ! Parameter: pj = An indirect pointer to store column indices.
     !type(ja_ptr)    ,pointer,dimension(:) :: pj ! Parameter: pj = An indirect pointer to store column indices.
     type(submat_type),dimension(:), allocatable :: pa ! parameter: pa = An indirect pointer to store corresponding nonzero entries. 
     !type(submat_type),pointer,dimension(:), allocatable  :: pa ! parameter: pa = An indirect pointer to store corresponding nonzero entries. 
     !type(submat_ptr),dimension(:) :: pa ! parameter: pa = An indirect pointer to store corresponding nonzero entries. 

  end type cs_type
  !|--------------------------------------------------------------------*/
  type :: csptr
     type(cs_type),pointer ::ptr
  end type csptr
  !|--------------------------------------------------------------------*/

  !typedef struct parms_vcsr SparMat;
  !typedef parms_vcsr csptr;
  !typedef double *BData;
  
  type :: vbs_type
     integer :: n       !/* the block row dimension of the matrix      */
     integer :: nsubrow
     integer :: nsubcol
     integer,dimension(:),allocatable :: bsz !the row/col of the first element of each diagonal block
     integer,dimension(:),allocatable :: nnzrow  ! length of each row
     !type(ja_type),    dimension(:), allocatable :: ja !/* pointer-to-pointer to store column indices */
     type(ja_type),    pointer,dimension(:) :: ja !/* pointer-to-pointer to store column indices */
     !type(submat_type),dimension(:), allocatable :: ba !/* pointer-to-pointer to store nonzero blocks */
     type(submat_type),pointer,dimension(:) :: ba !/* pointer-to-pointer to store nonzero blocks */
     !type(submat_type),dimension(:), allocatable :: D!/* to store inversion of diagonals  */
     !real(dp),dimension(:,:), allocatable :: D!/* to store inversion of diagonals  */
     real(dp),pointer,dimension(:,:) :: D!/* to store inversion of diagonals  */
  end type vbs_type
  !|--------------------------------------------------------------------*/
  type :: vbsptr
     type(vbs_type),pointer::ptr
  end type vbsptr
  !|--------------------------------------------------------------------*/

  type :: vbilu_type
     integer n
     integer :: nsubrow
     integer :: nsubcol
     integer,dimension(:),allocatable :: bsz !the row/col of the first element of each diagonal block
     !type(submat_type),dimension(:), allocatable :: D !/* diagonal blocks  */
     real(dp),pointer,dimension(:) :: D !/* diagonal blocks  */
     type(vbs_type),pointer :: L     !/* L part blocks                              */
     type(vbs_type),pointer :: U     !/* U part blocks                              */
     !integer,dimension(:,:),allocatable :: work ! working buffer 
     integer,pointer,dimension(:,:) :: work ! working buffer 
     !BData bf     !/* buffer of a temp block    */
     integer DiagOpt  !/* Option for diagonal inversion/solutiob     */
     !/* opt =  1 -->> call luinv                   */
     !/* opt == 2 -->> block inverted call dgemv    */
     !} VBILUSpar, *vbiluptr; 
  end type vbilu_type
  !|--------------------------------------------------------------------*/
  type :: vbiluptr
     type(vbilu_type),pointer::ptr
  end type vbiluptr
  !|--------------------------------------------------------------------*/

  type :: ilu_type ! ILUfac {
     integer n
     integer :: nsubrow
     integer :: nsubcol
     type(cs_type),pointer :: L ! L part elements
     type(cs_type),pointer :: U      ! U part elements
     !type(submat_type),dimension(:), allocatable :: D ! diagonal elements
     !real(dp),dimension(:,:), allocatable :: D ! diagonal elements
     real(dp),pointer,dimension(:,:) :: D ! diagonal elements
     !integer,dimension(:,:),allocatable :: work ! working buffer 
     integer,pointer,dimension(:,:) :: work ! working buffer 
     !} ILUSpar, LDUmat, *iluptr;
  end type ilu_type
  !|--------------------------------------------------------------------*/
  type :: iluptr
     type(ilu_type),pointer::ptr
  end type iluptr
  !|--------------------------------------------------------------------*/



  !/* -------------------------------------------------------------------*/
  !typedef struct ILUTfac *ilutptr;
  type :: ilut_type !ILUTfac !{
     !/*------------------------------------------------------------
     !| struct for storing data related to the last schur complement 
     !| we need to store the C matrix associated with the last block
     !| and the ILUT factorization of the related Schur complement.
     !| 
     !| n       = size of C block = size of Schur complement
     !| C       = C block of last level matrix. 
     !| L, U    = ILU factors of last schur complement. 
     !|
     !| meth[4] = parameters for defining variants in factorization 
     !|           - see function readin for details
     !| rperm    = row permutation used for very nonsymmetric matrices 
     !|            [such as bottleneck transversal] -- NOT IN THIS VERSION
     !| perm2     = unsymmetric permutation (columns) - used primarily
     !|           for the ILUTP version of ILUT/.. 
     !| D1, D2  = diagonal matrices (left, right) used for scaling
     !|           if scaling option is turned on. Note that the 
     !|           method works by scaling the whole matrix first
     !|           (at any level) before anything else is done. 
     !| wk     = a work vector of length n needed for various tasks
     !|            [reduces number of calls to malloc]           
     !|-----------------------------------------------------------*/
     integer :: n                  
     integer :: nsubrow
     integer :: nsubcol
     ! /*-------------------- C matrix of last block */
     type(cs_type),pointer:: C
     !  /* LU factorization       */
     type(cs_type),pointer:: L
     type(cs_type),pointer:: U
     ! /*--------------------  related to variants and methods */
     ! /*    int meth[4];   */
     !integer,pointer,dimension(:),allocatable :: rperm ! row single-sinded permutation
     integer,pointer,dimension(:) :: rperm ! row single-sinded permutation
     !integer,pointer,dimension(:),allocatable :: perm  ! col permutation 
     integer,pointer,dimension(:) :: perm  ! col permutation 
     !integer,pointer,dimension(:),allocatable :: perm2 ! column permutation coming from pivoting in ILU
     integer,pointer,dimension(:) :: perm2 ! column permutation coming from pivoting in ILU
     !type(submat_type),dimension(:), allocatable :: D1 !/* diagonal scaling row    */ 
     real(dp),pointer,dimension(:,:) :: D1 !/* diagonal scaling row    */ 
     !type(submat_type),dimension(:), allocatable :: D2 !/* diagonal scaling columns*/
     real(dp),pointer,dimension(:,:) :: D2 !/* diagonal scaling columns*/
     !type(submat_type),dimension(:), allocatable :: wk !/* work array*/
     !type(submat_type),pointer,dimension(:,:) :: wk !/* work array*/
     real(dp),pointer,dimension(:,:) :: wk !/* work array*/
     !} IluSpar;
  end type ilut_type
  !|--------------------------------------------------------------------*/
  type :: ilutptr
     type(ilut_type),pointer::ptr
  end type ilutptr
  !|--------------------------------------------------------------------*/





  !typedef struct PerMat4 *p4ptr;
  type :: p4_type !{
     !/*------------------------------------------------------------
     !| struct for storing the block LU factorization 
     !| contains all the block factors except the 
     !| data related to the last block. 
     !| n       = size of current block
     !| symperm = whether or not permutations are symmetric.
     !|           used only in cleanP4..
     !| nB      = size of B-block
     !| L, U    = ILU factors of B-block
     !| F, E    = sparse matrices in (1,2) and (2,1) 
     !|           parts of matrix. 
     !| perm    = (symmetric) permutation used in factorization
     !|           comes from the independent set ordering
     !| rperm   = unsymmetric permutation (rows) not used in this
     !|           version -- but left here for compatibility..
     !| D1, D2  = diagonal matrices (left, right) used for scaling
     !|           if scaling option is turned on. Note that the 
     !|           method works by scaling the whole matrix first
     !|           (at any level) before anything else is done. 
     !| wk     = a work vector of length n needed for various tasks
     !|            [reduces number of calls to malloc]           
     !|----------------------------------------------------------*/ 
     integer :: n      
     integer :: nsubrow
     integer :: nsubcol            
     integer :: nB 
     integer :: symperm
     integer :: mylev
     !/*   LU factors  */
     !type(csptr) :: L
     type(cs_type),pointer :: L
     !type(csptr) :: U
     type(cs_type),pointer :: U
     !/* E, F blocks   */
     !type(csptr) :: E
     type(cs_type),pointer :: E
     !type(csptr) :: F
     type(cs_type),pointer :: F
     !integer,dimension(:),allocatable :: rperm ! /* row permutation         */ 
     integer,pointer,dimension(:) :: rperm ! /* row permutation         */ 
     !integer,dimension(:),allocatable :: perm ! /* col permutation         */ 
     integer,pointer,dimension(:) :: perm ! /* col permutation         */ 
     !real(dp),dimension(:,:), allocatable :: D1 !/* diagonal scaling row    */ 
     real(dp),pointer,dimension(:,:) :: D1 !/* diagonal scaling row    */ 
     !real(dp),dimension(:,:), allocatable :: D2 !/* diagonal scaling columns*/
     real(dp),pointer,dimension(:,:) :: D2 !/* diagonal scaling columns*/
     !type(submat_type),dimension(:), allocatable :: wk !/* work array*/
     !type(submat_type),pointer,dimension(:,:) :: wk !/* work array*/
     real(dp),pointer,dimension(:,:) :: wk !/* work array*/
     !/* pointer to next and previous struct         */
     type(p4_type),pointer :: prev
     type(p4_type),pointer :: next
  end type p4_type
!/*-----------------------------------------------------------------------
  type :: p4ptr
     type(p4_type),pointer::ptr
  end type p4ptr
!/*-----------------------------------------------------------------------





  !#define MAX_BLOCK_SIZE   100
  !!/* FORTRAN style vblock format, compatible for many FORTRAN routines */
  !#define DATA(a,row,i,j)  (a[(j)*(row)+(i)])
  !!/* the dimension of ith Block */
  !#define B_DIM(bs,i)      (bs[i+1]-bs[i])
  !typedef struct parms_vcsr SparMat;
  !typedef parms_vcsr csptr;
  !typedef double *BData;


  type :: arms_data_type !typedef struct parms_arms_data {
     integer :: n !/* dimension of matrix */
     integer :: nsubrow
     integer :: nsubcol
     integer :: nlev !/* number of levels */
     integer :: schur_start
     integer :: ind
     integer :: nnz_mat
     integer :: nnz_prec
     type(ilut_type),pointer :: ilus
     type(p4_type),pointer :: levmat
     integer,dimension(18) :: ipar
     real(dp),dimension(2):: pgfpar
  end type arms_data_type
!/*-----------------------------------------------------------------------
  type :: arms_dataptr
     type(arms_data_type),pointer::ptr
  end type arms_dataptr
!/*-----------------------------------------------------------------------

  !typedef struct __CompressType
  !{
  !  int grp;   !/* -1: begin new group, >=0: belong to grp-th row */
  !  int count; !/* block size, valid only if grp = -1 */
  !} CompressType;
  


  type :: SMat_type !typedef struct _SMat {
     !  /*-------------------- 3 types of matrices so far */
     integer :: n 
     integer :: nsubrow
     integer :: nsubcol
     integer :: Mtype           !/*--  type 1 = CSR, 2 = VBCSR, 3 = LDU    */
     type(cs_type),pointer :: CSR           !/* place holder for a CSR/CSC type matrix */
     type(ilu_type),pointer :: LDU          !/* struct for an LDU type matrix          */
     type(vbs_type),pointer :: VBCSR        !/* place holder for a block matrix        */
     !  void (*matvec)(struct _SMat*, double *, double *);
     !} SMat, *SMatptr;
  end type SMat_type
  !|--------------------------------------------------------------------*/
  type :: SMatptr
     type(SMat_type),pointer::ptr
  end type SMatptr
  !|--------------------------------------------------------------------*/

  type :: SPre_type !typedef struct _SPre {
     !  !/*-------------------- 3 types of matrices so far */
     integer :: Ptype           !/*-- Ptype 1 = ILU, 2 = VBILU, 3 = Crout */
     type(ilu_type),pointer :: ILU        !/* struct for an ILU type preconditioner */
     type(vbilu_type),pointer :: VBILU      !/* struct for a block     preconditioner */
     type(arms_data_type),pointer :: ARMS           !/* struct for a block     preconditioner */
     !  int (*precon) (double *, double *, struct _SPre*); 
     !} SPre, *SPreptr;
  end type SPre_type
  !|--------------------------------------------------------------------*/
  type :: SPreptr
     type(SPre_type),pointer::ptr
  end type Spreptr
  !|--------------------------------------------------------------------*/








!/*!
!  \file   parms_mat_impl.h
!  \brief  The matrix structure and functions
!  \author zzli
!  \date   2006-05-08
!*/

!#ifndef _PARMS_MAT_IMPL_H_
!#define _PARMS_MAT_IMPL_H_

!#include "parms_mat.h"
!#include "parms_operator.h"
!#include "parms_map_impl.h"

!#define PARMS_MAT_FORMAT_KIND    0x00f00000
!#define PARMS_MAT_SHIFT          20
!#define PARMS_MAT_GET_KIND(a)       (((a)&PARMS_MAT_FORMAT_KIND) >> PARMS_MAT_SHIFT)
!#define PARMS_MAT_SET_KIND(a, kind) ((a)|(kind) << PARMS_MAT_SHIFT)

!/*
!typedef  int (*PARMS_ILU)(parms_Mat self, parms_FactParam param, void
!			  *data, parms_Operator *op);  
!*/


type :: parms_Map
!struct parms_Map_ {
  !integer :: ref
  !parms_Table table;		!/**< contains pair (gindex,lindex) */
  !MPI_Comm    comm;		!/**< MPI communicator */
  !integer         pid;		!/**< processor ID */
  !integer         npro;		!/**< number of processors */
  !integer         lsize;		!/**< number of local variables */
  !integer         gsize;		!/**< number of global variables */
  !/*! numbering style used.
  ! *
  ! *  \f{tabular}{lcl} 
  ! *   FORTRAN &-& 1 \\
  ! *   C       &-& 0 
  ! *   \f}
  ! */
  integer :: start
  !   The number of variables associated with each vertex.
  integer :: dof!if square submat this would be nsubrow	
  !/*! style of labelling variables u_i,v_i associated with the i-th vertex.
  ! * - NONINTERLACED u_1,v_1,u_2,v_2.
  ! * - INTERLACED u_1,u_2,\cdots,v_1,v_2.
  !VARSTYPE    vtype;	
  !/*! array of size lsize, stores local variables in global labels */
  !integer,dimension(:) :: lvars		
  integer:: isserial	!BOOL	!//!< data layout on one processor?
  !/*! variables are permuted or not.
  ! *  -true: interior variables labelled first followed by interface variables */
  integer:: isperm !BOOL		
  integer:: isvecperm !BOOL	!/* used to check if vector has been permuted in gmres */	
  integer:: ispermalloc !BOOL	!/* permutation array allocated? */
  integer,dimension(:),allocatable :: perm		!//!< permutation array of size lsize.
  integer,dimension(:),allocatable :: iperm		!//!< inverse permutation array. 
  integer:: nint		!//!< number of interior variables.
  integer:: ninf		!//!< number of interface variables.
  !  \param schur_start  start of the local Schur complement which
  !  may be lesser than nbnd for the matrix with unsymmetric pattern 
  integer:: schur_start	

  !integer:: ninf_send!	//!< number of variables to be sent.
  ! param vsend  array of size ninf_send which stores variables in local indices 
  !integer,dimension(:) :: vsend		
  ! param vstable  stores variables to be sent to adjacent processors in pairs
  ! 	 (local_index, index_in_vsend) 
  !parms_Table vstable
end type parms_Map
  !|--------------------------------------------------------------------*/






!from parms_sys.h
type :: parms_FactParam !{
  !/*
  !  mc = multicoloring or not (0-not, 1-yes).
  !  ipar[1:18]  = integer array to store parameters for both
  !     arms construction (arms2) and iteration (armsol2).
  !  
  !     ipar[1]:=nlev.  number of levels (reduction processes).  see
  !     also "on return" below.
!
!       ipar[2]:= level-reordering option to be used.  
!                 if ipar[2]==0 ARMS uses a block independent set ordering
!                  with a sort of reverse cutill Mc Kee ordering to build 
!                  the blocks. This yields a symmetric ordering. 
!                  in this case the reordering routine called is indsetC
!                 if ipar[2] == 1, then a nonsymmetric ordering is used.
!                  In this case, the B matrix is constructed to be as
!                  diagonally dominant as possible and as sparse as possble.
!                  in this case the reordering routine called is ddPQ.
!    
!       ipar[3]:=bsize. Dimension of the blocks. In this version, bsize
!       is only a target block size. The size of each block can vary
!       and is >= bsize.
!    
!       ipar[4]:=iout if (iout > 0) statistics on the run are printed
!       to FILE *ft
!
!       The following are not used by arms2 -- but should set before
!       calling the preconditionning operation armsol2:
!
!       ipar[5]:= Krylov subspace dimension for last level 
!         ipar[5] == 0 means only backward/forward solve is performed
!            on last level.                   
!       ipar[6]:=  maximum # iterations on last level
!    
!       ipar[7-10] NOT used [reserved for later use] - must be 
!       set to zero. [see however a special use of ipar[7] in fgmresC.] 
!    
!       The following set method options for arms2. Their default
!       values can all be set to zero if desired.  
!    
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
!             scaling,..) are applied to the matrix before processing
!             it   
!    
!    droptol  = Threshold parameters for dropping elements in ILU
!       factorization. 
!      droptol[1:5] = related to the multilevel  block factorization
!      droptol[6:6] = related to ILU factorization of last block.
!        This flexibility is more than is really needed. one can use a
!        single parameter for all. it is preferable to use one value
!        for droptol[1:5] and another (smaller) for droptol[5:6] 
!      droptol[1] = threshold for dropping  in L [B]. See piluNEW.c:
!      droptol[2] = threshold for dropping  in U [B].
!      droptol[3] = threshold for dropping  in L^{-1} F 
!      droptol[4] = threshold for dropping  in E U^{-1} 
!      droptol[5] = threshold for dropping  in Schur complement
!      droptol[6] = threshold for dropping  in L in last block [see
!        ilutpC.c] 
!      droptol[7] = threshold for dropping  in U in last block [see
!        ilutpC.c] 
!    
!    lfil     = lfil[0:6] is an array containing the fill-in parameters.
!      similar explanations as above, namely:
!      lfil[0] = amount of fill-in kept  in L [B]. 
!      lfil[1] = amount of fill-in kept  in U [B].
!      etc..
!    
!    tolind   = tolerance parameter used by the indset function. 
!      a row is not accepted into the independent set if the *relative*
!      diagonal tolerance is below tolind. see indset function for
!      details. Good values are between 0.05 and 0.5 -- larger values
!      tend to be better for harder problems.
!    
!    Note:   The comments above are for arms. The first element for every 
!    array is for other preconditioner: ilut, iluk
!  */
  integer::    schur_start!		/* start position of the Schur complement */
  integer::    start 	!		/* start row index of the submatrix in
			!	   the whole matrix */
  integer::    n		!	/* the size of the local matrix */
  integer :: nsubrow
  integer :: nsubcol
  integer::    mc		
  integer,dimension(7)::  lfil	
  integer,dimension(18):: ipar	
  real(dp),dimension(7):: droptol
  real(dp),dimension(2):: pgfpar
  real(dp) ::tolind
  integer:: isalloc !BOOL   isalloc;
!} *parms_FactParam;
end type parms_FactParam
  !|--------------------------------------------------------------------*/



type :: parms_Operator_ops! parms_ilu_sol_vptr = {
   !pointers to procedures 
   integer:: tmp
   !  parms_ilu_sol_vcsr,
   !  parms_ilu_lsol_vcsr,
   !  parms_ilu_invs_vcsr,
   !  parms_ilu_ascend_vcsr,
   !  parms_ilu_getssize_vcsr,
   !  parms_ilu_nnz_vcsr,
   !  parms_ilu_free,
   !  parms_ilu_view
end type parms_Operator_ops
  !|--------------------------------------------------------------------*/

type :: parms_Operator!_ {
  integer :: ref
  type(parms_Operator_ops) :: ops
  type(arms_data_type),pointer::data
end type parms_Operator !};
  !|--------------------------------------------------------------------*/



! \struct parms_Mat_ops.
!  \brief struct parms_mat_ops.
type :: parms_Mat_ops
   !typedef struct parms_Mat_ops {
   !  /*! function parameter: apply
   !    \brief performs \f$y = self \times x\f$. 
   !  */
   !  int (*apply)(parms_Mat self, FLOAT *x, FLOAT *y);
   !  /*! 
   !    Divide the local equation into diagonal part and off-diagonal
   !    part. Set up communcation handler for matrix-vector product.
   !   */
   !  int (*setup)(parms_Mat self);	
   !  /*! Set the communication type
   !    \f{tabular}{ll}
   !    P2P & Copying data into a temporary buffer \\
   !    DERIVED & Creating a derived datatype as in the old version of
   !    pARMS. 
   !    \f}
   !   */
   !  int (*setcommtype)(parms_Mat self, COMMTYPE ctype);
   !
   !  /*! function parameter: mvpy
   !    \brief Performs z = beta*y + alpha*self*x.
   !  */
   !  int (*mvpy)(parms_Mat self, REAL alpha, FLOAT *x, REAL beta,
   !	      FLOAT *y, FLOAT *z);
   !
   !  /*! function parameter: getdiag
   !    \brief Get the diagonal part of the local matrix.
   !   */
   !  int (*getdiag)(parms_Mat self, void **mat);
   !  /*! function parameter: getlmat
   !    \brief Get the local matrix including diagonal and off-diagonal
   !    matrix.  
   !   */
   !  int (*getlmat)(parms_Mat self, void **mat);
   !  /*! function parameter: extend
   !    \brief Extend the submatrix mat by including equations that
   !    correspond to immediate neighbouring variables.
   !   */
   !  int (*extend)(parms_Mat self, parms_Comm handler, int start, void
   !		*mat, int *n, void **ext_mat);
   !  /*! function Parameter: mvoffd
   !    \brief The matrix-vector product for off-diagonal part.
   !   */
   !  int (*mvoffd)(parms_Mat self, FLOAT *x, FLOAT *y, int pos);
   !
   !  /*! function Parameter: matfree
   !    \brief Free the matrix mat.
   !   */
   !  int (*matfree)(parms_Mat self, void *mat);
   !
   !  /*! function Parameter: gethandler
   !      \brief Get the communication handler for mv product.
   !   */
   !  int (*gethandler)(parms_Mat self, parms_Comm *handler);
   integer::tmp
   !} *parms_Mat_ops;
end type parms_Mat_ops

  !/*! \struct parms_vcsr.
  !  \ Description: struct parms_vcsr.
!ceb identical to cs_type
  !type :: parms_vcsr !typedef struct parms_vcsr 
  !   integer::      n ! the dimension of the matrix
  !   integer:: nsubrow
  !   integer:: nsubcol 
  !   integer,dimension(:),allocatable:: nnzrow !   //!< the length of each row 
  !   integer,dimension(:),allocatable:: space !    //!<  length of space ( a work array)
  !   !! Parameter: pj = An indirect pointer to store column indices.
  !   !!  /*! parameter: pa = An indirect pointer to store corresponding nonzero entries.  
  !   type(ja_type)    ,dimension(:), allocatable :: pj ! Parameter: pj = An indirect pointer to store column indices.
  !   type(submat_type),dimension(:), allocatable :: pa ! parameter: pa = An indirect pointer to store corresponding nonzero entries. 
  !   ! *parms_vcsr;
  !end type parms_vcsr






!/*! \struct parms_Mat_
!  \Description: struct parms_Mat_
!struct parms_Mat_ {
  type :: parms_Mat
     integer:: ref
     !  parms_Mat_ops ops;
     !  void        *data;
     integer :: isserial;!BOOL
     integer :: issetup;!BOOL
     integer :: isperm;!BOOL
     integer :: isalloc;!BOOL
     !MATTYPE     type;
     !PCILUTYPE   ilutype;
     integer :: m
     integer :: n
     integer :: MM
     integer :: NN
     type(parms_Map):: is
     !  parms_vcsr  aux_data;
  end type parms_Mat
!};


  
  !static int arms_free_vcsr(parms_Operator *self)
  !!{
  !  parms_arms_data data;
  !  data = (parms_arms_data)(*self)%data;
  !  cleanARMS(data);
  !  PARMS_FREE(data);
  !  return 0;
  !!}
  
  !static int arms_view_vcsr(parms_Operator self, parms_Viewer v)
  !!{
  !!/* So far, only available for viewing last level operator */
  !  parms_arms_data data;
  !  ilutptr LU;
  !  int i, j, n, nnz, *pj;
  !  FLOAT *pa;
  !  FILE *fp;
  !  parms_ViewerGetFP(v, &fp);
  !  data = (parms_arms_data)self%data;
  !  LU = data%ilus;
  !  n = LU%n;
  !!/* L part */
  !  fprintf(fp, "L part of the last level matrix:\n");
  !  fprintf(fp, "n = %d\n", n);
  !  !for (i = 0; i < n; i++) !{
  !  do i = 0,n-1 !{
  !    nnz = LU%L%nnzrow(i);
  !    pj  = LU%L%pj(i);
  !    pa  = LU%L%pa(i);
  !    fprintf(fp, "nnzrow(%d)=%d\n", i, nnz);
  !    !for (j = 0; j < nnz; j++) !{
  !    do j = 0,nnz-1 !{
  !      fprintf(fp, "(%d,%d,%f) ", i, pj(j), pa(j));
  !    !}
  !   end do
  !  !}
  !   end do
  !!/* U part */
  !  fprintf(fp, "U part of the matrix:\n");
  !  fprintf(fp, "n = %d\n", n);
  !  !for (i = 0; i < n; i++) !{
  !  do i = 0,n-1 !{
  !    nnz = LU%U%nnzrow(i);
  !    pj  = LU%U%pj(i);
  !    pa  = LU%U%pa(i);
  !    fprintf(fp, "nnzrow(%d)=%d\n", i, nnz);
  !    !for (j = 0; j < nnz; j++) !{
  !    do j = 0,nnz-1 !{
  !      fprintf(fp, "(%d,%d,%f) ", i, pj(j), pa(j));
  !    !}
  !   end do
  !  !}
  !   end do
  !  parms_ViewerStoreFP(v, fp);
  !  return 0;
  !!}
  
  !!static void parms_arms_nnz(parms_Operator self, int *nnz_mat, int *nnz_pc)
  !subroutine parms_arms_nnz(self,nnz_mat,nnz_pc)
  !  type(parm_Operator) :: self
  !  integer,intent(inout) :: nnz_mat
  !  integer,intent(inout) :: nnz_pc
  !  type(parms_arms_data),pointer ::data
  !  data = self%data
  !  nnz_mat = data%nnz_mat
  !  nnz_pc  = data%nnz_prec
  !end subroutine parms_arms_nnz
  !!|--------------------------------------------------------------------*/



  !static int parms_arms_lsol_vcsr(parms_Operator self, FLOAT *y, FLOAT *x)
  !!{
  !parms_arms_data data;
  !p4ptr levmat;
  !ilutptr ilus;
  !int *ipar, schur_start, nlev, n;
  !data = (parms_arms_data)self%data;
  !levmat = data%levmat;
  !ilus   = data%ilus;
  !ipar   = data%ipar;
  !schur_start = data%schur_start;
  !if (ipar(0) == 0) !{
  !  Lsolp(schur_start, ilus%L, y, x); 
  !  return 0;
  !!}
  !nlev = data%nlev;
  !n = data%n;
  !PARMS_MEMCPY(x, y, n);
  !Lvsol2(x, nlev, levmat, ilus, 0) ;
  !return 0;
  !!}
  !!|--------------------------------------------------------------------*/


  !static int parms_arms_sol_vcsr(parms_Operator self, FLOAT *y,
  !			  FLOAT *x) 
  !!{
  !  parms_arms_data data;
  !  int n;
  !  data = (parms_arms_data)self%data;
  !  n = data%n;
  !  PARMS_MEMCPY(x, y, n);
  !  armsol2(x, data);
  !  return 0;
  !!}
  !|--------------------------------------------------------------------*/
  
  
  !static int parms_arms_invs_vcsr(parms_Operator self,  FLOAT *y, FLOAT *x) 
  !!{
  !  parms_arms_data data;
  !  ilutptr ilus;
  !  int *ipar, schur_start, n;
  !  data = (parms_arms_data)self%data;
  !  schur_start = data%schur_start;
  !  ilus   = data%ilus;
  !  ipar   = data%ipar;
  !  n = ilus%n;
  !  if (y == NULL || x == NULL) return 0;
  !  if (ipar(0) == 0) then
  !     invsp(schur_start, ilus, y, x);
  !  else
  !     PARMS_MEMCPY(x, y, n);
  !     SchLsol(ilus, x);
  !     SchUsol(ilus, x);
  !  end if
  !  return 0;
  !!}
  !!|--------------------------------------------------------------------*/
  
  !static int parms_arms_ascend_vcsr(parms_Operator self, FLOAT *y, FLOAT *x)  
  !!{
  !  parms_arms_data data;
  !  p4ptr levmat=NULL, last=NULL;
  !  ilutptr ilus;
  !  int *ipar, schur_start, n, nloc, lenB, first;
  !  data = (parms_arms_data)self%data;
  !  levmat = data%levmat;
  !  ilus   = data%ilus;
  !  schur_start = data%schur_start;
  !  ipar   = data%ipar;
  !  n      = data%n;
  !  if (ipar(0) == 0) then
  !    Usolp(schur_start, ilus%U, y, x);
  !    return 0;
  ! end if
  !
  !  while (levmat) !{
  !    last = levmat;
  !    levmat = levmat%next;
  !  !}
  !  levmat = last;
  !
  !  nloc=levmat%n; 
  !  lenB=levmat%nB; 
  !  first = n - nloc; 
  !  !/*-------------------- last level                                 */
  !  first += lenB; 
  !  !/*-------------------- other levels                               */
  !  while (levmat) !{
  !    nloc = levmat%n; 
  !    first -= levmat%nB;
  !    if (levmat%n) 
  !      ascend(levmat, &x(first),&x(first));
  !    !/*-------------------- right scaling */
  !    if (levmat%D2 !=  NULL) 
  !      dscale(nloc, levmat%D2, &x(first), &x(first)) ;
  !    levmat = levmat%prev; 
  !  !}
  !  return 0;
  !!}
  !|--------------------------------------------------------------------*/
  
  
  !static int parms_arms_getssize_vcsr(parms_Operator self)
  !!{
  !  parms_arms_data data;
  !  data = (parms_arms_data)self%data;
  !  return data%schur_start;
  !!}
  !!|--------------------------------------------------------------------*/
  

!/* external function protos */
!extern int parms_MatCreate_vcsr(parms_Mat self);
!extern int parms_MatCreate_dvcsr(parms_Mat self);
!extern int parms_MatFree_dvcsr(parms_Mat *self);
!extern int parms_MatView_vcsr(parms_Mat self, parms_Viewer v);
!extern int parms_MatView_dvcsr(parms_Mat self, parms_Viewer v);

end module kinddefs
