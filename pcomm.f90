!module contains everything needed for parallel communication
module pcomm
 use kinddefs, only : dp,intarray,realarray
 use blas
 use mpi
 implicit none

!interface PostSends
!module function PostSends_r,PostSends_i
!end interface

 public :: PComm_type

  type:: PComm_type 
     !protected:
     !  void CachePointers(void**ptrs,integer*strides,integer nptrs,MPI_Datatype dtype);
     integer :: isparallel
     integer :: np                ! number of processes (assumed == number of subdomains)
     integer :: my_rank           ! rank of this process
     integer :: np_universe       ! rank of this process in MPI_COMM_WORLD
     integer :: my_rank_universe  ! rank of this process in MPI_COMM_WORLD
     integer :: my_color          ! the color that we chose to take in the communicator
     integer :: my_subdomain      ! ID of this subdomain
     integer :: ishost            ! flag to denote whether this process is host
     integer, pointer,dimension(:) :: domain2rank      ! mapping of subdomain to rank
     integer, pointer,dimension(:) :: rank2domain      ! mapping of rank to subdomain
     integer :: m_verbose         ! how much to print out in normal operation
     integer :: nphnodes             ! number of phantom cv in subdomain connectivity map
     integer, pointer,dimension(:) :: subdomain_conn   ! subdomain connectivity map
     integer, pointer,dimension(:) :: nexpect          ! number of cv's expected from each subdomain
     !ptr to buffer?
     !integer, pointer,dimension(:,:) :: expect          ! list of actual cv's expected from each subdomain
     type(intarray), pointer,dimension(:) :: expect          ! list of actual cv's expected from each subdomain
     integer,pointer,dimension(:) :: nsend            ! number of cv's to send to each subdomain
     !ptr to buffer?
     !integer,pointer,dimension(:,:) :: send            ! list of actual cv's to send to each subdomain
     type(intarray),pointer,dimension(:) :: send            ! list of actual cv's to send to each subdomain
     integer:: comm_world  !type(MPI_COMM)!all mpi types are integer in fortran   ! the communicator that we will utilize.  This is the world we live in.
     
     ! support for persistent buffers to use for communication.  The reason for this is to 
     ! be able to support traditional point-to-point communication at the same time as
     ! supporting the one-sided MPI-2 RMA operations.
     integer :: buffer_bytes_max  ! maximum number of bytes allocated so far.
     
     !ceb we may want pointer to buffers here (array of buffers)
     !integer,pointer,dimension(:,:):: send_buffers     ! preallocated send buffers
     type(realarray),pointer,dimension(:):: send_buffers     ! preallocated send buffers
     !integer,pointer,dimension(:,:):: recv_buffers     ! preallocated receive buffers.
     type(realarray),pointer,dimension(:):: recv_buffers     ! preallocated receive buffers.
     
     !Timer* pack_timer;
     !Timer* unpack_timer;
     !Timer* msg_timer;
     !Timer *recv_timer;
     !Timer *send_timer;
     !Timer *wait_timer;
     ! custom datatype
     !mutable integer idebug;
     !mutable MPI_Datatype custom_type_floatcomplex;
     !mutable MPI_Datatype custom_type_doublecomplex;
     !mutable MPI_Datatype custom_type_floatdual;
     !mutable MPI_Datatype custom_type_doubledual;
     ! support variables for message send/recv
     integer :: nreq_isend
     integer :: nreq_irecv
     !   void**expect_buffer;
     !   void**send_buffer;
     integer,pointer,dimension(:):: pending_srequests !type(MPI_REQUEST)
     integer,pointer,dimension(:):: pending_rrequests !type(MPI_REQUEST)
     ! pointers to cache during a message send/recv
     integer :: m_nptrs
     type(intarray),pointer,dimension(:):: m_ptrs_i
     type(realarray),pointer,dimension(:):: m_ptrs_r
     !void**m_ptrs
     integer,pointer,dimension(:) :: m_strides
     integer:: m_cached_datatype !type(MPI_Datatype)
     ! integer AllocatePersistentBuffers(const integer nbytes_unit);
  end type PComm_type
  



contains
  
  




  !public:
  
  !PCommCreate(pcomm,const integer verb=1)
  !virtual PCommDelete(pcomm)
  !integer Initialize(pcomm,integer color = -1)
  !static integer StartMPI(pcomm,int *argc = NULL,char **argv() = NULL)
  !static integer GetNPAll(pcomm)
  !static integer GetRankAll(pcomm)
  !integer StopMPI(pcomm)
  !integer Build(pcomm,integer nnodes,integer nphcv,integer*subdomain_conn_,const integer*old2new=NULL/*optional remapping*/,PComm *pcomm_base=NULL);
  !integer DirectSet_SubdomainConn(pcomm,integer * subdomain_conn_, integer nphcv_);
  !void Modify_SubdomainConn(pcomm,integer * subdomain_conn_phn, integer iphcv);!ceb
  !integer*sn2phn(pcomm,integer nnodes);
  
  
  integer function GetComm(pcomm)
    type(PComm_type),intent(in) :: pcomm
    GetComm=pcomm%comm_world
    return 
  end function GetComm

  integer function GetDomain(pcomm)
    type(PComm_type),intent(in) :: pcomm
    GetDomain=pcomm%Rank2Domain(pcomm%my_rank)
    return 
  end function GetDomain
  
  function GetSubdomainID(pcomm) result(retval)
    type(PComm_type),intent(in) :: pcomm
    integer:: retval
    retval=pcomm%my_subdomain
    return 
  end function GetSubdomainID
  
  integer function IsHost(pcomm)
    type(PComm_type),intent(in) :: pcomm
    IsHost=pcomm%ishost
    return 
  end function IsHost
  
  integer function GetHost(pcomm)
    type(PComm_type),intent(in) :: pcomm
    GetHost=pcomm%np-1 !??
    return 
  end function GetHost
  
  integer function GetRank(pcomm)
    type(PComm_type),intent(in) :: pcomm
    GetRank=pcomm%my_rank
    return 
  end function GetRank
  
  integer function GetColor(pcomm)
    type(PComm_type),intent(in) :: pcomm
    GetColor=pcomm%my_color
    return 
  end function GetColor
  
  integer function GetNP(pcomm)
    type(PComm_type),intent(in) :: pcomm
    GetNP=pcomm%np
    return 
  end function GetNP
  
  function GetNP_Universe(pcomm) result(retval)
    type(PComm_type),intent(in) :: pcomm
    integer::retval
    retval=pcomm%np_universe
    return 
  end function GetNP_Universe

  function GetRank_Universe(pcomm) result(retval)
    type(PComm_type),intent(in) :: pcomm
    integer::retval
    retval=pcomm%my_rank_universe
    return 
  end function GetRank_Universe

  function IsParallel(pcomm) result(retval)
    type(PComm_type),intent(in) :: pcomm
    integer::retval
    retval=pcomm%isparallel
    return 
  end function IsParallel

  function GetNumberOfPhantomNodes(pcomm) result(retval)
    type(PComm_type),intent(inout) :: pcomm
    integer::retval
    retval=pcomm%nphnodes
    return 
  end function GetNumberOfPhantomNodes
  
  function GetSubdomainConnectivity(pcomm) result(retarray)
    type(PComm_type),intent(inout) :: pcomm
    integer,pointer,dimension(:)::retarray
    retarray=>pcomm%subdomain_conn
    return
  end function GetSubdomainConnectivity

  !integer *GetReceiveCounts(integer*sendcnt,integer np) const;
  !integer GetCountsAndDisplacements(integer cnt,integer *cnts,integer*displs=NULL,integer displ_offset=0) const;
  !void SetDebug(integer idebug_) const {idebug = idebug_;};
  !void RemapAndReorder(const integer nnodes,const integer *old2new,const integer istart);
  !integer Transfer_Status_Flags_to_Phantom(integer *stat_cv, const integer nnodes, const integer status_flags); 
  !uxreal GetLoadBalanceEfficiency(const integer nnodes);
  !integer ReportLoadBalanceEfficiency(const integer nnodes, const char * tag);
  !integer Dump() const;


  integer function SyncAll(pcomm)
    type(PComm_Type),intent(inout)::pcomm
    integer :: ierr
    ierr=0
    if (IsParallel(pcomm).ne.0) then 
       call MPI_Barrier(pcomm%comm_world,ierr)
    end if
    SyncAll=ierr
    return 
  end function SyncAll






  subroutine PCommCreate(pcomm,verb)
    !PComm(const integer verb)
    type(PComm_Type),intent(inout)::pcomm
    integer,intent(in)::verb
    !local vars
    integer :: istat
    pcomm%nphnodes = 0
    pcomm%isparallel = 0
    pcomm%np = 1
    pcomm%np_universe = 1
    pcomm%my_subdomain = 1
    pcomm%my_rank = 0
    pcomm%my_rank_universe = 0
    pcomm%my_color = -1
    pcomm%comm_world = MPI_COMM_NULL
    pcomm%ishost = 1
    !pcomm%idebug = 0
    pcomm%m_verbose = verb
    nullify(pcomm%send)! = NULL
    nullify(pcomm%nsend)! = NULL
    nullify(pcomm%expect)! = NULL
    nullify(pcomm%nexpect)! = NULL
    nullify(pcomm%domain2rank)! = NULL
    nullify(pcomm%rank2domain)! = NULL
    nullify(pcomm%subdomain_conn)! = NULL
    nullify(pcomm%send_buffers)! = NULL
    nullify(pcomm%recv_buffers)! = NULL
    pcomm%buffer_bytes_max = 0
    allocate(pcomm%domain2rank(2));!leaves dirty flag from ctse
    allocate(pcomm%rank2domain(1));!leaves dirty flag from ctse
    !pcomm%rank2domain(0) = 1
    pcomm%rank2domain(1) = 1
    !pcomm%domain2rank(1) = 0
    pcomm%domain2rank(2) = 0
    pcomm%nreq_isend = 0
    pcomm%nreq_irecv = 0
    nullify(pcomm%pending_srequests)! = NULL
    nullify(pcomm%pending_rrequests)! = NULL
    !  send_buffer = NULL;
    !   expect_buffer = NULL;
    pcomm%m_nptrs = 0
    nullify(pcomm%m_ptrs_r)! = NULL
    nullify(pcomm%m_ptrs_i)! = NULL
    nullify(pcomm%m_strides)! = NULL
    !custom_type_floatcomplex = 0;
    !custom_type_doublecomplex = 0;
    !custom_type_floatdual = 0;
    !custom_type_doubledual = 0;
    !pack_timer = new Timer(1, "Pack");
    !unpack_timer = new Timer(2, "Unpack");
    !msg_timer = new Timer(3, "Messaging");
    !recv_timer = new Timer(4, "Receive");
    !send_timer = new Timer(5, "Send");
    !wait_timer = new Timer(6, "Wait");
  end subroutine PCommCreate
  
  
  subroutine PCommClean(pcomm)
    type(PComm_Type),intent(inout)::pcomm
    integer i,istat
    deallocate(pcomm%nsend)
    if (associated(pcomm%send))then
       do i = 1,pcomm%np 
          deallocate(pcomm%send(i)%array,STAT=istat)
       end do
       deallocate(pcomm%send,STAT=istat)
    end if
    if (associated(pcomm%send_buffers)) then
       do i = 1,pcomm%np 
          deallocate(pcomm%send_buffers(i)%array,STAT=istat)
       end do
       deallocate(pcomm%send_buffers,STAT=istat)
       
    end if
    if (associated(pcomm%recv_buffers))then
       do i = 1,pcomm%np
          deallocate(pcomm%recv_buffers(i)%array,STAT=istat)
       end do
       deallocate(pcomm%recv_buffers,STAT=istat)
    end if
    
    deallocate(pcomm%nexpect,STAT=istat)
    if (associated(pcomm%expect))then
       do i = 1,pcomm%np 
          deallocate(pcomm%expect(i)%array,STAT=istat)
       end do
       deallocate(pcomm%expect,STAT=istat)
    end if
    deallocate(pcomm%pending_srequests,STAT=istat)
    deallocate(pcomm%pending_rrequests,STAT=istat)
    
    deallocate(pcomm%domain2rank,STAT=istat)
    deallocate(pcomm%rank2domain,STAT=istat)
    deallocate(pcomm%subdomain_conn,STAT=istat)
    
    if (pcomm%m_verbose.ne.0)then
       !msg_timer->Report();
       !pack_timer->Report();
       !unpack_timer->Report();
       !recv_timer->Report();
       !send_timer->Report();
       !wait_timer->Report();
    end if
    
    !if (custom_type_floatcomplex .ne. 0)MPI_Type_free(&custom_type_floatcomplex);!ceb causing segfaults
    !if (custom_type_doublecomplex .ne. 0) MPI_Type_free(&custom_type_doublecomplex);!ceb causing segfaults
    
    !delete(pack_timer)
    !delete(unpack_timer)
    !delete(msg_timer)
    !delete(recv_timer)
    !delete(send_timer)
    !delete(wait_timer)
  end subroutine PCommClean
  
  
  
  
  !========================================================================
  !
  !  PComm::StartMPI
  !  This is a static function in PComm, and should only be called
  !  once.
  !
  !------------------------------------------------------------------------
  !integer function StartMPI(pcomm,int* argc, char ** argv())
  function StartMPI(pcomm)  result(retval)
    type(PComm_Type),intent(inout)::pcomm
    integer ::iretval,ierr,retval,prov
    logical::bretval
    write(6,*)"PCOMM: performing one-time initialization."
    ! first, we need to make sure that MPI_Init() doesn't get called
    ! multiple times.
    call MPI_INITIALIZED(bretval,ierr)
    if (iretval.ne.0)then
       write(6,*)"PCOMM: preparing to initialize MPI."
       call MPI_INIT(iretval)
       if (iretval == MPI_SUCCESS) then
          write(6,*)"PCOMM: MPI initialized successfully: "
          call MPI_Query_thread(prov,ierr)
          if (prov == MPI_THREAD_MULTIPLE)then
             write(6,*)"MPI_THREAD_MULTIPLE provided."
          else if (prov == MPI_THREAD_SERIALIZED)then
             write(6,*)"MPI_THREAD_SERIALIZED provided."
          else if (prov == MPI_THREAD_FUNNELED)then
             write(6,*)"MPI_THREAD_FUNNELED provided."
          else if (prov == MPI_THREAD_SINGLE)then
             write(6,*)"MPI_THREAD_SINGLE provided."
          else
             write(6,*)"unable to determine threading model."
          end if
       else
          write(6,*)"PCOMM: Unable to initialize MPI."
          retval=1
          return
       end if
    end if
    retval=0
    ! return success
    return
  end function StartMPI




  ! gets the number of processes in our universe
  function GetNPAll() result(retval)
    integer ::iretval,ierr,np,retval
    logical::bretval
    call MPI_INITIALIZED(bretval,ierr)
    if (iretval.eq.0) then
       retval=1
       return 
    else
       call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
       retval=np
    end if
    return
  end function GetNPAll
  
  
  ! gets our rank in the universe
   function GetRankAll() result(retval)
    integer ::iretval,ierr,rank,retval
    logical::bretval
    call MPI_INITIALIZED(bretval,ierr)
    if (iretval.eq.0)then
       retval=0
       return
    else
       call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
       retval=rank
       return
    end if
    return
  end function GetRankAll



  
  !========================================================================
  !
  !  PComm::Initialize
  !
  !  If pcomm is set, parallel properties are inherited from it.  If not,
  !  a new communicator is created, split with the color parameter.
  !
  !------------------------------------------------------------------------
  !integer function InitializePComm(PComm* pcomm,integer icolor)
  function InitializePComm(pcomm,icolor) result(retval)
    type(PComm_Type),intent(inout)::pcomm
    integer,intent(inout)::icolor
    integer:: iretval,ierr,i,istat,retval
    logical::bretval
    write(6,*)"PCOMM: initializing."
    
    ! check to see if MPI has been initialized.
    call MPI_INITIALIZED(bretval,ierr)
    
    if (.not.bretval)then
       ! signify that MPI_Init() has been successfully called at some point.
       pcomm%isparallel = 1
       ! get the rank of this process.
       call MPI_COMM_RANK(MPI_COMM_WORLD,pcomm%my_rank_universe,iretval,ierr)
       if (iretval .ne. MPI_SUCCESS)then
          write(6,*)"PCOMM: unable to get MPI rank!"
          pcomm%my_rank_universe = -1
       end if
       
       ! get the number of processes that we will be
       ! communicating with.
       call MPI_COMM_SIZE(MPI_COMM_WORLD,pcomm%np_universe,iretval,ierr)
       if (iretval .ne. MPI_SUCCESS)then
          write(6,*)"PCOMM: unable to get MPI communicator size!"
          pcomm%np_universe = -1
       end if
       
       !if (.not.pcomm)then!??
          ! now break the communicator into colors
          pcomm%my_color = icolor
          !iretval = MPI_Comm_split(MPI_COMM_WORLD,my_color,my_rank_universe,&comm_world);
          call MPI_Comm_split(MPI_COMM_WORLD,pcomm%my_color,&
                              pcomm%my_rank_universe,pcomm%comm_world,ierr)
          if (iretval .ne. MPI_SUCCESS)then
              write(6,*)"PCOMM: unable to break MPI_COMM_WORLD into a split communicator!"
           end if

          ! get our rank in the split communicator
           call MPI_COMM_RANK(pcomm%comm_world,pcomm%my_rank,iretval,ierr)
          if (iretval .ne. MPI_SUCCESS)then
              write(6,*)"PCOMM: unable to get MPI rank!"
              pcomm%my_rank = -1
           end if

           ! get the number of processes that we will be
           ! communicating with, in the split communicator
           !
           call MPI_COMM_SIZE(pcomm%comm_world,pcomm%np,iretval,ierr)
           if (iretval .ne. MPI_SUCCESS)then
              write(6,*)"PCOMM: unable to get MPI communicator size!"
              pcomm%np = -1
           end if
        !else
           ! if we are inheriting our properties from another PComm object,
           ! do that here.
           !isparallel = pcomm->IsParallel();
           !my_color = pcomm->GetColor();
           !comm_world = pcomm->GetComm();
           !my_rank = pcomm->GetRank();
           !np = pcomm->GetNP();
        !end if

        ! for now, the subdomain that we will deal with is 
        ! the same as the process rank.
        pcomm%my_subdomain = pcomm%my_rank + 1
  
        ! generate the domain2rank map.
        deallocate(pcomm%domain2rank,STAT=istat)
        allocate(pcomm%domain2rank(pcomm%np + 1),STAT=istat)
        !do i = 1,pcomm%np
        do i = 2,pcomm%np
           pcomm%domain2rank(i) = i - 1
        end do
      
        ! generate the rank2domain map
        deallocate(pcomm%rank2domain,STAT=istat)
        allocate(pcomm%rank2domain(pcomm%np),STAT=istat)
        !do i = 0,pcomm%np-1
        do i = 1,pcomm%np
           pcomm%rank2domain(i) = i + 1
        end do
                ! determine if we will act as host or not
                if (pcomm%my_rank == GetHost(pcomm))then
         pcomm%ishost = 1
      else
         pcomm%ishost = 0
      end if
   else
      ! if MPI is not running by this point, we have to assume that we are running serial, so leave all values as set up in the constructor.
   end if

   ! return success
   retval=0
   return
 end function InitializePComm




 function StopMPI(pcomm) result(retval)
   type(PComm_Type),intent(inout)::pcomm
   integer:: ierr,retval
   call MPI_COMM_FREE(pcomm%comm_world)
   call MPI_FINALIZE(ierr)
   retval=0
   return 
 end function StopMPI
 

 function Dump(pcomm) result(retval)
   type(PComm_Type),intent(inout)::pcomm
   integer ::i,j,retval
   write(6,*)"PCOMM DUMP:"
   write(6,*)" np = ",pcomm%np
   write(6,*)" expect lists:"
   do j = 1,pcomm%np
      write(6,*)"Subdomain ",j,pcomm%nexpect(j)
      do i = 1,pcomm%nexpect(j)
         write(6,*)" ",pcomm%expect(j)%array(i)	
      end do
      write(6,*)""
   end do
   write(6,*)" send lists:"
   do j = 1,pcomm%np
      write(6,*)"Subdomain ",j,pcomm%nsend(j)
      do i = 1,pcomm%nsend(j)
         write(6,*)" ",pcomm%send(j)%array(i)
      end do
      write(6,*)""
   end do
   retval=0
   return
 end function Dump
 
 
 !integer(:) function sn2phn(integer nnodes)
 !  type(PComm_Type),intent(inout)::pcomm
 !  integer,intent(in)::nnodes
 !  integer,dimension(:),allocatable::map! = NULL;
 !  integer i,inode,istat
 !  allocate(map(nnodes+1),STAT=istat)
 !  do i = 1,pcomm%nphnodes
 !     inode = subdomain_conn(3*i + 0)
 !     map(inode) = i
 !  end do
 !  sn2phn=map
 !  return 
 !end function sn2phn
 
 
 
 
 
 
 
 !
 ! this routine allocates storage for the buffers that we will use for 
 ! sending and receiving.
 !
 ! In the case of one-sided RMA operations, the memory windows will be
 ! established here.
 !
 function AllocatePersistentBuffers(pcomm,nbytes_unit) result(retval)
   type(PComm_Type),intent(inout)::pcomm
   integer,intent(in) :: nbytes_unit
   integer :: idomain,istat,retval,size
   if (.not.associated(pcomm%recv_buffers)) then
      allocate(pcomm%recv_buffers(pcomm%np + 1),STAT=istat)
   end if
   if (.not.associated(pcomm%send_buffers)) then
      allocate(pcomm%send_buffers(pcomm%np + 1),STAT=istat)
   end if
   if (nbytes_unit > pcomm%buffer_bytes_max)then
      pcomm%buffer_bytes_max = nbytes_unit
      do idomain = 1,pcomm%np
         if (pcomm%nexpect(idomain) > 0)then
            size=pcomm%nexpect(idomain)*pcomm%buffer_bytes_max
            allocate(pcomm%recv_buffers(idomain)%array(size),STAT=istat)
         end if
         if (pcomm%nsend(idomain) > 0)then
            size=pcomm%nsend(idomain)*pcomm%buffer_bytes_max
            allocate(pcomm%send_buffers(idomain)%array(size),STAT=istat)
         end if
      end do
   end if
   retval=0
   return
 end function AllocatePersistentBuffers
 
 
 
 
 
 function Handle_Reordering(pcomm,nnodes,old2new) result(retval)
   type(PComm_Type),intent(inout)::pcomm
   integer,intent(in):: nnodes
   integer,intent(inout),pointer,dimension(:)::old2new
   !local vars
   integer idomain,i,icv,istat,retval,ierr
   integer,pointer,dimension(:):: old2new_other! = tmemdup(old2new,nnodes+1);
   integer :: nptrs
   !real(dp),pointer,dimension(:):: ptrs
   integer,pointer,dimension(:):: ptrs
   integer,pointer,dimension(:)::strides
   !old2new_other => tmemdup(old2new,nnodes+1)
   
   ! we need to remap the subdomain_conn list.  This one is a little difficult, because it references
   ! nodes in other subdomains, whose numbers have also changed.
   
   ! perform a message swap in order to update each cv with the old2new entry from the appropriate
   ! subdomain.
   ptrs => old2new_other
   nptrs = 1
   !strides() = {1}
   !ierr = Aswap(pcomm,ptrs,strides,nptrs) !uses the expect/send lists, which are not remapped yet.
   ierr = Aswap_i(pcomm,ptrs,strides,nptrs) !uses the expect/send lists, which are not remapped yet.
   
   ! now remap entries 0 and 2 of subdomain_conn.
   do i = 1,pcomm%nphnodes
      icv = pcomm%subdomain_conn(3*i + 0)
      !_ASSERT(old2new(icv) > 0)
      !_ASSERT(old2new(icv) <=  nnodes)
      !_ASSERT(old2new_other(icv) > 0)
      pcomm%subdomain_conn(3*i + 0) = old2new(icv)
      pcomm%subdomain_conn(3*i + 2) = old2new_other(icv)
   end do
   
   deallocate(old2new_other,STAT=istat)
   
   ! since both the send and expect lists hold references to cv's
   ! LOCAL TO THIS SUBDOMAIN, we can just directly remap both of
   ! them without knowledge of the reordering in other subdomains.
   ! Their ordering in relation to each other is unchanged.
   
   ! first, directly remap the send lists
   do idomain = 1,pcomm%np
      !_ASSERT(send(idomain) .ne. NULL);
      do i = 1,pcomm%nsend(idomain)
!         icv = send(idomain)(i);
         icv = pcomm%send(idomain)%array(i)
         !_ASSERT(icv > 0);
         !if (nnodes >= 0){_ASSERT(icv <= nnodes);}
         !_ASSERT(old2new(icv) > 0);
         !if (nnodes >= 0){_ASSERT(old2new(icv) <= nnodes);}
         pcomm%send(idomain)%array(i) = old2new(icv)
      end do
   end do
   
   ! next, directly remap the expect lists
   do idomain = 1,pcomm%np
      !_ASSERT(expect(idomain) .ne. NULL);
      do i = 1,pcomm%nexpect(idomain)
         icv = pcomm%expect(idomain)%array(i)
         !_ASSERT(icv > 0);
         !if (nnodes >= 0) {_ASSERT(icv <= nnodes);}
         !_ASSERT(old2new(icv) > 0);
         !if (nnodes >= 0) {_ASSERT(old2new(icv) <= nnodes);}
         pcomm%expect(idomain)%array(i) = old2new(icv)
      end do
   end do
   retval=0
   return
 end function Handle_Reordering
 
 
 
 
! !
! ! return receive counts given the send counts
! !
! ! annoyance: np is a redundant argument here.  Hides the member variable np.
! integer(:) function GetReceiveCounts(pcomm,sendcnt,np)
!   type(PComm_Type),intent(inout)::pcomm
!   integer,dimension(:) :: sendcnt! = NULL;
!   integer,dimension(:),allocatable :: recvcnt! = NULL;
!   integer :: idomain,nreq_r,nreq_s,istat,ierr
!   MPI_Request,dimension(:),allocatable :: srequest! = NULL;
!   MPI_Request,dimension(:),allocatable :: rrequest! = NULL;
!   MPI_Status, dimension(:),allocatable :: statuses! = NULL;
!   nreq_r=0
!   nreq_s=0
!   ! now that we have the send counts, we can notify the receivers what they
!   ! will get
!   allocate(recvcnt(np+1),STAT=istat)
!   if (.not.pcomm%isparallel)then
!      recvcnt(1) = sendcnt(1)
!   else
!      allocate(srequest(np),STAT=istat)
!      allocate(rrequest(np),STAT=istat)
!      allocate(statuses(np),STAT=istat)
!      do idomain = 1,pcomm%np
!         if (idomain .ne. pcomm%my_subdomain)then
!            call MPI_Isend(sendcnt(idomain),1,MPI_INT,pcomm%domain2rank(idomain),pcomm%my_subdomain,pcomm%comm_world,srequest(nreq_s),ierr)
!            call MPI_Irecv(recvcnt(idomain),1,MPI_INT,pcomm%domain2rank(idomain),idomain,pcomm%comm_world,rrequest(nreq_r),ierr)
!            nreq_r=nreq_r+1
!            nreq_s=nreq_s+1
!         else
!            recvcnt(idomain) = sendcnt(idomain)
!         end if
!      end do
!      
!      call MPI_Waitall(nreq_r,rrequest,statuses,ierr)
!      call MPI_Waitall(nreq_s,srequest,statuses,ierr)  
!      
!      deallocate(srequest,STAT=istat)
!      deallocate(rrequest,STAT=istat)
!      deallocate(statuses,STAT=istat)
!   end if
!   GetReceiveCounts=recvcnt
!   return
! end function GetReceiveCounts
 



 !
 ! this is a routine to bypass the normal Build() function and set subdomain_conn directly.  Should
 ! only be used when utilizing PComm as a storage device, since you cannot communicate without the
 ! expect/send lists being built (as Build() does).
 !
 function DirectSet_SubdomainConn(pcomm,subdomain_conn_, nphnodes_) result(retval)
   type(PComm_Type),intent(inout)::pcomm
   integer,intent(in),dimension(:)::subdomain_conn_
   integer,intent(in) :: nphnodes_
   !_ASSERT(subdomain_conn == NULL);
   !_ASSERT(nphnodes == 0);
   integer::retval
   pcomm%subdomain_conn = subdomain_conn_
   pcomm%nphnodes = nphnodes_
   retval=0
   return
 end function DirectSet_SubdomainConn

 !This function allows modification of a specofic phantom node entry in the 
 ! subdomain connectivity array
 subroutine Modify_SubdomainConn(pcomm, subdomain_conn_phn, iphcv)
   type(PComm_Type),intent(inout)::pcomm
   integer,intent(in),dimension(:)::subdomain_conn_phn !index from 0
   integer,intent(in)::iphcv
   !pcomm%subdomain_conn(3*iphcv+0)=subdomain_conn_phn(0)
   !pcomm%subdomain_conn(3*iphcv+1)=subdomain_conn_phn(1)
   !pcomm%subdomain_conn(3*iphcv+2)=subdomain_conn_phn(2)
   pcomm%subdomain_conn(3*(iphcv-1)+1)=subdomain_conn_phn(1)
   pcomm%subdomain_conn(3*(iphcv-1)+2)=subdomain_conn_phn(2)
   pcomm%subdomain_conn(3*(iphcv-1)+3)=subdomain_conn_phn(3)
 end subroutine Modify_SubdomainConn









 ! Build the expect and send lists
 ! A word on "universe" communicators...subdomain_conn has subdomain ID's based specifically on the
 ! "world" that they saw themselves in at partition time.  In the case that we want to build a pcomm
 ! based on that subdomain_conn, then we need to have a way of mapping the subdomain ID's from the
 ! world to the universe, so that subdomain_conn can be modified.
 !
 ! To facilitate this, we pass in the pcomm that subdomain_conn_ is based on so that we can do the
 ! remapping.
 !
 integer function Build(pcomm,nnodes,nphnodes_,subdomain_conn_,old2new,pcomm_base)
   type(PComm_type),intent(inout)::pcomm
   integer,intent(inout)::nnodes,nphnodes_
   integer,intent(inout),pointer,dimension(:)::subdomain_conn_
   integer,intent(inout),pointer,dimension(:)::old2new!/*optional remapping*/
   type(PComm_type),intent(inout),optional::pcomm_base!/* PComm that subdomain_conn_ is actually based on */
   !local vars
   integer :: i,idomain,istat,domain,nreq_r,nreq_s,irank,ierr
   type(intarray),pointer,dimension(:):: expect_other
   integer,dimension(:),allocatable::srequest,rrequest!type(MPI_REQUEST)
   integer,dimension(:),allocatable::statuses!type(MPI_STATUS)
   integer,dimension(:),allocatable::uranks
   
   allocate(pcomm%subdomain_conn(3*(nphnodes_+1)),STAT=istat)
   if (associated(subdomain_conn_))then
      ierr=tmemcpy(pcomm%subdomain_conn,subdomain_conn_,int(3*(nphnodes_+1)))
   else
      ierr=tblank(pcomm%subdomain_conn,int(3*(nphnodes_+1)))
   end if
      
   !nullify(srequest)
   !nullify(rrequest)
   !nullify(statuses)
   allocate(srequest(pcomm%np),STAT=istat)
   allocate(rrequest(pcomm%np),STAT=istat)
   allocate(statuses(pcomm%np),STAT=istat)
   
   ! set nphnodes explicitly
   pcomm%nphnodes = nphnodes_
   
   ! take care of remapping required for universe connectivity.  subdomain_conn
   ! is given with respect to pcomm_base, and we need to adjust the subdomain
   ! numbers so that it works for this pcomm.
   if (present(pcomm_base) .and. IsParallel(pcomm).ne.0 .and. associated(subdomain_conn_))then
      allocate(uranks(GetNP(pcomm)),STAT=istat)
      call MPI_Allgather(pcomm%my_rank_universe,1,MPI_INTEGER,uranks,&
                         1,MPI_INTEGER,GetComm(pcomm_base),ierr)
      do i = 1,pcomm%nphnodes
         idomain = pcomm%subdomain_conn(3*i + 1)
         irank = pcomm_base%Domain2Rank(idomain)
         irank = uranks(irank)
         idomain = pcomm%Rank2Domain(irank)
         pcomm%subdomain_conn(3*i + 1) = idomain
      end do
      deallocate(uranks,STAT=istat)
   end if
   
   
   ! allocate the expect lists
   if (associated(pcomm%expect))then
      !do i=0,pcomm%np
      do i=1,pcomm%np
         deallocate(pcomm%expect(i)%array,STAT=istat)
      end do
   end if
   deallocate(pcomm%expect,STAT=istat)
   allocate(pcomm%expect(pcomm%np + 1),STAT=istat)
   
   if (associated(expect_other))then
      !do i=0,pcomm%np
      do i=1,pcomm%np
         deallocate(expect_other(i)%array,STAT=istat)
      end do
   end if
   deallocate(expect_other)
   allocate(expect_other(pcomm%np + 1),STAT=istat)
   
   deallocate(pcomm%nexpect,STAT=istat)
   allocate(pcomm%nexpect(pcomm%np + 1),STAT=istat)
   
   deallocate(pcomm%nsend,STAT=istat)
   allocate(pcomm%nsend(pcomm%np + 1),STAT=istat)
   
   do idomain = 1,pcomm%np
      allocate(pcomm%expect(idomain)%array(pcomm%nphnodes+1))
      allocate(expect_other(idomain)%array(pcomm%nphnodes+1))
   end do
   
   ! build the expect lists
   do i = 1,pcomm%nphnodes
      !idomain = pcomm%subdomain_conn(3*i + 1)
      idomain = pcomm%subdomain_conn(3*(i-1) + 2)
      !_ASSERT(idomain > 0);
      !_ASSERT(idomain <= np);
      pcomm%nexpect(idomain)=pcomm%nexpect(idomain)+1
      ! record the node number in this domain that this particular
      ! phantom node references
      !pcomm%expect(idomain)%array(pcomm%nexpect(idomain)) = pcomm%subdomain_conn(3*i + 0)
      pcomm%expect(idomain)%array(pcomm%nexpect(idomain)) = pcomm%subdomain_conn(3*(i-1) + 1)
      ! record the node number in idomain that this particular
      ! phantom node referenes
      !expect_other(idomain)%array(pcomm%nexpect(idomain)) = pcomm%subdomain_conn(3*i + 2)
      expect_other(idomain)%array(pcomm%nexpect(idomain)) = pcomm%subdomain_conn(3*(i-1) + 3)
   end do

   ! shrink the expect lists
   do idomain = 1, pcomm%np
      !_ASSERT(nexpect(idomain) >= 0);
      !_ASSERT(nexpect(idomain) <= nphnodes);
      allocate(pcomm%expect(idomain)%array(pcomm%nexpect(idomain) + 1),STAT=istat)
      allocate(expect_other(idomain)%array(pcomm%nexpect(idomain) + 1),STAT=istat)
   end do

   ! exchange messages to determine nsend()
   nreq_r = 0
   nreq_s = 0
   do idomain = 1,pcomm%np
      if (idomain .ne. pcomm%my_subdomain) then
         call MPI_Isend(pcomm%nexpect(idomain),1,MPI_INTEGER, &
                        pcomm%domain2rank(idomain),pcomm%my_subdomain, &
                        pcomm%comm_world,srequest(nreq_s),ierr)
         call MPI_Irecv(pcomm%nsend(idomain),1, MPI_INTEGER, &
                        pcomm%domain2rank(idomain),idomain,pcomm%comm_world, &
                        rrequest(nreq_r),ierr)
         nreq_r=nreq_r+1
         nreq_s=nreq_s+1
      else
         pcomm%nsend(idomain) = pcomm%nexpect(idomain)
      end if
   end do
   
   if (nreq_r > 0) call MPI_Waitall(nreq_r,rrequest,statuses,ierr)
   if (nreq_s > 0) call MPI_Waitall(nreq_s,srequest,statuses,ierr)
   
   ! allocate memory for the send lists
   !deallocate(send);!ceb
   if(.not.associated(pcomm%send)) allocate(pcomm%send(pcomm%np + 1),STAT=istat)
   do idomain = 1,pcomm%np
      if(.not.allocated(pcomm%send(idomain)%array))then
         allocate(pcomm%send(idomain)%array(pcomm%nsend(idomain)+1),STAT=istat)
      else
         allocate(pcomm%send(idomain)%array(pcomm%nsend(idomain)+1),STAT=istat)
      end if
   end do
   
   ! now, exchange messages to ship the send lists to each
   ! appropriate process
   nreq_r = 0
   nreq_s = 0
   do idomain = 1,pcomm%np
      if (pcomm%nexpect(idomain) > 0) then
         if (idomain .ne. pcomm%my_subdomain)then
            call MPI_ISEND(expect_other(idomain)%array,pcomm%nexpect(idomain)+1,&
                           MPI_INTEGER,pcomm%domain2rank(idomain),&
                           pcomm%my_subdomain,pcomm%comm_world,&
                           srequest(nreq_s),ierr)
            nreq_s=nreq_s+1
         else
            !_ASSERT(nexpect(idomain) == nsend(idomain));
            ierr=tmemcpy(pcomm%send(idomain)%array, expect_other(idomain)%array, pcomm%nexpect(idomain)+1)
         end if
      end if
      if (pcomm%nsend(idomain) > 0)then
         if (idomain .ne. pcomm%my_subdomain)then
            call MPI_IRECV(pcomm%send(idomain)%array, pcomm%nsend(idomain)+1, &
                           MPI_INTEGER, pcomm%domain2rank(idomain), &
                           idomain,pcomm%comm_world,rrequest(nreq_r),ierr)
            nreq_r=nreq_r+1
         end if
      end if
   end do
   
   if (nreq_r > 0) call MPI_WAITALL(nreq_r,rrequest,statuses,ierr)
   if (nreq_s > 0) call MPI_WAITALL(nreq_s,srequest,statuses,ierr)
   
   !if (associated(old2new)) then
   !   do idomain = 1,pcomm%np
   !      do i = 1,pcomm%nsend(idomain)
   !         !_ASSERT(send(idomain)(i) > 0);
   !      end do
   !   end do
   !end if

   ! clean up
   do idomain = 1,pcomm%np
      deallocate(expect_other(idomain)%array,STAT=istat)
   end do
   deallocate(expect_other,STAT=istat)
   deallocate(statuses,STAT=istat)
   deallocate(rrequest,STAT=istat)
   deallocate(srequest,STAT=istat)
   Build=0
   return
 end function Build






 function RebuildSubdomainConnectivity(pcomm) result(retval)
   type(PComm_Type),intent(inout)::pcomm
   !local vars
   type(intarray),pointer,dimension(:) :: send_other ! = NULL;
   integer::idomain, nreq_r,nreq_s,i,cnt,retval,istat
   integer,pointer,dimension(:):: srequest ! = NULL;!type(MPI_Request)
   integer,pointer,dimension(:):: rrequest ! = NULL;!type(MPI_Request)
   integer,pointer,dimension(:):: statuses ! = NULL;!type(MPI_Status)
   
   !if (.not.IsMultiprocess(pcomm))then
      write(6,*)"PCOMM:RebuildSubdomainConnectivity aborted; single process."
      retval=0
      return
   !end if
   
   allocate(srequest(pcomm%np),STAT=istat)
   allocate(rrequest(pcomm%np),STAT=istat)
   allocate(statuses(pcomm%np),STAT=istat)
   
   ! first, determine the number of phantom cv's we have.
   pcomm%nphnodes = 0
   do idomain = 1,pcomm%np
      if (idomain .ne. pcomm%my_subdomain)then
         pcomm%nphnodes = pcomm%nphnodes + pcomm%nexpect(idomain)
      end if
   end do
   
   ! allocate the send_other lists
   allocate(send_other(pcomm%np + 1),STAT=istat)
   do idomain = 1,pcomm%np
      allocate(send_other(idomain)%array(pcomm%nexpect(idomain)+1),STAT=istat)
   end do
   
   ! now, exchange messages to ship the send lists to each
   ! appropriate process
   nreq_r = 0
   nreq_s = 0
   do idomain = 1,pcomm%np
      if (pcomm%nsend(idomain) > 0) then
         call MPI_ISEND(pcomm%send(idomain)%array,pcomm%nsend(idomain)+1,MPI_INTEGER,&
                        pcomm%domain2rank(idomain),pcomm%my_subdomain,&
                        pcomm%comm_world,srequest(nreq_s))
         nreq_s=nreq_s+1
      end if
      if (pcomm%nexpect(idomain) > 0)then
         call MPI_IRECV(send_other(idomain)%array,pcomm%nexpect(idomain)+1,MPI_INTEGER,&
                        pcomm%domain2rank(idomain),idomain,pcomm%comm_world,&
                        rrequest(nreq_r))
         nreq_r=nreq_r+1
      end if
   end do
   call MPI_WAITALL(nreq_r,rrequest,statuses)
   call MPI_WAITALL(nreq_s,srequest,statuses)
   
   ! allocate memory for subdomain_conn
   deallocate(pcomm%subdomain_conn,STAT=istat)
   allocate(pcomm%subdomain_conn(3*(pcomm%nphnodes+1)),STAT=istat)
   
   ! fill in subdomain_conn
   cnt = 0
   do idomain = 1,pcomm%np
      do i = 1,pcomm%nexpect(idomain)
         cnt=cnt+1
         pcomm%subdomain_conn(3*(cnt-1) + 1) = pcomm%expect(idomain)%array(i)
         pcomm%subdomain_conn(3*(cnt-1) + 2) = idomain
         pcomm%subdomain_conn(3*(cnt-1) + 3) = send_other(idomain)%array(i)	  
      end do
   end do
   !_ASSERT(cnt == nphnodes);
   
   ! clean up
   do idomain = 1,pcomm%np
      deallocate(send_other(idomain)%array,STAT=istat)
   end do
   deallocate(send_other,STAT=istat)
   deallocate(statuses,STAT=istat)
   deallocate(rrequest,STAT=istat)
   deallocate(srequest,STAT=istat)
   retval=0
   return
 end function RebuildSubdomainConnectivity
 
 
 
 
 
 function ReceiveAny() result(retval)
   integer:: ai,retval
!   real :: af
!   real(dp) :: ad
!
!   char ac = 0
!   ai=0
!   af=0
!   ad=0
!   !complex<float> afc=0.;
!   !complex<double> adc=0.;
!   !Dualn<float,DUAL_SIZE> afd;
!   !Dualn<double,DUAL_SIZE> add;
!   !Dual1<float>  afd = 0.0;
!   !Dual1<double> add = 0.0;
!   !_ASSERT(m_cached_datatype .ne. 0);
!   if (m_cached_datatype == MPI_CHARACTER) retval= ReceiveAny(ac)
!   if (m_cached_datatype == MPI_INTEGER) return ReceiveAny(ai);
!   if (m_cached_datatype == MPI_DOUBLE_PRECISION) return ReceiveAny(ad);
!   !if (m_cached_datatype == MPI_FLOAT) return ReceiveAny(af);
!   !if (m_cached_datatype == custom_type_floatcomplex) return ReceiveAny(afc);
!   !if (m_cached_datatype == custom_type_doublecomplex) return ReceiveAny(adc);
!   !if (m_cached_datatype == custom_type_floatdual) return ReceiveAny(afd);
!   !if (m_cached_datatype == custom_type_doubledual) return ReceiveAny(add);
!   
!   write(6,*)"FATAL ERROR in RecieveAny(): unrecognized type."
   retval=-1
   return
 end function ReceiveAny
 


 
 ! make sure that pending sends and receives are finished
 function FinishPendingComms(pcomm) result(retval)
   type(PComm_Type),intent(inout)::pcomm
   integer :: pending,ierr,retval 
   integer,pointer,dimension(:):: statuses ! = NULL;!type(MPI_STATUS)
   integer :: i,istat
   pending=1
   !wait_timer->On();
   ! if we have no pointers cached, we have no communications pending.
   if (pcomm%m_nptrs .eq. 0) then
      !_ASSERT(nreq_irecv == 0);
      pending = 0
   end if
   ! finish out the receives
   do!(pending)
      !pending = ReceiveAny(pcomm)
      ! something really bad has happened.
      if (pending < 0)then
         write(6,*)"FATAL ERROR in FinishPendingComms....the communications were not finished."
      end if
      if (pending.eq.0) exit !break
   end do
   !msg_timer->On();
   ! finish out the sends
   if (pcomm%nreq_isend > 0)then
      !_ASSERT(pending_srequests .ne. NULL);
      allocate(statuses(pcomm%nreq_isend),STAT=istat)
      call MPI_WAITALL(pcomm%nreq_isend,pcomm%pending_srequests,statuses,ierr)
      deallocate(statuses,STAT=istat)
   end if
   
   ! free all of the send bufffers
   if (associated(pcomm%send_buffers))then
      do i = 1,pcomm%np
         deallocate(pcomm%send_buffers(i)%array,STAT=istat)
      end do
      deallocate(pcomm%send_buffers,STAT=istat)
   end if

   ! free the pending requests
   deallocate(pcomm%pending_srequests,STAT=istat)
   pcomm%nreq_isend = 0
   ! free the cached pointers
   deallocate(pcomm%m_ptrs_r,STAT=istat)
   deallocate(pcomm%m_ptrs_i,STAT=istat)
   deallocate(pcomm%m_strides,STAT=istat)
   pcomm%m_nptrs = 0
   pcomm%m_cached_datatype = 0
   !wait_timer->Off();
   !msg_timer->Off();
   retval=0
   return
 end function FinishPendingComms
  
 
 
 function GetCountsAndDisplacements(pcomm,cnt,cnts,displs,displ_offset) result(retval)
   type(PComm_Type),intent(inout)::pcomm
   integer,intent(inout):: cnt,displ_offset
   integer,intent(inout),pointer,dimension(:)::cnts,displs !(index from 0)
   !internal vars
   integer i,ierr,retval
   !if (IsMultiprocess(pcomm)) then
   if (IsParallel(pcomm).eq.1) then
      call MPI_Allgather(cnt,1,MPI_INTEGER,cnts,1,MPI_INTEGER,pcomm%comm_world,ierr)
   else
      cnts(1) = cnt
   end if
   !if (displs .ne. NULL)then
   if (associated(displs))then
      !displs(0) = displ_offset;
      displs(1) = displ_offset
      !do i = 1,pcomm%np 
      do i = 2,pcomm%np+1 
         displs(i) = displs(i-1) + cnts(i-1)
      end do
   end if
   retval=0
   return
 end function GetCountsAndDisplacements
 
 
 
 function RemapAndReorder(pcomm,nnodes,old2new,istart) result(ierr)
   type(PComm_Type),intent(inout)::pcomm
   integer,intent(inout):: nnodes,istart
   integer,intent(inout),pointer,dimension(:)::old2new
   integer::ierr
   ierr=Handle_Reordering(pcomm,nnodes,old2new)
   return
 end function RemapAndReorder
 
 




!
! compute the load imbalance, significant on the host process only.
!
  function GetLoadBalanceEfficiency(pcomm,nnodes) result(retval)
    type(PComm_type),intent(inout) :: pcomm
    integer,intent(in) ::nnodes
    integer:: nnodes_this,np,nnodes_max,nnodes_sum,i,ierr,istat
    integer,dimension(:),allocatable:: nnodes_array! = NULL;
    real(dp)::li,retval
    nnodes_this=nnodes
    ! if we are not in parallel mode, we have perfect load balancing
    if (IsParallel(pcomm).eq.0) then
       retval=1.0
       return
    end if
    li = -1.0
    np = GetNP(pcomm)
    allocate(nnodes_array(pcomm%np),STAT=istat);
    
    call MPI_GATHER(nnodes_this,1,MPI_INTEGER,nnodes_array,1,MPI_INTEGER,&
                    GetHost(pcomm),pcomm%comm_world,ierr)
    
    if (IsHost(pcomm).eq.1) then
       nnodes_max = 0
       nnodes_sum = 0
       !do i = 0,pcomm%np-1 
       do i = 1,pcomm%np 
          nnodes_max = max(nnodes_max,nnodes_array(i))
       end do
       !do i = 0,pcomm%np-1 
       do i = 1,pcomm%np 
          nnodes_sum =nnodes_sum + nnodes_array(i)
       end do
       li = nnodes_sum/(pcomm%np*nnodes_max)!make sure this is real eval not int
    end if
    deallocate(nnodes_array,STAT=istat)
    retval=li
    return
  end function GetLoadBalanceEfficiency



  
  ! wait for particular requests to complete
  function CompleteRequests(pcomm,nreq,requests) result(retval)
    type(PComm_Type),intent(inout)::pcomm
    integer,intent(inout):: nreq
    integer,intent(inout),dimension(:)::requests!type(MPI_REQUEST)
    !local vars
    integer,dimension(:),allocatable::statuses!type(MPI_STATUS)
    integer ierr,istat,retval
    ! if we are not in parallel mode, return as there is nothing to be done
    if (IsParallel(pcomm).eq.0) then
       retval=0
       return
    end if
    allocate(statuses(nreq),STAT=istat)
    call MPI_WAITALL(nreq,requests,statuses,ierr)
    deallocate(statuses,STAT=istat)
    retval=0
    return
  end function CompleteRequests


  function Die(icode) result(retval)
    integer,intent(inout):: icode
    integer:: flag,ierr,retval
    logical::bflag
    call MPI_INITIALIZED(bflag,ierr)
    if (.not.bflag) then 	
       call MPI_FINALIZE(ierr)
    end if
    !exit(icode)
    retval=icode
    return
  end function Die












!=========================================================================================




  
!====================================================================


  !integer AllocatePersistentBuffers(const integer nbytes_unit);
 
  
  ! return receive counts given the send counts
  !integer(:) function GetReceiveCounts(pcomm,*sendcnt,integer my_rank,integer np,MPI_Comm*comm)
  function GetReceiveCounts(pcomm,sendcnt,my_rank,np) result (retarray)
    type(PComm_Type),intent(inout)::pcomm
    integer,intent(inout),dimension(:)::sendcnt
    integer,intent(in)::my_rank,np
    integer,pointer,dimension(:)::retarray
    !// annoyance: np is a redundant argument here.  Hides the member variable np.
    !local vars
    integer :: idomain,nreq_r,nreq_s,ierr,istat
    integer ,pointer,dimension(:) :: recvcnt
    integer,dimension(:),allocatable::srequest,rrequest!type(MPI_Request)  
    integer,dimension(:),allocatable::statuses!type(MPI_Status)
    
    nreq_r = 0
    nreq_s = 0
    
    ! now that we have the send counts, we can notify the receivers what they will get
    allocate(recvcnt(np+1),STAT=istat)
    if (IsParallel(pcomm).eq.0) then
       recvcnt(1) = sendcnt(1)
    else
       allocate(srequest(np),STAT=istat)
       allocate(rrequest(np),STAT=istat)
       allocate(statuses(np),STAT=istat)
       do idomain = 1,np
          if (idomain .ne. pcomm%my_subdomain)then
             call MPI_ISEND(sendcnt(idomain),1,MPI_INT,pcomm%domain2rank(idomain),&
                  pcomm%my_subdomain,pcomm%comm_world,srequest(nreq_s),ierr)
             call MPI_IRECV(recvcnt(idomain),1,MPI_INT,pcomm%domain2rank(idomain),&
                  idomain,pcomm%comm_world,rrequest(nreq_r),ierr)
             nreq_r=nreq_r+1
             nreq_s=nreq_s+1
          else
             recvcnt(idomain) = sendcnt(idomain)
          end if
       end do
       call MPI_WAITALL(nreq_r,rrequest,statuses)
       call MPI_WAITALL(nreq_s,srequest,statuses)  
       deallocate(srequest)
       deallocate(rrequest)
       deallocate(statuses)
    end if
    retarray=>recvcnt
    return
  end function GetReceiveCounts
  
  
  
  

!
! performs a send operation
!
! note: the sendcnt() array is indexed by rank. , sndcnt(i) is the number of 
!       data items to be sent to rank i (and remember that i is zero-based).
!       likewise for sendbuf(i).
!

!template <class TheType>
!integer GenericSend(TheType **sendbuf,integer*sendcnt,
!		  integer**recvcnt_,TheType***recvbuf_,
!		  MPI_Datatype mpi_type, integer blocksize,
!		  MPI_Comm *comm)
!integer function GenericSend(sendbuf,sendcnt,recvcnt_,recvbuf_,mpi_type,blocksize,comm)
!    type(MPI_COMM)::comm
!    type(MPI_DATATYPE)::mpi_type
!    integer::blocksize
!    integer,dimension(:)::recvcnt_,sendcnt
!    real(dp),dimension(:,:)::sendbuf,recvbuf_
!    !local vars
!    integer :: np,my_rank,i,nreq_r,nreq_s,irank,ierr,istat
!    integer,dimension(:),allocatable::recvcnt
!    real(dp),dimension(:,:),allocatable::recvbuf
!    type(MPI_Request),dimension(:),allocatable::srequest,rrequest
!    type(MPI_Status), dimension(:),allocatable::statuses
!    real(dp):: a
!    nreq_r = 0
!    nreq_s = 0
!
!    call MPI_COMM_SIZE(comm,np,ierr)
!    call MPI_Comm_rank(comm,my_rank,ierr)
!  ! TODO: we should try to detect the type of mpi_type is set to zero!
!  ! Now that we have the send counts, we can notify the receivers what they
!  ! will get
!  recvcnt = GetReceiveCounts(sendcnt,my_rank,np,comm);
!  ! allocate memory for the receive buffers.  The blocksize parameter is used when each "unit" of
!  ! the data buffer is a multiple of the basic datatype (TheType)
!    allocate(recvbuf(np),STAT=istat)
!    do i = 0,np-1
!       allocate(&recvbuf(i),recvcnt(i)*blocksize);
!       ! initiate the sends and receives
!       if (np == 1) then
!          ! if we are not running in parallel, do a direct memory copy.
!          tmemcpy(recvbuf(0),sendbuf(0),sendcnt(0)*blocksize);
!       else
!          allocate(&statuses,np);
!          allocate(&rrequest,np);
!          allocate(&srequest,np);
!          do irank = 0; irank < np; irank++)
!             {
!             if (recvcnt(irank) > 0)then
!                call MPI_IRECV(recvbuf(irank),recvcnt(irank)*blocksize,mpi_type,irank,
!                irank,*comm,&rrequest(nreq_r++));
!             end if
!             if (sendcnt(irank) > 0)then
!                call MPI_ISEND(sendbuf(irank),sendcnt(irank)*blocksize,mpi_type,irank,
!                my_rank,*comm,&srequest(nreq_s++));
!	     end if
!   end if
!   !
!   ! wait for messages to complete
!   !
!   MPI_Waitall(nreq_r,rrequest,statuses);
!   MPI_Waitall(nreq_s,srequest,statuses);  
!   !
!   ! free up temp memory
!   !
!   deallocate(srequest);
!   deallocate(rrequest);
!   deallocate(statuses);
!end if
!
!    ! return info back to sender
!    if (recvcnt_ .ne. NULL) then
!          *recvcnt_ = recvcnt;
!      else 
!          deallocate(recvcnt);
!     end if
!  if (recvbuf_ .ne. NULL) then
!          *recvbuf_ = recvbuf;
!      else 
!          write(6,*)"ERROR in GenericSend: receive buffer is NULL!"
!      do i = 0,np-1
!         deallocate(recvbuf(i),STAT=istat)
!      end do
!      deallocate(recvbuf);
!     end if
! return(0);
!}




!  integer function Aswap(ptr,stride) 
!    TheType *ptrs(1) = {ptr}
!    integer strides(1) = {stride};
!    integer nptrs = 1
!    return Aswap(ptrs,strides,nptrs)
!  end function Aswap








  ! post the send messages for a specific set of pointers
  function PostSends_r(pcomm,ptrs,strides,nptrs) result(retval)
  use kinddefs, only : dp
    type(PComm_Type),intent(inout)::pcomm
    real(dp),intent(inout),pointer,dimension(:)::ptrs
    integer,intent(in),pointer,dimension(:)::strides
    integer,intent(in)::nptrs
    integer::retval
    !local vars
    integer :: unit,i,idomain,p,icv,k,istat,ierr
    integer ,pointer,dimension(:)::soffset
    real(dp),pointer,dimension(:)::dataptr
    real(dp):: a
    unit=0
    a=0
    nullify(soffset)
    
    ! this is a special bit of functionality that I think will be convenient,
    ! but we might want to change it up later.  Here, if we attempt to post sends
    ! while there are still receives and/or sends pending, we want to let those
    ! be fulfilled before we move on.
    !
    ! This is for efficiency, so that we can update some quantity by doing a 
    ! PostSends(), doing work, and then allowing the sends to finish later via
    ! yet another PostSends() or a FinishPendingComms().
    !
    ! The downside is that the user MUST remember that the quantity being updated
    ! will not fully be updated until the second PostSends() or FinishPendingComms()
    ! is encountered in the code.
    
    ierr = FinishPendingComms(pcomm)
    
    ! Timers switched off by FinishPendingComms() need to be restarted
    ! Therefore, they have been moved here
    
    !msg_timer->On();
    !send_timer->On();
    
    ! make sure that the send buffer is currently being unused
    
    !_ASSERT(send_buffer == NULL);
    !_ASSERT(pending_srequests == NULL);
    !_ASSERT(nreq_isend == 0);
    !_ASSERT(m_ptrs == NULL);
    !_ASSERT(m_nptrs == 0);
    !_ASSERT(m_strides == NULL);
    
    !CachePointers((void**)ptrs,strides,nptrs,GetMPIDatatype(a))
    
    allocate(pcomm%pending_srequests(pcomm%np),STAT=istat)
    
    ! determine the size of a unit (per cv amount of data)
    !do i = 0,nptrs-1 !note index starts at zero
    do i = 1,nptrs !note index starts at zero
       unit = unit + strides(i)
    end do

    !call AllocatePersistentBuffers(pcomm,unit*sizeof(TheType));!allocates send and receive buffers
    
    ! generate a small offset map which tells us, within a unit, where to
    ! place data
    allocate(soffset(nptrs),STAT=istat)
    do i = 1,nptrs-1 
       soffset(i) = soffset(i-1) + strides(i-1)
    end do
    
    !   ! allocate a buffer that will be undisturbed until the posted sends complete.
    
    !   !tpcalloc((TheType***)&send_buffer,np+1);
    !   allocate((TheType***)&send_buffer,np+1);
    do idomain = 1,pcomm%np
       if (pcomm%nsend(idomain) > 0) then
          !pack_timer->On()
          ! 	  ! allocate a buffer that we will copy data into
          ! 	  allocate((TheType**)&send_buffer(idomain),unit*nsend(idomain));
          ! first, we need to pack the data
          do i = 1,pcomm%nsend(idomain)
             !do p = 0,nptrs-1
             do p = 1,nptrs
                dataptr = ptrs(p)
                icv = pcomm%send(idomain)%array(i)
                !do k = 0,strides(p)-1
                do k = 1,strides(p)
                   !((TheType**)send_buffers)(idomain)(unit*(i-1) + soffset(p) + k) = dataptr(icv*strides(p) + k);
                   pcomm%send_buffers(idomain)%array(unit*(i-1) + soffset(p) + k) = &
                        dataptr(icv*strides(p) + k)
                end do
             end do
          end do
          
          ! dump the message that we are about to send out to file.
          !if (idebug == 2)then
          !   FILE *fp = NULL;
          !   char fname(1024);
          !   sprintf(fname,"uxmessage.%d_to_%d",my_rank,domain2rank(idomain));
          !   fp = fopen(fname,"w");
          !   if (fp)then
          !      fwrite(&nsend(idomain),sizeof(integer),1,fp);
          !      fwrite(&unit,sizeof(integer),1,fp);
          !      fwrite(send_buffers(idomain),sizeof(TheType),unit*nsend(idomain),fp);
          !      fclose(fp);
          !   end if
          !end if
          !pack_timer->Off();
          
          ! now, post a send
          if (idomain .ne. pcomm%my_subdomain)then
             !MPI_Isend(send_buffers(idomain),unit*nsend(idomain),GetMPIDatatype(a),domain2rank(idomain),my_subdomain,comm_world,&pending_srequests(nreq_isend)); 
             call MPI_ISEND(pcomm%send_buffers(idomain)%array,unit*pcomm%nsend(idomain), &
                            MPI_DOUBLE_PRECISION,pcomm%domain2rank(idomain), &
                            pcomm%my_subdomain,pcomm%comm_world, &
                            pcomm%pending_srequests(pcomm%nreq_isend))
             pcomm%nreq_isend=pcomm%nreq_isend+1
          end if
       end if
    end do
    deallocate(soffset,STAT=istat)
    !send_timer->Off();
    !msg_timer->Off();
    retval=0
    return
  end function PostSends_r



  ! post the send messages for a specific set of pointers
  function PostSends_i(pcomm,ptrs,strides,nptrs) result(retval)
    type(PComm_Type),intent(inout)::pcomm
    integer,intent(inout),pointer,dimension(:)::ptrs
    integer,intent(in),pointer,dimension(:)::strides
    integer,intent(in)::nptrs
    integer::retval
    !local vars
    integer :: unit,i,idomain,p,icv,k,istat,ierr
    integer ,pointer,dimension(:):: soffset
    integer,pointer,dimension(:)::dataptr
    integer:: a
    unit=0
    a=0
    nullify(soffset)
    
    ! this is a special bit of functionality that I think will be convenient,
    ! but we might want to change it up later.  Here, if we attempt to post sends
    ! while there are still receives and/or sends pending, we want to let those
    ! be fulfilled before we move on.
    !
    ! This is for efficiency, so that we can update some quantity by doing a 
    ! PostSends(), doing work, and then allowing the sends to finish later via
    ! yet another PostSends() or a FinishPendingComms().
    !
    ! The downside is that the user MUST remember that the quantity being updated
    ! will not fully be updated until the second PostSends() or FinishPendingComms()
    ! is encountered in the code.
    
    ierr = FinishPendingComms(pcomm)
   
    !CachePointers((void**)ptrs,strides,nptrs,GetMPIDatatype(a))
    
    allocate(pcomm%pending_srequests(pcomm%np),STAT=istat)
    
    ! determine the size of a unit (per cv amount of data)
    !do i = 0,nptrs-1 
    do i = 1,nptrs 
       unit = unit + strides(i)
    end do
    !AllocatePersistentBuffers(unit*sizeof(TheType));
    
    ! generate a small offset map which tells us, within a unit, where to
    ! place data
    allocate(soffset(nptrs),STAT=istat)
    !do i = 1,nptrs-1 
    do i = 2,nptrs 
       soffset(i) = soffset(i-1) + strides(i-1)
    end do
    
    !   ! allocate a buffer that will be undisturbed until the posted sends complete.
    do idomain = 1,pcomm%np
       if (pcomm%nsend(idomain) > 0) then
          ! 	  ! allocate a buffer that we will copy data into
          ! 	  allocate((TheType**)&send_buffer(idomain),unit*nsend(idomain));
          ! first, we need to pack the data
          do i = 1,pcomm%nsend(idomain)
             !do p = 0,nptrs-1
             do p = 1,nptrs
                dataptr = ptrs(p)
                icv = pcomm%send(idomain)%array(i)
                !do k = 0,strides(p)-1
                do k = 1,strides(p)
                   pcomm%send_buffers(idomain)%array(unit*(i-1) + soffset(p) + k) = &
                        dataptr(icv*strides(p) + k)
                end do
             end do
          end do
          

          ! now, post a send
          if (idomain .ne. pcomm%my_subdomain)then
             call MPI_ISEND(pcomm%send_buffers(idomain)%array,unit*pcomm%nsend(idomain), &
                            MPI_INTEGER,pcomm%domain2rank(idomain), &
                            pcomm%my_subdomain,pcomm%comm_world, &
                            pcomm%pending_srequests(pcomm%nreq_isend))
             pcomm%nreq_isend=pcomm%nreq_isend+1
          end if
       end if
    end do
    deallocate(soffset,STAT=istat)
    retval=0
    return
  end function PostSends_i


  
  ! post the receive messages for a specific set of pointers
  function PostReceives_r(pcomm,ptrs,strides,nptrs) result(retval)
  use kinddefs, only : dp
    type(PComm_Type),intent(inout)::pcomm
    real(dp),intent(inout),pointer,dimension(:)::ptrs
    integer,intent(in),pointer,dimension(:)::strides
    integer,intent(in)::nptrs
    integer::retval
    !local vars
    integer :: unit,i,idomain,ierr,istat
    real(dp)::a
    unit=0
    a=0
    !msg_timer->On();
    !recv_timer->On();
    
    ! make sure that no other requests are pending
    !_ASSERT(pending_rrequests == NULL);
    !  _ASSERT(recv_buffers == NULL);
    !_ASSERT(nreq_irecv == 0);
    
    ! make sure that data matches up
    !_ASSERT(m_ptrs .ne. NULL);
    !_ASSERT(m_nptrs > 0);
    !_ASSERT(m_strides .ne. 0);
    !_ASSERT(m_nptrs == nptrs);
    !do i = 0,m_nptrs-1
    !   _ASSERT(m_strides(i) == strides(i));
    !   _ASSERT(m_ptrs(i) == ptrs(i));
    !end do
    
    allocate(pcomm%pending_rrequests(pcomm%np),STAT=istat)
    !do i = 0,nptrs-1 
    do i = 1,nptrs 
       unit = unit + strides(i)
    end do
    !AllocatePersistentBuffers(unit*sizeof(TheType));
    !  allocate((TheType***)&recv_buffer,np + 1);
    do idomain = 1,pcomm%np
       if (pcomm%nexpect(idomain) > 0)then
          ! allocate((TheType**)&expect_buffer(idomain),nexpect(idomain)*unit);
          if (idomain .ne. pcomm%my_subdomain)then
             ! post a receive
             call MPI_IRECV(pcomm%recv_buffers(idomain)%array,&
                  unit*pcomm%nexpect(idomain),&
!                  GetMPIDatatype(a),&
                  MPI_DOUBLE_PRECISION,&
                  pcomm%domain2rank(idomain),&
                  idomain,pcomm%comm_world,&
                  pcomm%pending_rrequests(pcomm%nreq_irecv)) 
             pcomm%nreq_irecv=pcomm%nreq_irecv+1
          end if
       end if
    end do
    !recv_timer->Off();
    !msg_timer->Off();
    retval=0
    return
  end function PostReceives_r


  
  ! post the receive messages for a specific set of pointers
  function PostReceives_i(pcomm,ptrs,strides,nptrs) result(retval)
    type(PComm_Type),intent(inout)::pcomm
    integer,intent(inout),pointer,dimension(:)::ptrs
    integer,intent(in),pointer,dimension(:)::strides
    integer,intent(in)::nptrs
    integer::retval
    !local vars
    integer :: unit,i,idomain,ierr,istat
    integer::a
    unit=0
    a=0
    
    allocate(pcomm%pending_rrequests(pcomm%np),STAT=istat)
    !do i = 0,nptrs-1 
    do i = 1,nptrs 
       unit = unit + strides(i)
    end do
    !AllocatePersistentBuffers(unit*sizeof(TheType));
    !  allocate((TheType***)&recv_buffer,np + 1);
    do idomain = 1,pcomm%np
       if (pcomm%nexpect(idomain) > 0)then
          ! allocate((TheType**)&expect_buffer(idomain),nexpect(idomain)*unit);
          if (idomain .ne. pcomm%my_subdomain)then
             ! post a receive
             call MPI_IRECV(pcomm%recv_buffers(idomain)%array,&
                  unit*pcomm%nexpect(idomain),&
                  MPI_INTEGER,&
                  pcomm%domain2rank(idomain),&
                  idomain,pcomm%comm_world,&
                  pcomm%pending_rrequests(pcomm%nreq_irecv)) 
             pcomm%nreq_irecv=pcomm%nreq_irecv+1
          end if
       end if
    end do
    retval=0
    return
  end function PostReceives_i
  



 
! subroutine CachePointers_r(pcomm,ptrs,strides,nptrs)
!  use kinddefs, only : dp
!   type(PComm_Type),intent(inout)::pcomm
!   integer,intent(in)::nptrs
!   type(realarray),intent(inout),dimension(:)::ptrs
!   integer,intent(in),dimension(:)::strides
!   integer ::i,istat
!   if (pcomm%m_nptrs == 0)then
!      allocate(pcomm%m_ptrs_r(nptrs),STAT=istat)
!      allocate(pcomm%m_strides(nptrs),STAT=istat)
!      pcomm%m_nptrs = nptrs
!      !do i = 0,nptrs-1
!      do i = 1,nptrs
!         pcomm%m_ptrs_r(i) => ptrs(i)
!         pcomm%m_strides(i) = strides(i)
!      end do
!      pcomm%m_cached_datatype = MPI_DOUBLE_PRECISION
!   else
!!      !_ASSERT(pcomm%m_nptrs == nptrs);
!   end if
!   return
! end subroutine CachePointers_r
 
! subroutine CachePointers_i(pcomm,ptrs,strides,nptrs)
!   type(PComm_Type),intent(inout)::pcomm
!   integer,intent(in)::nptrs
!   type(intarray),intent(inout),pointer,dimension(:)::ptrs
!   integer,intent(in),dimension(:)::strides
!   integer ::i,istat
!   if (pcomm%m_nptrs == 0)then
!      allocate(pcomm%m_ptrs_i(nptrs),STAT=istat)
!      allocate(pcomm%m_strides(nptrs),STAT=istat)
!      pcomm%m_nptrs = nptrs
!      !do i = 0,nptrs-1
!      do i = 1,nptrs
!         pcomm%m_ptrs_i(i) => ptrs(i)
!         pcomm%m_strides(i) = strides(i)
!      end do
!      pcomm%m_cached_datatype = MPI_INTEGER
!   else
!!      !_ASSERT(pcomm%m_nptrs == nptrs);
!   end if
!   return
! end subroutine CachePointers_i
 


 function PostSendsAndReceives_r(pcomm,ptrs,strides,nptrs) result(retval)
  use kinddefs, only : dp
   type(PComm_Type),intent(inout)::pcomm
   real(dp),intent(inout),pointer,dimension(:)::ptrs
   integer,intent(in),pointer,dimension(:)::strides
   integer,intent(in)::nptrs
   integer::retval,ierr
   ierr = PostSends_r(pcomm,ptrs,strides,nptrs)
   ierr = PostReceives_r(pcomm,ptrs,strides,nptrs)
   ierr = FinishPendingComms(pcomm)
   retval=0
   return
 end function PostSendsAndReceives_r
 
 function PostSendsAndReceives_i(pcomm,ptrs,strides,nptrs) result(retval)
   type(PComm_Type),intent(inout)::pcomm
   integer,intent(inout),pointer,dimension(:)::ptrs
   integer,intent(in),pointer,dimension(:)::strides
   integer,intent(in)::nptrs
   integer::retval,ierr
   ierr = PostSends_i(pcomm,ptrs,strides,nptrs)
   ierr = PostReceives_i(pcomm,ptrs,strides,nptrs)
   ierr = FinishPendingComms(pcomm)
   retval=0
   return
 end function PostSendsAndReceives_i


  function Aswap_r(pcomm,ptrs,strides,nptrs) result(retval)
  use kinddefs, only : dp
    type(PComm_Type),intent(inout)::pcomm
    real(dp),intent(inout),pointer,dimension(:)::ptrs
    integer,intent(in),pointer,dimension(:)::strides
    integer,intent(in)::nptrs
    integer::ierr,retval
    ierr = PostSendsAndReceives_r(pcomm,ptrs,strides,nptrs)
    ierr = FinishPendingComms(pcomm)
    retval=0
    return
  end function Aswap_r

  function Aswap_i(pcomm,ptrs,strides,nptrs) result(retval)
    type(PComm_Type),intent(inout)::pcomm
    integer,intent(inout),pointer,dimension(:)::ptrs
    integer,intent(in),pointer,dimension(:)::strides
    integer,intent(in)::nptrs
    integer::ierr,retval
    ierr = PostSendsAndReceives_i(pcomm,ptrs,strides,nptrs)
    ierr = FinishPendingComms(pcomm)
    retval=0
    return
  end function Aswap_i



 function GetMaxScalar_r(pcomm,local_val,global_val) result(retval)
  use kinddefs, only : dp
   type(PComm_Type),intent(inout)::pcomm
   real(dp),intent(in) ::local_val
   real(dp),intent(out)::global_val
   integer::retval
   !local vars
   integer:: ierr
   integer::op!type(MPI_OP)
   integer::DataType!type(MPI_DATA_TYPE)
   if(IsParallel(pcomm).eq.0) then
      !call MPI_Op_create(DUAL_MPI_MAX,1,op,ierr) 
      call MPI_Op_create(MPI_MAX,1,op,ierr) 
      !DataType = GetMPIDatatype(local_val)
      !call MPI_Allreduce(local_val,global_val,1,DataType,op,comm_world,ierr)
      call MPI_Allreduce(local_val,global_val,1,MPI_DOUBLE_PRECISION,op,&
                         pcomm%comm_world,ierr)
      call MPI_Op_free(op)
   else
      global_val = local_val
   end if
   retval=0
   return
 end function GetMaxScalar_r


 function GetMaxScalar_i(pcomm,local_val,global_val) result(retval)
   type(PComm_Type),intent(inout)::pcomm
   integer,intent(in) ::local_val
   integer,intent(out)::global_val
   integer::retval
   !local vars
   integer:: ierr
   integer::op!type(MPI_OP)
   integer::DataType!type(MPI_DATA_TYPE)
   if(IsParallel(pcomm).eq.1) then
      !call MPI_Op_create(DUAL_MPI_MAX,1,op,ierr) 
      call MPI_Op_create(MPI_MAX,1,op,ierr) 
      !DataType = GetMPIDatatype(local_val)
      !call MPI_Allreduce(local_val,global_val,1,DataType,op,comm_world,ierr)
      call MPI_Allreduce(local_val,global_val,1,MPI_INTEGER,op,&
                         pcomm%comm_world,ierr)
      call MPI_Op_free(op)
   else
      global_val = local_val
   end if
   retval=0
   return
 end function GetMaxScalar_i


 function GetMinScalar_r(pcomm,local_val, global_val) result(retval)
  use kinddefs, only : dp
   type(PComm_Type),intent(inout)::pcomm
   real(dp),intent(in) ::local_val
   real(dp),intent(out)::global_val
   integer::retval
   !local vars
   integer:: ierr
   integer::op!type(MPI_OP)
   integer::DataType!type(MPI_DATA_TYPE)
   if(IsParallel(pcomm).eq.1)then
      !call MPI_Op_create(DUAL_MPI_MIN,1,op); 
      call MPI_Op_create(MPI_MIN,1,op,ierr) 
      !DataType = GetMPIDatatype(local_val);
      !call MPI_Allreduce(local_val,global_val,1,DataType,op,pcomm%comm_world);
      call MPI_Allreduce(local_val,global_val,1,MPI_DOUBLE_PRECISION,op,&
                         pcomm%comm_world,ierr)
      call MPI_Op_free(op) 
   else
      global_val = local_val
   end if
   retval=0
   return
 end function GetMinScalar_r

 function GetMinScalar_i(pcomm,local_val, global_val) result(retval)
   type(PComm_Type),intent(inout)::pcomm
   integer,intent(in) ::local_val
   integer,intent(out)::global_val
   integer::retval
   !local vars
   integer:: ierr
   !type(MPI_OP)::op
   !type(MPI_DATA_TYPE)::DataType
   if(IsParallel(pcomm).eq.1)then
      !call MPI_Op_create(DUAL_MPI_MIN,1,op); 
      !call MPI_Op_create(MPI_MIN,1,op,ierr) 
      !DataType = GetMPIDatatype(local_val);
      !call MPI_Allreduce(local_val,global_val,1,DataType,op,pcomm%comm_world);
      call MPI_Allreduce(local_val,global_val,1,MPI_INTEGER,MPI_MIN,&
                         pcomm%comm_world,ierr)
      !MPI_Op_free(op) 
   else
      global_val = local_val
   end if
   retval=0
   return
 end function GetMinScalar_i
 

 function GetScalarSum_r(pcomm,local_val,global_val) result(retval)
   use kinddefs, only : dp
  type(PComm_Type),intent(inout)::pcomm
   real(dp),intent(in) ::local_val
   real(dp),intent(out) ::global_val
   integer::retval
   !local vars
   integer:: j,ierr
   integer::op!type(MPI_OP)
   integer::DataType!type(MPI_DATA_TYPE)
   if(IsParallel(pcomm).eq.1)then
      call MPI_Allreduce(local_val,global_val,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
                         pcomm%comm_world)
   else
      global_val = local_val
   end if
   retval=0
   return
 end function GetScalarSum_r

 function GetScalarSum_i(pcomm,local_val,global_val) result(retval)
   type(PComm_Type),intent(inout)::pcomm
   integer,intent(in) ::local_val
   integer,intent(out) ::global_val
   integer::retval
   !local vars
   integer:: j,ierr
   integer::op!type(MPI_OP)
   integer::DataType!type(MPI_DATA_TYPE)
   if(IsParallel(pcomm).eq.1)then
   call MPI_Allreduce(local_val,global_val,1,MPI_INTEGER,MPI_SUM,&
        pcomm%comm_world)
   else
      global_val = local_val
   end if
   retval=0
   return
 end function GetScalarSum_i

 function GetVectorSum_r(pcomm,local_vec,global_vec,vector_size) result(retval)
   use kinddefs, only : dp
  type(PComm_Type),intent(inout)::pcomm
   real(dp),intent(in) ,dimension(:)::local_vec
   real(dp),intent(out),dimension(:)::global_vec
   integer,intent(in) :: vector_size
   !local vars
   integer:: ierr,i,j,retval
   integer::op!type(MPI_OP)
   integer::DataType!type(MPI_DATA_TYPE)
   if (vector_size == 0) then
      ! nothing to do.
      retval=0
      return
   end if
   if (vector_size < 0) then
      write(6,*)"FATAL: invalid vector size passed to GetVectorSum().not. vector_size =",vector_size
      retval=-1
      return
   end if
   if(IsParallel(pcomm).eq.1) then
      !*local_vector = NULL;
      !allocate(local_vector, (DUAL_SIZE+1)*vector_size);
      !*global_vector = NULL;
      !allocate(&global_vector, (DUAL_SIZE+1)*vector_size);	
      !do i = 0; i < vector_size; i++)then
      !{
      !local_vector((DUAL_SIZE+1)*i + 0) = local_vec(i).a(0);
      !do(j=0;j<=DUAL_SIZE;++j)local_vector((DUAL_SIZE+1)*i + j) = local_vec(i).a(j);
      !}
      !call MPI_Allreduce(local_vector,global_vector,vector_size,MPI_DOUBLE_PRECISION,MPI_SUM,pcomm%comm_world,ierr)
      call MPI_Allreduce(local_vec,global_vec,vector_size,&
                         MPI_DOUBLE_PRECISION,MPI_SUM,pcomm%comm_world,ierr)
      !do (i = 0; i < vector_size; i++)
      !do j=0,DUAL_SIZE;++j)
      !   global_vec(i).a(j) = global_vector(i*(DUAL_SIZE+1)+j);
      !end do
      !end do
      !deallocate(local_vector);
      !deallocate(global_vector);
   else
      !do i = 0,vector_size-1
      do i = 1,vector_size
         global_vec(i) = local_vec(i)
      end do
   end if
   retval=0
   return
 end function GetVectorSum_r


 function GetVectorSum_i(pcomm,local_vec,global_vec,vector_size) result(retval)
   type(PComm_Type),intent(inout)::pcomm
   integer,intent(in) ,dimension(:)::local_vec
   integer,intent(out),dimension(:)::global_vec
   integer,intent(in) :: vector_size
   !local vars
   integer:: ierr,i,j,retval
   integer::op!type(MPI_OP)
   integer::DataType!type(MPI_DATA_TYPE)
   if (vector_size == 0) then
      ! nothing to do.
      retval=0
      return
   end if
   if (vector_size < 0) then
      write(6,*)"FATAL: invalid vector size passed to GetVectorSum().not. vector_size =",vector_size
      retval=-1
      return
   end if
   if(IsParallel(pcomm).eq.1) then
      !*local_vector = NULL;
      !allocate(local_vector, (DUAL_SIZE+1)*vector_size);
      !*global_vector = NULL;
      !allocate(&global_vector, (DUAL_SIZE+1)*vector_size);	
      !do i = 0; i < vector_size; i++)then
      !{
      !local_vector((DUAL_SIZE+1)*i + 0) = local_vec(i).a(0);
      !do(j=0;j<=DUAL_SIZE;++j)local_vector((DUAL_SIZE+1)*i + j) = local_vec(i).a(j);
      !}
      !call MPI_Allreduce(local_vector,global_vector,vector_size,MPI_DOUBLE_PRECISION,MPI_SUM,pcomm%comm_world,ierr)
      call MPI_Allreduce(local_vec,global_vec,vector_size,&
                         MPI_INTEGER,MPI_SUM,pcomm%comm_world,ierr)
      !do (i = 0; i < vector_size; i++)
      !do j=0,DUAL_SIZE;++j)
      !   global_vec(i).a(j) = global_vector(i*(DUAL_SIZE+1)+j);
      !end do
      !end do
      !deallocate(local_vector);
      !deallocate(global_vector);
   else
      do i = 1,vector_size
         global_vec(i) = local_vec(i)
      end do
   end if
   retval=0
   return
 end function GetVectorSum_i




!===========================================================================================
end module pcomm

