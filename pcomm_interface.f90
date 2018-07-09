!this module contains interfaces for pcomm function to allow
! one function name for multiple datatypes
module pcomm_interface
  use kinddefs
  use pcomm



 interface PostSends
    function PostSends_r(pcomm,ptrs,strides,nptrs) result(retval)
      use kinddefs, only : dp
      type(PComm_Type),intent(inout)::pcomm
      real(dp),intent(inout),pointer,dimension(:)::ptrs
      integer,intent(in),pointer,dimension(:)::strides
      integer,intent(in)::nptrs
      integer::retval
    end function PostSends_r
    function PostSends_i(pcomm,ptrs,strides,nptrs) result(retval)
      type(PComm_Type),intent(inout)::pcomm
      integer,intent(inout),pointer,dimension(:)::ptrs
      integer,intent(in),pointer,dimension(:)::strides
      integer,intent(in)::nptrs
      integer::retval
    end function PostSends_i
 end interface

 interface PostReceives
  function PostReceives_r(pcomm,ptrs,strides,nptrs) result(retval)
      use kinddefs, only : dp
      type(PComm_Type),intent(inout)::pcomm
      real(dp),intent(inout),pointer,dimension(:)::ptrs
      integer,intent(in),pointer,dimension(:)::strides
      integer,intent(in)::nptrs
      integer::retval
    end function PostReceives_r
  function PostReceives_i(pcomm,ptrs,strides,nptrs) result(retval)
      type(PComm_Type),intent(inout)::pcomm
      integer,intent(inout),pointer,dimension(:)::ptrs
      integer,intent(in),pointer,dimension(:)::strides
      integer,intent(in)::nptrs
      integer::retval
    end function PostReceives_i
 end interface


! interface CachePointers
!    subroutine CachePointers_r(pcomm,ptrs,strides,nptrs)
!      use kinddefs, only : dp
!      type(PComm_Type),intent(inout)::pcomm
!      real(dp),intent(inout),pointer,dimension(:)::ptrs
!      integer,intent(in),pointer,dimension(:)::strides
!      integer,intent(in)::nptrs
!    end subroutine CachePointers_r
!    subroutine CachePointers_i(pcomm,ptrs,strides,nptrs)
!      type(PComm_Type),intent(inout)::pcomm
!      integer,intent(inout),pointer,dimension(:)::ptrs
!      integer,intent(in),pointer,dimension(:)::strides
!      integer,intent(in)::nptrs
!    end subroutine CachePointers_i
! end interface

 interface PostSendsAndReceives
    function PostSendsAndReceives_r(pcomm,ptrs,strides,nptrs) result(retval)
      use kinddefs, only : dp
      type(PComm_Type),intent(inout)::pcomm
      real(dp),intent(inout),pointer,dimension(:)::ptrs
      integer,intent(in),pointer,dimension(:)::strides
      integer,intent(in)::nptrs
      integer::retval
    end function PostSendsAndReceives_r
    function PostSendsAndReceives_i(pcomm,ptrs,strides,nptrs) result(retval)
      type(PComm_Type),intent(inout)::pcomm
      integer,intent(inout),pointer,dimension(:)::ptrs
      integer,intent(in),pointer,dimension(:)::strides
      integer,intent(in)::nptrs
      integer::retval
    end function PostSendsAndReceives_i
 end interface


 interface Aswap
  function Aswap_r(pcomm,ptrs,strides,nptrs) result(retval)
      use kinddefs, only : dp
      type(PComm_Type),intent(inout)::pcomm
      real(dp),intent(inout),pointer,dimension(:)::ptrs
      integer,intent(in),pointer,dimension(:)::strides
      integer,intent(in)::nptrs
      integer::retval
    end function Aswap_r
  function Aswap_i(pcomm,ptrs,strides,nptrs) result(retval)
      type(PComm_Type),intent(inout)::pcomm
      integer,intent(inout),pointer,dimension(:)::ptrs
      integer,intent(in),pointer,dimension(:)::strides
      integer,intent(in)::nptrs
      integer::retval
    end function Aswap_i
 end interface

 interface GetMaxScalar
    function GetMaxScalar_r(pcomm,local_val,global_val) result(retval)
      use kinddefs, only : dp
      type(PComm_Type),intent(inout)::pcomm
      real(dp),intent(in) ::local_val
      real(dp),intent(out)::global_val
      integer::retval
    end function GetMaxScalar_r
    function GetMaxScalar_i(pcomm,local_val,global_val) result(retval)
      type(PComm_Type),intent(inout)::pcomm
      integer,intent(in) ::local_val
      integer,intent(out)::global_val
      integer::retval
    end function GetMaxScalar_i
 end interface

 interface GetMinScalar
    function GetMinScalar_r(pcomm,local_val,global_val) result(retval)
      use kinddefs, only : dp
      type(PComm_Type),intent(inout)::pcomm
      real(dp),intent(in) ::local_val
      real(dp),intent(out)::global_val
      integer::retval
    end function GetMinScalar_r
    function GetMinScalar_i(pcomm,local_val,global_val) result(retval)
      type(PComm_Type),intent(inout)::pcomm
      integer,intent(in) ::local_val
      integer,intent(out)::global_val
      integer::retval
    end function GetMinScalar_i
 end interface

 interface GetScalarSum
    function GetScalarSum_r(pcomm,local_val,global_val) result(retval)
      use kinddefs, only : dp
      type(PComm_Type),intent(inout)::pcomm
      real(dp),intent(in) ::local_val
      real(dp),intent(out)::global_val
      integer::retval
    end function GetScalarSum_r
    function GetScalarSum_i(pcomm,local_val,global_val) result(retval)
      type(PComm_Type),intent(inout)::pcomm
      integer,intent(in) ::local_val
      integer,intent(out)::global_val
      integer::retval
    end function GetScalarSum_i
 end interface

 interface GetVectorSum
    function GetVectorSum_r(pcomm,local_vec,global_vec,vector_size) result(retval)
      use kinddefs, only : dp
      type(PComm_Type),intent(inout)::pcomm
      real(dp),intent(in) ,dimension(:)::local_vec
      real(dp),intent(out),dimension(:)::global_vec
      integer::retval
    end function GetVectorSum_r
    function GetVectorSum_i(pcomm,local_vec,global_vec,vector_size) result(retval)
      type(PComm_Type),intent(inout)::pcomm
      integer,intent(in) ,dimension(:)::local_vec
      integer,intent(out),dimension(:)::global_vec
      integer::retval
    end function GetVectorSum_i
 end interface

contains

end module pcomm_interface
