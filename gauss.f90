program gauss
      implicit none
      !------------------------------------------
      real                      :: Pi = 3.1415926
      integer,parameter         :: DP = 8
      complex(DP),allocatable   :: psi(:)
      complex(DP),allocatable   :: phi(:)
      complex                   :: cj=(0.,1.)
      integer(DP)               :: npts=100, i
      complex(DP),external             :: begin
      real(DP)                  :: x0=-50, xf=50, x
      !------------------------------------------
      allocate(psi(0:npts))
      x=x0
      do i=0, npts
            psi(i) = begin(x,Pi,cj)
            x=x+1
      end do
      !do i=0, npts
      !      psi(i)=abs(psi(i))
      !end do
      write(*,*)psi














end program


function begin(x,Pi,cj)
      implicit none
      real(8)           ::  x
      complex(8)        ::  begin
      real              ::  Pi
      complex           ::  cj 
      begin=(0.2**(0.5))*(1/(Pi**(0.25)))*Exp(cj*0.5*x)*Exp((-x**2*0.04)/2)
end function

