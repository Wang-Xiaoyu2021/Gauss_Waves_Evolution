program test
      implicit none
      integer(8) :: x
      real(8)    :: Pi=3.1415926
      complex(8) :: cj=(0.,1.), cjn
      do x=1, 10
            cjn = (0.2**(0.5))*(1/(Pi**(0.25)))*Exp(cj*0.5*x)*Exp((-x**2*0.04)/2)
            write(*,*)cjn
      end do
end program
