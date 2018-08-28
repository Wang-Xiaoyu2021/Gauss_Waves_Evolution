program main
      implicit none
      integer(kind=4)   :: a
      character(9)      :: c1

      do a = 1, 19
      write(c1,'(i3)')a
      enddo
      write(*,*)c1
end program main
