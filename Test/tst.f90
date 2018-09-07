program Test
!本程序旨在证明我想要做的scale步骤是对的，我用的maxval函数的用法也是对的，mmp
    implicit none
    real,allocatable  :: a(:,:), b(:,:)
    integer              :: i, j
    allocate(a(10,10))
    do i=1, 10
        do j=1, 10
            a(i,j)=i*j
        end do
    end do
    !write(*,*)a
    do j=1,10
        do i=1,10
            a(i,j)=a(i,j)/maxval(a(1:10,j))
        end do
    end do
    open(1,file="out.dat")
        do i=1,10
            do j=1, 10
                write(1,'(f4.1)',advance='no')a(i,j)
            end do
            write(1,*)''
        end do
    close(1)

end program

