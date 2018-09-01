program gauss
    implicit none
    !------------------------------------------
    real                      :: Pi = 3.1415926
    integer,parameter         :: DP = 8
    complex(DP),allocatable   :: psi(:,:)
    complex(DP),allocatable   :: phi(:,:)
    complex                   :: cj=(0.,1.), A
    complex(DP),allocatable   :: matrx(:,:), matrx2(:,:), B(:), V(:)
    integer(DP)               :: npts=100, i, j
    complex(DP),external      :: begin
    real(DP)                  :: x0=-50, xf=50, x, deltax, t, deltat=1
    deltax=(xf-x0)/npts
    !------------------------------------------
    allocate(psi(0:npts,0:0))
    x=x0
    do i=0, npts
        psi(i,0) = begin(x,Pi,cj)
        x=x+1
    end do
    !------------------------------------------
    !do i=0, npts
    !      psi(i)=abs(psi(i))
    !end do
    !write(*,*)psi
    !------------------------------------------
    allocate(matrx(0:npts,0:100))
    matrx(:,:)=0
    !write(*,*)matrx
    A=cj*0.25*deltat/(deltax**2)
    do i=0, npts-1
        j=i+1
        matrx(i,j)=A
    end do
    !do i=0, 5
    !      write(*,*)matrx(i,0:5)
    !end do
    do i=1, npts
        j=i-1
        matrx(i,j)=A
    end do 
    !do i=0, 5
    !      write(*,*)matrx(i,0:5)
    !end do
    allocate(B(0:npts))
    allocate(V(0:100))
    V(:)=0
    V(10:15)=1
    V(85:90)=1
    do i=0, npts
        B(i)=cj*0.5*deltat*V(i)
    end do
    !write(*,*)B
    do i=0, npts
        j=i
        matrx(i,j)=1-2*A-B(i)
    end do
    !do i=0, 5 
    !      write(*,*)matrx(i,0:5)
    !end do
    !-------------------------------------------
    allocate(matrx2(0:npts,0:100))
    matrx2(:,:)=0
    do i=0, npts-1
        j=i+1
        matrx2(i,j)=-A
    end do
    do i=1, npts
        j=i-1
        matrx2(i,j)=-A
    end do 
    !write(*,*) B
    do i=0, npts
        j=i
        matrx2(i,j)=1+2*A+B(i)
    end do
    !do i=0,5
    !      write(*,*)matrx2(i,0:5)
    !end do
    !-------------------------------------------
    allocate(phi(0:npts,0:0))
    phi=matmul(matrx,psi) 
    !write(*,*)phi
    














end program


function begin(x,Pi,cj)
    !函数在定义时只是再次说明一下变量类型，不可以再次赋值
    implicit none
    real(8)           ::  x
    complex(8)        ::  begin
    real              ::  Pi
    complex           ::  cj 
    begin=(0.2**(0.5))*(1/(Pi**(0.25)))*Exp(cj*0.5*x)*Exp((-x**2*0.04)/2)
end function


