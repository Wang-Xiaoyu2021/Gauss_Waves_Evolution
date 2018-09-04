program gauss
    implicit none
    !------------------------------------------
    real                      :: Pi = 3.1415926
    integer,parameter         :: DP = 8
    complex(DP),allocatable   :: psi(:,:), psinx(:,:), output(:,:)
    complex(DP),allocatable   :: phi(:,:), temp(:,:)
    complex                   :: cj=(0.,1.), A
    complex(DP),allocatable   :: matrx(:,:), matrx2(:,:),inv_matrx2(:,:), B(:), V(:)
    integer(DP)               :: npts=100, i, j
    complex(DP),external      :: begin
    real(DP)                  :: x0=-50, xf=50, x, deltax, t, deltat=2
    real(DP),allocatable      :: realoutput(:,:)
    character(DP)             :: temp1
    deltax=(xf-x0)/npts
    !------------------------------------------
    !定义初始时刻的波包形状
    allocate(psi(0:npts,0:0))
    allocate(temp(0:npts,0:0))
    x=x0
    do i=0, npts
        psi(i,0) = begin(x,Pi,cj)
        x=x+1
    end do
    temp=psi
    !------------------------------------------
    !do i=0, npts
    !      psi(i)=abs(psi(i))
    !end do
    !write(*,*)psi
    !------------------------------------------
    !定义第一个矩阵，用于计算phi的，在这期间定义了势能形状
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
    !定义第二个矩阵，用于求逆矩阵进而求解下一时刻的psi
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
    !--------------------------------------------
    !求解phi
    allocate(output(0:npts,0:npts))
    allocate(phi(0:npts,0:0))
    allocate(inv_matrx2(0:npts,0:npts))
    allocate(psinx(0:npts,0:0))
    do i=1,100
        phi=matmul(matrx,psi) 
        !write(*,*)phi
    !--------------------------------------------
        !print*,"test"
        call inv_mat(matrx2,inv_matrx2)
        !print*,"test2"
        inv_matrx2=matmul(matrx2,inv_matrx2)
        !do i=0, 5 
        !      write(*,*)inv_matrx2(i,0:5)
        !end do
        psinx=matmul(inv_matrx2,phi)
        !do i=0,5
        !    write(*,*)psinx(i,0:5)
        !end do 
        output(:,i:i)=psinx
        psi=psinx
    end do
    output(0:npts,0:0)=temp
    output=abs(output)
    allocate(realoutput(0:npts,0:npts))
    realoutput=real(output)**2
    do i=0,100
        realoutput(i,0:0)=i*deltax
    end do
    !write(*,*) realoutput
    open(1,file="realoutput.dat")
    write(1,'(a15)',advance='no')'x'
    do i=1,100
        write(temp1,'(a,i0)') 't',i
        write(1,'(a15)',advance='no')trim(temp1)
    end do
    write(1,*)'' 
    do i=0,100
        do j=0,100
            write(1,'(f15.7)',advance='no')realoutput(i:i,j:j)
        end do
            write(1,*)''
    end do
    close(1)
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


