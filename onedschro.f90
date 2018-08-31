	! program for calculating a 1-dimensional wave function
	! atomic units used throughout
	! everything on an equally spaced grid x_j=x_0+dx*j
	! http://www.physics.purdue.edu/~robichf/qmmovies/theory.htm
	!
program onedschreq
	implicit none
	!----------------------------------------------------------------------------
	integer    , parameter   :: DP=8
	integer                  :: fid, ilog, ipsi, ipsit, irecl
	integer                  :: iout, it, j, igo, ites
	integer                  :: npts,nt,nmod, nptfin
	real(DP)                 :: ak,x0,xf,xm,wi,dt
	real(DP)                 :: dx, tfin, t, tes, pref, xnorm, s, fsum
	real(DP)                 :: x, dx2m1, potential
	real(DP)   , allocatable :: xt(:,:)
	real(DP)   , allocatable :: vpot(:)
	complex(DP)              :: offd, II, c
	complex(DP), allocatable :: psi(:)
	complex(DP), allocatable :: psit(:,:)
	complex(DP), allocatable :: phi(:)
	complex(DP), allocatable :: phit(:,:)
	complex(DP), allocatable :: diag(:)
	character(10)            :: ctime
	real                     :: time
	real                     :: btime(10)
	real                     :: etime(10)
	real                     :: dtime(10)
	character(200)           :: finp, flog, fout, fpsi, fpsit, ctmp
	logical                  :: ifout
	!----------------------------------------------------------------------------
	namelist /control/ ak, x0, xf, xm, wi, npts, dt, nt, &
		nmod, flog, nptfin, fout
	!----------------------------------------------------------------------------
	! vpot will contain the potential + diagonal term from kinetic energy
	! vpot 将包含来自动能的潜在+对角线项
	! psi  will contain the wave function
	! psi  将包含波函数
	! phi  will contain the inhomogeneous term
	! phi  将包含非齐次项
	! diag will contain the diagonal term of the matrix
	! diag 将包含矩阵的对角线项
	! offd will contain the offdiagonal term of the matrix, since
	!      this does not depend on the index for this problem it
	!      has been changed to a number
	! offd 将包含矩阵的非对角线项，因为这个并不依赖于这个问题相关的参数，所以被定义为数
	! the matrix is assumed to be symmetric
	! 假设矩阵是对称的？（这矩阵本来就是对称的呀）
	!----------------------------------------------------------------------------
	! ak   = initial momentum   初始动量
	! x0   = minimum of x   x的最小值
	! xf   = maximum of x   x的最大值
	! xm   = middle of packet   包的中间
	! wi   = width of initial wave function   初始波函数的宽度
	! npts = number of spatial points   空间点数
	! dt   = time step size   时间步长
	! nt   = number of time steps to do propagation   传播的时间步数
	! nmod = how many steps to skip before printing out wave function 
	! 打印波函数之前需要进行的步数
	!----------------------------------------------------------------------------
	call cpu_time(time)
	btime(1) = time
	!应该是用于测试代码段执行的时间的
	!----------------------------------------------------------------------------
	fid   = 25
	ilog  = 26
	ipsi  = 27
	ipsit = 28
	ifout = .false.
	II    = cmplx(0.D0, 1.D0, DP)
	! 赋值复数双精度i
	fpsi  = "psi.dat"
	fpsit = "psit.dat"
	! 给字符串变量赋文件名称
	!----------------------------------------------------------------------------
	! Example of input file:
	! &control
	!   ak = 0
	!   x0 = 0
	!   xf = 10
	!   xm = 5
	!   wi =  npts, dt, nt, &
	!   nmod, flog, nptfin, fout
	! \
	!----------------------------------------------------------------------------
	call getarg(1, finp)
	open(fid,file=finp,status='old')
		read(fid,control)
	close(fid)
	!----------------------------------------------------------------------------
	if(npts .ge. nptfin) stop 'increase the size of nptfin'
	!----------------------------------------------------------------------------
	dx   = (xf-x0)/dfloat(npts+1)
	tfin = dt*dfloat(nt)
	iout = nt/nmod
	!----------------------------------------------------------------------------
	allocate(vpot(  nptfin   ))
	allocate(psi (0:nptfin   ))
	allocate(psit(0:nptfin,nt))
	allocate(phi (0:nptfin   ))
	allocate(phit(0:nptfin,nt))
	allocate(xt  (0:nptfin,nt))
	allocate(diag(0:nptfin   ))
	!----------------------------------------------------------------------------
	open(ilog, file=flog)
	!----------------------------------------------------------------------------
	! initialize the wave function
	!
	!----------------------------------------------------------------------------
	write(ilog,'(2x, "dx   = ", f15.5)') dx
	write(ilog,'(2x, "tfin = ", f15.5)') tfin
	write(ilog,'(2x, "number of output wave functions = ", i15)') iout+1
	!----------------------------------------------------------------------------
	dx2m1       = 1.d0/(dx*dx)
	fsum        = 0.d0
	psi(0)      = (0.d0, 0.d0)
	psi(npts+1) = (0.d0, 0.d0)
	!----------------------------------------------------------------------------
	do j=1, npts
		x=x0+dx*j
		s=(x-xm)/wi
		psi(j)=exp(-s*s+II*ak*x)
		!-------------------------------------------------------------------------
		! if the potential is time dependent need to put this in the 4 loop
		!
		potential = -x/100.D0
		vpot(j)   = potential+dx2m1
		fsum      = fsum + abs(psi(j))**2
	end do
	!----------------------------------------------------------------------------
	! normalize the wave function
	!
	xnorm=1.d0/dsqrt(fsum)
	psi=psi*xnorm
	!----------------------------------------------------------------------------
	! output the initial wave function
	!
	if(ifout) then
	open(ipsi, file=fpsi)
	write(ipsi,'(a20, a20)') "x", "rho"
	do j=0, npts+1,10
		x=x0+dx*dfloat(j)
		write(ipsi,'(f20.10, es20.10)') x, abs(psi(j))**2
	end do
	close(ipsi)
	end if
	!----------------------------------------------------------------------------
	! the 3 loop does the time propagation
	!
	pref=dt*dx2m1/4.d0
	!----------------------------------------------------------------------------
	! set the off diagonal term
	!
	offd=-II*pref
	call cpu_time(time)
	btime(1) = time
	do it=1, nt
		t=dt*it
		!-------------------------------------------------------------------------
		! calculate phi and initialize the diagonal and off diagonal elements
		!
		do j=1, npts
			phi (j) = (1.-(0.,.5)*dt*vpot(j))*psi(j)+(0.,1.)*pref*(psi(j+1)+psi(j-1))
			diag(j) = (1.+(0.,.5)*dt*vpot(j))
		end do
		!-------------------------------------------------------------------------
		! sweep down the index to do the LU decomposition
		!
		do j=2, npts
			c       = offd/diag(j-1)
			phi (j) = phi (j) - phi(j-1)*c
			diag(j) = diag(j) - offd*c
		end do
		!-------------------------------------------------------------------------
		! sweep up the index to do the back substitution
		!
		do j=npts, 1, -1
			psi(j)=(phi(j)-offd * psi(j+1)) / diag(j)
		end do
		!psit(:,it) = psi
		!-------------------------------------------------------------------------
	end do
	call cpu_time(time)
	etime(1) = time
	dtime(1) = etime(1) - btime(1)
	write(*,'("Time",x,f8.2,x,"s")') dtime(1)
	!----------------------------------------------------------------------------
	if (ifout) then
	irecl=(nt+1) * 20
	open(ipsit, file=fpsit, recl=irecl)
		!-------------------------------------------------------------------------
		write(ipsit, '(a20, $)') "x"
		do it=1, nt
			t=dt*it
			if(mod(it,nmod) == 0) then
				write(ctmp,'(f10.4)') t
				write(ipsit, '(a20, $)') "f("//trim(adjustl(ctmp))//")"
			end if
		end do
		write(ipsit,*)
		!-------------------------------------------------------------------------
		do j=0, npts+1, 10
			x=x0+dfloat(j)*dx
			write(ipsit, '(f20.10, $)') x
			do it=1, nt
				if(mod(it,nmod) == 0) then
					write(ipsit, '(es20.10, $)') abs(psit(j,it))**2
				end if
			end do
			write(ipsit,*)
		end do
		!-------------------------------------------------------------------------
	close(ipsit)
	end if
	call cpu_time(time)
	etime(1) = time
	dtime(1) = etime(1) - btime(1)
	write(*,'("Time",x,f8.2,x,"s")') dtime(1)
	!----------------------------------------------------------------------------
 900   format(2f20.10)
	!----------------------------------------------------------------------------
	deallocate(vpot)
	deallocate(xt  )
	deallocate(psi )
	deallocate(psit)
	deallocate(phi )
	deallocate(phit)
	deallocate(diag)
	!----------------------------------------------------------------------------
	close(ilog)
	!----------------------------------------------------------------------------
	stop
end

