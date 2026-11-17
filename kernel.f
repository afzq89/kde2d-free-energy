program kernel_estimator

!2D Kernel estimator by Mahmoud Moradi
!Update and example provided by Ashkan Fakhar
!Convert scattered 2D samples into a smooth, regular grid using Gaussian kernels.
!K(X,Y) = exp[-(Y-y)^2/2h_y^2-(X-x)^2/2h_x^2]

!Syntax `./kernel < monitor > BFE.txt` 
!monitor file includes the points; the first line must give the number of points
!input is like x y w
!"kernel.in" is the input file of the fortran kernel program
!and needs to be present. The first 9 lines are read as follows
!(1)x_min the range of x
!(2)y_min the range of y
!(3)x_max the range of x
!(4)y_max the range of x
!(5)n_x number of output data
!(6)n_y number of output data
!(7)h_x kernel width; Larger values produce smoother maps
!(8)h_y kernel width; Smaller values preserve more local variation.
!(9)Delta_h to avoid unnecessary calculations (data farther than Delta_h*h_x,y is not considered)
! dh determines how many kernel widths around each point are evaluated.
!(10) temp
!The output files are:
!kernel.out
!progress.out
!histogram.txt
!the derivatives may not be accurate!
!BFE.txt that includes 5 columns {col1:x col2:y col3:f(x,y) col4:df/dx col5:dr/dy}


    implicit none
    integer, parameter :: dp = kind(1.d0)
    real(dp), parameter :: kB = 8.31441E-3/4.184;
    real(dp), dimension(:,:), allocatable:: z, zx, zy
    integer :: n, k, i, j, ii, jj, hi, hj, nx, ny, n100, i0, j0, im, jm, dh
    real(dp) :: x0, y0, xm, ym, hx, hy, dx, dy, BFE, dBFEx, dBFEy, Nz, dz, dzx, dzy, sc
    real(dp) :: temp, x, y, w, sx, sy, xii, yjj, xi, yj, ex, dzm
    real(dp) :: kBT;
    
    open (1,file='kernel.in')
    read(1,*) x0
    read(1,*) y0
    read(1,*) xm
    read(1,*) ym
    read(1,*) nx
    read(1,*) ny
    read(1,*) hx
    read(1,*) hy
    read(1,*) dh
    read(1,*) temp
    close(1)

    kBT = kB * temp;
    
    open (2,file='kernel.out')
    write(2,*) 'x0 = ',x0, ' y0 = ',y0
    write(2,*) 'xm = ',xm, ' ym = ',ym
    write(2,*) 'nx = ',nx, ' ny = ',ny
    write(2,*) 'hx = ',hx, ' hy = ',hy,' dh = ',dh
    
    dx = (xm - x0) / nx
    dy = (ym - y0) / ny
    
    hi = int((hx/dx) * dh)
    hj = int((hy/dy) * dh)
    
    sx = 2 * hx * hx
    sy = 2 * hy * hy

    allocate(z(0:nx, 0:ny))
    allocate(zx(0:nx, 0:ny))
    allocate(zy(0:nx, 0:ny))
    z = 0.0
    zx = 0.0
    zy = 0.0
    sc = 0.0
    read(*,*) n
    write(2,*)'dx = ',dx,' dy =',dy
    write(2,*)'number of samples = ',n
    close(2)
    open(3, file='progress.out')
    write(3,*)'sampling started...'
    n100 = n/100
    do k = 1, n
        read(*,*) x, y, w
        if (k == n100*(k/n100)) write(3,*)k,' samples out of ', n, ' samples'
	if ((x > x0).and.(x < xm).and.(y > y0).and.(y < ym)) then
	    i = int((x - x0) / dx)
	    j = int((y - y0) / dy)
	    if (i-hi<0) then
		i0 = 0	
	    else 
		i0=i-hi
	    end if
	    if (j-hj<0) then
		j0 = 0	
	    else 
		j0=j-hj
	    end if
	    if (i+hi>nx-1) then
		im = nx-1
	    else 
		im=i+hi 
	    end if
	    if (j+hj>ny-1) then
		jm = ny-1
	    else 
		jm=j+hj 
	    end if
	    dzm = 1.0
	    do ii = i0, im
    		xii = x0 + ii * dx
	    	ex = w*exp(-(xii - x)**2 / sx)
	    	do jj = j0, jm
	    	    yjj = y0 + jj * dy
	    	    dz = exp(-(yjj - y)**2 / sy) * ex
	    	    dzx = xii * dz
	    	    dzy = yjj * dz
	    	    if (dz<dzm) dzm = dz
	    	    z(ii, jj) = z(ii, jj) + dz
	    	    zx(ii, jj) = zx(ii, jj) + dzx
	    	    zy(ii, jj) = zy(ii, jj) + dzy
	    	    sc = sc + dz
	    	end do
	    end do
        end if
    end do
    close(3)

    open(4, file='histogram.txt')

    do i = 0, nx-1
	do j = 0, ny-1
	    xi = x0 + i * dx
	    yj = y0 + j * dy
	    if (z(i,j) == 0) then 
		Nz = 0
		BFE = -kBT * log(dzm/sc)
		dBFEx = 0
		dBFEy = 0
	    else 
		Nz = z(i,j)
		BFE = -kBT * log(Nz/sc)
		dBFEx = (-2 * kBT / sx) * (xi - zx(i,j) / Nz)
		dBFEy = (-2 * kBT / sy) * (yj - zy(i,j) / Nz)
	    end if
	    write(*,'(E13.6E2," ",E13.6E2," ",E13.6E2," ",E13.6E2," ",E13.6E2)\n') xi, yj, BFE, dBFEx, dBFEy
	    write(4,*) xi, yj, Nz
	end do
    end do
    close(4)
end program kernel_estimator
