      program driver
      implicit none

      common /path/ output, input, name_prof
      
      integer N, sufi, suff, step
      character*120 output, input, name_prof
      real :: start, finish

 
      N = IARGC()
      
      if (N .ne. 3) then
         write(*,'(A)') 'Number of given arguments should be 3!'
         write(*,'(I1.1, A)') N, ' is given. Try the following order:'
         write(*,'(A)')
     &        './KERNELS_OUT [INPUT] [OUTPUT] [VISC_PROFILE]'
         stop
      endif

      call getarg(1,input)       ! Obtain the path+prefix
      call getarg(2,output)      ! char1 is sufi
      call getarg(3,name_prof)   ! char1 is sufi
      
      call prop_matrix

      end program driver

!  -------------------------------------------------
contains
! --------------------------------------------------

      subroutine prop_matrix

      implicit none

      include "size.h" 
      include "parameter.h"
      include "grid.h"
      
      common /rad/ r, rho0
      common /vis/ rvsc
      common /path/ output, input, name_prof
      common /ResCorr/ Corr, surmax, geomax

      integer i, j, k, m, n, ierr
      real kappa, mf
      real rho0(nr+1),r(nr+1),rvsc(nr+1)
      real a(6,6), prop(6,6,1:nr+2)
      real eigen_i(6), eigen_r(6)
      real c(8,8), b(8), u(8)
      real dump, fac, grav, Corr, surmax, geomax
      real oDT(grid, grid2)
      real oG(grid, grid2)
      real oCMB(grid, grid2)
      real :: start, finish
      character*80 output, input, name_prof

      real kernel(l_min:l_max,nr+1), velkernel(l_min:l_max,nr+1)
      real topsurfkernel(l_min:l_max,nr+1)
      real gravanomkernel(l_min:l_max,nr+1)
      real molnarkernel(l_min:l_max,nr+1)
      real topcmbkernel(l_min:l_max,nr+1)
	  
      call system('mkdir -p '//trim(output))

c..   reading the density and redius
      call initialize
      call viscinit

c..   constructing the matrix

      do j= l_min, l_max,1
	 ! have removed carriagecontrol='fortran') as gfortran can no longer handle 
         ! this even with -fdec compile flag
         ! Subroutine progress generates a loading bar on the screen
         ! Don't need to write out progress in inversion
!          call progress(int(dble(j-l_min+1)/dble(l_max-l_min+1)*10)*3) 
!         write(*,'(A, I2.2, A, I2.2)')
!     &       'Kernel Computation of order: ', j, '/', l_max
         ! Initialisation of a 6*6 matrix with zero values

         call matunitinit(prop(:,:,nr+2), 6,6)
         do i=nr+1,2,-1
            call aconstruct(a, i, real(j))
            call eigen(a, 6, eigen_r, eigen_i)
            call propmat(a,prop(:,:,i),eigen_r,eigen_i,r(i)/r(i-1))
            prop(:,:,i) = matmul(prop(:,:,i+1), prop(:,:,i))
         enddo
c..      prop(:,:,2) represents prop{6*6}(c,a)
         call cconstruct(C, j, prop(:,:,2))
         do i=nr+1,1,-1
            call bconstruct(i, prop(:,:,i+1), j, b)
            call solve_lgs(C, b, 8, u)
!           The  Geoid Kernel
            kernel(j,i)       =
     &        (visc0/(bgrho*r(nr+1)*grav(r(nr+1))))*u(3)
!           Gravity Anomaly Kernel
            gravanomkernel(j,i)       =
     &        1e5*((G*Emass)/r(nr+1)**3)*(j-1)
     &        *(visc0/(bgrho*r(nr+1)*grav(r(nr+1))))*u(3)
!           The top topography kernel
            topsurfkernel(j,i)=
!     &       -(visc0/((rho0(nr+1)-rho_w)*r(nr+1)*grav(r(nr+1))))
     &       -(visc0/((bgrho-rho_w)*r(nr+1)*grav(r(nr+1))))
     &         *(u(1)+(rho_w/bgrho)*u(3))
!           Lateral Velocity kernel 
            velkernel(j,i)=u(2)*sqrt(real(j)*real(j+1.0))
!           The CMB topography kernel
            topcmbkernel(j,i)=
     &       - (visc0/((rho0(1)-rho_c)*r(1)*grav(r(1))))
     &         *(u(5)+(rho_c/bgrho)*u(7))
!           Impedence Kernel
            molnarkernel(j,i) = 
     &           gravanomkernel(j,i)/(-abs((topsurfkernel(j,i)*1e-3)))
         enddo
!        Molranrkernel on CMB is a 0/0 which could go to inf.
         molnarkernel(:,nr+1) = 0.0
         molnarkernel(:,   1) = 0.0 
! don't need to write out kernels in inversion

          call kernel_out(kernel(j,:), 1, j)
          call kernel_out(velkernel(j,:), 2, j)
          call kernel_out(topcmbkernel(j,:), 3, j)
          call kernel_out(topsurfkernel(j,:), 4, j)
          call kernel_out(gravanomkernel(j,:), 5, j)
          call kernel_out(molnarkernel(j,:), 6, j)
      enddo

      call fldcalc(kernel,velkernel,topsurfkernel,topcmbkernel)

      end subroutine prop_matrix

c eigen
      subroutine eigen(a_in, nx, wr_out, wi_out)
      implicit none
      
      integer n
      parameter (n=6)
!      integer lda, ldvr, lwork
!      parameter (lda=n, ldvr=n, lwork=(2+nb)*n)
!      parameter (lda=n, ldvr=n, lwork=(2+nb)*n)
      integer  i,j,nx
      integer  INFO
!      integer INFO, LWORK

      integer          lda, ldvl, ldvr
      parameter        (lda=n, ldvl=n, ldvr=n)
!      integer          lwmax
!      parameter        (lwmax = 1000)
      integer          lwork
      parameter        (lwork = 204)
      real a_in(nx,nx), wr_out(nx), wi_out(nx)
      double precision a(lda, n), dummy(1,1), vr(ldvr, n),
     &               vl(ldvl, n), wi(n), wr(n), work(lwork) !, work(lwork)
!      double precision A(lda,n),vl(ldvl, n),
!     &	wr(n), wi(n), work(lwmax), vr(ldvr,n)
!      data             A/
!     & -1.01d+0, 3.98d+0, 3.30d+0, 4.43d+0, 7.31d+0,
!     &  0.86d+0, 0.53d+0, 8.26d+0, 4.96d+0,-6.43d+0, 
!     & -4.60d+0,-7.04d+0,-3.89d+0,-7.66d+0,-6.16d+0,
!     &  3.31d+0, 5.29d+0, 8.20d+0,-7.33d+0, 2.47d+0,
!     & -4.81d+0, 3.55d+0,-1.51d+0, 6.18d+0, 5.58d+0
!     &                  /
      external dgeev


      do i=1,nx
         do j=1,nx
            a(i,j)=a_in(i,j)
         enddo
      enddo

c..   Use Lapack
!      LWORK=-1
!     call dgeev('V','V',N, A, LDA, WR, WI, VL, LDVL,
!     &	VR, LDVR, WORK, lwork, INFO)
!      LWORK=min(lwmax,int(work(1)))
      call dgeev( 'N', 'N', N, A, LDA, WR, WI, dummy, 1,
     &            VR, LDVR, WORK, LWORK, INFO )

      if (info==0) then
         do i=1,nx
            wr_out(i) = wr(i)
            wi_out(i) = wi(i)
         enddo
      else
         write(*,'(a)') 'Error in dgeev!'
         stop
      endif

      end subroutine eigen



c aconstruct
c.. Constructs the A matrix for a given layer
c..   and a given degree l of SPH
      subroutine aconstruct(a, ind, l)
      implicit none

      include 'size.h'
      include 'parameter.h'
      common /rad/ r, rho0
      common /vis/ rvsc


c     type_in = 0; For A matrix in Panasyuk, Hager, Forte 1995
c             = 1; For A matrix as in Schaber Diploma Thesis

      real fac_visc, l, fac_rho0
      real rvsc(nr+1), a(6,6)
      integer ind, i, j
      real r(nr+1), rho0(nr+1)
      real kappa

c..   Assumption: Viscosity of each layer is the average
c..      of the viscosities above and at the bottom


c (WARNING) SHOULD THIS BE AN AVERAGE OR THE VALUE???
      fac_visc = (rvsc(ind))! + rvsc(ind+1) )/2
      fac_rho0 = (rho0(ind))! + rho0(ind+1) )/2

      do i=1,6,1
         do j=1,6,1
            a(i,j)=0.0
         enddo
      enddo

c..  See Corrieu, Thoraval & Ricard 1994 
c..   Line : 1    As a result of Continuity Equation
      a(1,1)=-2-kappa(ind)
      a(1,2)=l*(l+1)
c..   Line : 2    As a result of \theta component of 
c..               the continuity equation
      a(2,1)=-1
      a(2,2)=1
      a(2,4)=1/fac_visc
c..   Line : 3    \r component of the Stokes Equation
c..               combined with constitutive equation for 
c..               components \tau_{\theta \theta} and  \tau_{\phi \phi}
      a(3,1)=4*(3+kappa(ind))*fac_visc
      a(3,2)=-6*l*(l+1)*fac_visc
      a(3,3)=1
      a(3,4)=l*(l+1)
      a(3,6)=- fac_rho0/bgrho
c..   Line : 4    \theta component of the Stokes Equation
c..               combined with constitutive equation for 
c..               components \tau_{\theta \theta} and  \tau_{\phi \phi}
      a(4,1)=-2*(3+kappa(ind))*fac_visc
      a(4,2)=2*(2*l*(l+1)-1)*fac_visc
      a(4,3)=-1
      a(4,4)=-2
      a(4,5)=-fac_rho0/bgrho
c..   Line : 5 & 6   Simple laws of gravity and potential and the pertubation
c..                  In the potential field
      a(5,5)=1
      a(5,6)=1
      a(6,5)=l*(l+1)
      end subroutine aconstruct


! dk derivation
c..
      subroutine derivation(u, r,  u_der)
      implicit none

      include "size.h"      
      
      integer i
      real u(nr+1), r(nr+1)
      real u_der(nr+1)

      do i=2,nr,1
         u_der(i) = (u(i+1)-u(i-1))/(r(i+1)-r(i-1))
      enddo
      
      u_der(1) = 0.0
      u_der(nr+1) = 0.0


      end subroutine derivation


! dk kappa_calculation
c     Computes the kappa factor for the A matrix at each layer

      real function kappa(ind)
      implicit none

      include 'size.h'
      include 'parameter.h'

      common /rad/ r, rho0
      
      integer ind
      real rho0(nr+1), r(nr+1)

      if (kappaindex==1) then
         if (ind == 1 ) then
            kappa = 0.0
         else 
            kappa =  (r(ind)/(rho0(ind))) * ((rho0(ind)-rho0(ind-1))/
     &       (r(ind)-r(ind-1)))
         endif
      else 
         kappa = 0.0
      endif
      return
      end function

