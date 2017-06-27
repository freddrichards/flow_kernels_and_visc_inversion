*dk progress
      subroutine progress(j)
      implicit none
      integer(kind=4)::j,k
      character(len=37)::bar="???% |                              |"
      write(unit=bar(1:3),fmt="(i3)") 10*j/3
      do k=1, j
        bar(6+k:6+k)="*"
      enddo
      ! print the progress bar.
      write(unit=6,fmt="(a1,a1,a37)") '+',char(13), bar
      return
      end subroutine progress

*dk SpecReadIn
      subroutine SpecRef(Iunit,refgeo)
      implicit none
      
      include 'size.h'
      include 'parameter.h'
      integer l,m,i, Iunit
      real refgeo(0:l_max,-l_max:+l_max)
      do l=0,l_max
          do m=0,l
              if (m==0) then
                  read(Iunit,*) refgeo(l,m)
              else
                  read(Iunit,*) refgeo(l,+m),refgeo(l,-m)
              endif
              !??? Why: m=0,  l=2 and l=4 should be changed
              if(m==0.and.l==2) refgeo(l,m)=refgeo(l,m) +
     &  1072.618e-6/sqrt(2.0*l+1.0)     ! Nakiboglu 1982
              if(m==0.and.l==4) refgeo(l,m)=refgeo(l,m) -
     &  2.992e-6/sqrt(2.0*l+1.0)            ! Nakiboglu 1982
              !refgeo(l,m)=refgeo(l,m)*(l-1)
             !if(m/=0) refgeo(l,-m)=refgeo(l,-m)*(l-1)
          enddo
      enddo
      !refgeo=refgeo*Emass*G/Rearth/sqrt(2.0*pi)        ! geodetic normalization contained in EIGEN data
      refgeo=refgeo*Emass*G/Rearth        ! geodetic normalization contained in EIGEN data
      !refgeo=refgeo/R**3.0*G*Emass
!      read(Iunit,*)
!      do l=0,l_max
!         do m=0,l
!            read(Iunit,*) i,j,temp1, temp2
!            if (i .lt. l_max) then
!               if (m==0) then
!                  anom(l, m) = temp1
!               else
!                  anom(l, m) = temp2
!                  anom(l,-m) = temp1
!               endif
!            endif
!         enddo
!      enddo
      
      end subroutine SpecRef 


*dk initialize
c..   Using read_in_radial reads in the radial terms
      subroutine initialize
      implicit none

      include "size.h"
      include "parameter.h"

      common /rad/ r, rho0
      common /fld/ anom, degMOD
      common /path/ output, input

      real r(nr+1), rho0(nr+1)
      real anom(0:l_max,-l_max:+l_max,1:nr+1)
      integer dump,l, i, m, degMOD, counter
      character*120 output, input
      character*40 temp

      open(001, file= trim(input), status='old', action='read')
      read(001, '(A)') temp
      read(temp(8:11), *) degMOD
      read(temp(16:19), *) dump   

      if (l_max .ne. degMOD) then
         write(*,'(A)') 'l_max in size.h .ne. l_max in the input!'
         write(*,'(A)') 'Please set the l_max in size.h equal
     &  to l_max in the input, right now it cannot be done!'
         stop
      endif

      if ((dump .ne. nr) ) then
         write(*,'(A)') 
     & '# of layers in input file and size.h inconsistent'
      endif
     
      read(001,*) temp
      read(001,*) temp
      call read_in_radial(r,001)
      !Sia: To imitate what Andre does, Does not change a thing!
      ! There is no need to do that anymore, because it does not make
      !     any sense to me now!
!       do i=1,nr+1,1
!          r(i) = anint(r(i)/100)*100
!       enddo

      read(001,*) temp
      call read_in_radial(rho0,001)
      read(001,*) temp
  
!     TO DO :  To apply the limits of density applications
      do i=nr+1, 1, -1
         read(001,*) temp
         do l=0, l_max
            read(001,*) anom(l,0,i)
            anom(l,0,i)=anom(l,0,i)*sqt
!            anom(l,0,i)=anom(l,0,i)*rho0(i)*sqt
            if (l>0) then
               do m=1,l
                  read(001,*)  anom(l,-m,i),anom(l,+m,i)
                  anom(l, m,i)=anom(l, m,i)*sqt
!                  anom(l, m,i)=anom(l, m,i)*rho0(i)*sqt
                  anom(l,-m,i)=anom(l,-m,i)*sqt
!                  anom(l,-m,i)=anom(l,-m,i)*rho0(i)*sqt
                enddo
            endif
         enddo
         if ((r(nr+1)-r(i))<(top_depth*1000) .or.
     &  (r(i)-r(1))<(bot_depth*1000)   ) then
             anom(:,:,i)=0.0
         endif
      enddo
!     TO DO : HERE MAY COME OTHER READINGS
      close(001)
      end subroutine

*dk read_in_radial
c..   Read in radial profile u
      subroutine read_in_radial(u, in_unit)
      implicit none                       

      include "size.h"
      integer in_unit, i
   
      real u(nr+1)
c..   Be carfule of the indexing, increament
      do i=nr+1,1,-1
         read(in_unit, *) u(i)
      enddo
      end subroutine read_in_radial


*dk vector_out
c..   writes a ni vector on the iunit
      subroutine vector_out(A,ni, iunit)
      implicit none

      integer j
      integer ni, iunit
      real A(ni)
      write(iunit,'(A)', advance='yes') ' '
      do j=1,ni
         write(iunit, 10, advance='no') A(j)
      enddo
      write(iunit,'(A)', advance='yes') ' '

 10   format(f9.5,2x)
      end subroutine vector_out

*dk matrix_out
c..   writes a ni*nk matrix on the iunit
      subroutine matrix_out(A,ni,nk, iunit)
      implicit none

      integer j,k
      integer nk, ni, iunit
      real A(nk, ni)
      do j=1,ni
         do k=1,nk
            write(iunit, 10, advance='no') A(j,k)
         enddo
         write(iunit,'(A)', advance='yes') ' '
      enddo
 10   format(e10.3,2x)
      end subroutine matrix_out

*dk kernel_out
      subroutine kernel_out( u, io , l)
      implicit none

      include "size.h"

      common /rad/ r, rho0
      common /path/ output, input

      real r(nr+1), rho0(nr+1)
      real u(nr+1)
      integer l, i, io
      character*20 char1, char2
      character*120 output, input

      if (io .eq. 1) then
         char1 = 'Kernel'
      elseif (io .eq. 2) then
         char1 = 'VelKernel'
      elseif (io .eq. 3) then
         char1 = 'TopCMBKernel'
      elseif (io .eq. 4) then
         char1 = 'TopSurfKernel'
!      elseif (io .eq. )
      elseif (io .eq. 5) then
         char1 = 'GravAnoKernel'
      elseif (io .eq. 6) then
         char1 = 'MolnarKernel'
      else
         write(*,*) 'No valid io for kernel_out!'
      endif
      call system('mkdir -p '//trim(output)//trim(char1)//'/')

      write(char2, '(I3.3)') l

      open(333, file=trim(output)//trim(char1)//'/'//trim(char1)
     & //trim(char2), action = 'write', status='unknown')

      do i=1,nr+1,1
         write(333, '(f17.2,4X, es16.7e3)') r(i)/1000, u(i)
      enddo

      close(333)

      end subroutine kernel_out


      subroutine sphcout_3d(coef, char1)

      common /path/ output, input
      common /rad/ r, rho0

      include "parameter.h"
      include "size.h"
      include "grid.h"

      integer l, m
      real coef (0:l_max, -l_max:+l_max,1:nr+1)
      real r(nr+1)
      character*15 char1
      character*120 output, input

      open(002,file=trim(output)//'/'//trim(char1),action='write')

      write(002,'(A,I4,A,I4,A)') '# LMAX=',l_max, ', NR',nr,' 1'
      write(002,'(A,A)') '# Temperature: INPUT : ', trim(input)
      write(002,'(A,A)') '# Radius'

      do j= nr+1,1,-1
         write(002,'(f15.2)') r(j)
      enddo

      write(002,'(A,A)') '# Averages'
      do j= nr+1,1,-1
         write(002,'(f10.2)') coef(0,0,j)
      enddo

      write(002,'(A,A)') '# Spherical Harmonics'
      do j= nr+1,1,-1
         write(002,'(A,I4.1)') '# Radial Index:', j
          do l=0,l_max_calc
             do m=0,l
                if (m==0) then
                   write(002,'(e15.6, e15.6)')
     &   coef(l,m,j), 0.0
               else
                  write(002,'(e15.6, e15.6)')
     &   coef(l,-m,j), coef(l,+m,j)
                endif
             enddo
          enddo
      enddo

      close(002)

      end subroutine sphcout_3d


      subroutine sphcout(coef,char1)
      implicit none

      common /path/ output, input

      include "parameter.h"
      include "size.h"
      include "grid.h"

      integer l, m
      real coef (0:l_max, -l_max:+l_max)
      character*15 char1
      character*120 output, input

      open(444, file=trim(output)//trim(char1),
     & action='write', status='unknown')
!     The header part
      write(444,'(A)') '... HEADER'
      write(444,'(A,A)') '... INPUT : ', trim(input)
      write(444,'(A,A)') '... ', trim(output)
      write(444,'(A)') 'l ,m , Cnm(m>=0) , Snm(m<0)'
      if(l_max>1) then
         do l=2,l_max_calc
            do m=0,l
               if (m==0) then
                   write(444,'(A,I8,I8,e15.6, e15.6)')
     &   'SPHCO',l, m, coef(l,m), real(0)
               else
                  write(444,'(A,I8,I8,e15.6, e15.6)')
     &   'SPHCO',l, m, coef(l,-m), coef(l,+m)
               endif
            enddo
         enddo
      endif
      write(444,'(A)') '... END OF FILE!'
      close(444)

      end subroutine sphcout

      subroutine fldout(field, char1)

      implicit none

      include "parameter.h"
      include "size.h"
      include "grid.h"

      common /ang/ theta, phi
      common /path/ output, input

      real field(grid, grid2)
      real theta(grid), phi(grid2)
      integer i, j, fldoutid
      character*6 char1
      character*120 output, input


      open(645, file=trim(output)//char1, status = "unknown",
     &  action = 'write')
      

      do i=1,grid
         do j=1,grid2
            write(645,'(f7.2,3X,f7.2,3X,f15.7)')
     &  phi(j)*rad2deg,  (pi/2-theta(i))*rad2deg, field(i,j)
         enddo
      enddo
      

      close(645)

      end subroutine fldout



*dk summary 
      subroutine summary
      implicit none

      include "size.h"
      include "parameter.h"
      include "grid.h"

      common /rad/ r, rho0
      common /vis/ rvsc
      common /path/ output, input
      common /ResCorr/ Corr, surmax, geomax
      common /fld/ anom, degMOD

      integer i,l,m,degMOD
      real r(nr+1), rho0(nr+1)
      real anom(0:l_max,-l_max:+l_max,1:nr+1)
      real rvsc(nr+1), kappa
      real Corr, grav, surmax, geomax
      character*120 output, input

      open(111, file=trim(output)//'Summary', action = 'write',
     &     status = 'unknown' )

      write(111, '(5X,A,A)') 'Summary: Input File - ', trim(input)
      write(111, '(15X, A)') 'Saved to: ',trim(output)

      do i = 1,62
         write(111,'(A)', advance = 'no') '*'
      enddo
      write(111,*) 

      write(111, 22)nr*2, nr+1, l_min, l_max_calc, top_depth, bot_depth,
     & bdry, rho_c, lat1, lat2, lon1, lon2, Corr, geomax, surmax
      

 22   format
     &(' MT    =',I10  ,'| NR+1    =',I10  ,'| LMIN    = ',I10  ,/,
     & ' LMAX  =',I10  ,'| TOP_DPTH=',f10.2,'| BOT_DPTH= ',f10.2,/,
     & ' bdry  =',I10  ,'| RHO_C   =',e10.2,'| LAT1    = ',f10.2,/,
     & ' LAT2  =',f10.2,'| LON1    =',f10.2,'| LON2    = ',f10.2,/,
     & ' Corr  =',f10.2,'| MaxGeoid=',f10.2,'| MaxSurTo= ',f10.2 )

      do i = 1,62
         write(111,'(A)', advance = 'no') '*'
      enddo
      write(111,*) 

      write(111, '(3x,3A,A,A)')'Radius ',' Viscosity '
     &,' Density ',' CHI ',' Gravity ' 
      do i = nr+1, 1, -1        
         write(111,'(f15.0, e10.3, f8.2, e10.1, f7.2)')
     &            r(i),rvsc(i)*visc0, rho0(i), kappa(i), grav(r(i))
      enddo
      close(111)

      open(111,file=trim(output)//'Density', action = 'write',
     &     status = 'unknown')
      
      write(111,'(A)') ' l  m     Radius   Cnm   Snm'
      do i=nr+1,1,-1
         do l=0,l_max
            do m=0,l
               if (m==0) then
                  write(111,'(I2.2,3X,I2.2,3X,f11.3,3X,E12.4,3X,E12.4)')
     &      l,m,r(i),anom(l,m,i), 0.0        
               else                                    
                  write(111,'(I2.2,3X,I2.2,3X,f11.3,3X,E12.4,3X,E12.4)')
     &      l,m,r(i),anom(l,m,i),anom(l,-m,i)
               endif
            enddo
         enddo
      enddo

      close(111)
      end subroutine summary

*dk viscinit
      subroutine viscinit
      implicit none
      include 'size.h'

      common /path/ output, input, name_prof
      common /vis/  rvsc


      integer i
      character*120 input, output, name_prof
      real rvsc(nr+1)

      open(183, file=trim(name_prof),
     &  action='read', status='old')
      
      do i=nr+1,1,-1
         read(183,*) rvsc(i)
      enddo

      close(183)

      end subroutine viscinit

