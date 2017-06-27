      subroutine fldcalc(GeoK,VelK,SurK,CMBK)
      implicit none

      include "size.h"
      include "grid.h"
      include "parameter.h"
      
      common /ang/ theta, phi
      common /rad/ r, rho0
      common /fld/ anom, degMOD
      common /path/ output, input
      common /ResCorr/ Corr, surmax, geomax

      character*120 output, input
      real r(nr+1)      , rho0(nr+1)

      real theta(grid)  , theta2(grid)
      real phi(grid2)   , phi2(grid2)

      real GeoK(l_min:l_max, nr+1), VelK(l_min: l_max, nr+1)
      real SurK(l_min:l_max, nr+1), CMBK(l_min: l_max, nr+1)

      real sqr, powone, facto(0:2*l_max),facto2(0:2*l_max)
      real geomax, surmax

      real geo(grid,grid2), leg(0:l_max,0:l_max,grid)
      real anom(0:l_max, -l_max:l_max, 1:nr+1)
      real sphharm

      real geoidCO (0:l_max, -l_max:+l_max) 
      real SurTopCO(0:l_max, -l_max:+l_max) 
      real CMBTopCO(0:l_max, -l_max:+l_max) 
      real veloCO  (0:l_max, -l_max:+l_max)
      real RefGeo  (0:l_max, -l_max:+l_max)

      real geoidCO_3d  (0:l_max, -l_max:+l_max,1:nr+1)
      real SurTopCO_3d (0:l_max, -l_max:+l_max,1:nr+1)
      real CMBTopCO_3d (0:l_max, -l_max:+l_max,1:nr+1)

      real geoid (grid , grid2) 
      real SurTop(grid , grid2) 
      real CMBTop(grid , grid2) 
      real velo  (grid , grid2)
      real GeoFld(grid , grid2)

      real Corr, SumTemp1, SumTemp2, mfcalc


      real laysize

      integer i, l, m, j, degMOD

      laysize = (r(nr+1)- r(1))/(nr)

      ! convert to 'degrees', Round up to the 2 digits precission

      do i=1,grid
          theta(i)=(90.0-(lat1+(grid-i)*Vspac))*deg2rad
      enddo
      do i=1,grid2
          phi(i)=(lon1+(i-1)*Hspac)*deg2rad
      enddo

!     Legendre Polynomials
      do i=1,grid
          leg(0,0,i)=1
          leg(1,0,i)=cos(theta(i))
          leg(1,1,i)=-sqrt(1-leg(1,0,i)**2)
!          leg(1,1,i)=-sqrt(1-cos(theta(i))**2)
      enddo
      if(l_max>1) then
          do l=2,l_max
              do m=0,l
                  do i=1,grid
                      if (m==0) then
                          leg(l,m,i)=
     & ((2*l-1)*cos(theta(i))*leg(l-1,m,i)-(l-1)*leg(l-2,m,i))/l
                      else
                          leg(l,m,i)=
     & (leg(l-2,m,i)-(2*l-1)*sqrt(1-cos(theta(i))**2)*leg(l-1,m-1,i))
                      endif
                  enddo
              enddo
          enddo
      endif

!      open(22222, file='ThetaPhi2')
!      do i = 1,grid
!         do j = 1,grid2
!            write(22222, *) phi(j)*rad2deg, (pi/2-theta(i))*rad2deg
!         enddo
!      enddo
!      close(22222)
!      stop

! all square root factorials and square roots
      facto(:)=1.0
      facto2(:)=1.0
      do l=1,2*l_max
         if(l>170) then
             facto2(l)=l*facto2(l-1)
             facto(l)=facto(l-1)
         else
             facto(l)=l*facto(l-1)
         endif
      enddo
      facto(:)=sqrt(facto(:))
      facto2(:)=sqrt(facto2(:))


      geoidCO (:,:) = 0.0
      veloCO  (:,:) = 0.0
      SurTopCO(:,:) = 0.0
      CMBTopCO(:,:) = 0.0
      RefGeo  (:,:) = 0.0

      geoid (:,:) = 0.0
      velo  (:,:) = 0.0
      SurTop(:,:) = 0.0
      CMBTop(:,:) = 0.0
      GeoFld(:,:) = 0.0

      Corr = 0.0
      SumTemp1=0.0
      SumTemp2=0.0


      if (CorrFlag==1) then
         open(97,file='./INPUT/RefGeo/eigen5c2',
     &                  action='read',status='old')  
         do l=0,l_max
             do m=0,l
                 if (m==0) then
                     read(97,*) RefGeo(l,m)
                 else
                     read(97,*) RefGeo(l,-m),RefGeo(l,+m)
                 endif
      ! Sia: This section is copied from Andre's code, 
      !        These Corrections are just applied for the reference GEOID
                 if(m==0.and.l==2) RefGeo(l,m)=RefGeo(l,m)
     &  + 1072.618e-6/sqrt(2.0*l+1.0)         ! Nakiboglu 1982
                 if(m==0.and.l==4) RefGeo(l,m)=RefGeo(l,m)
     &  - 2.992e-6/sqrt(2.0*l+1.0)            ! Nakiboglu 1982
             enddo
         enddo
      ! Sia: geodetic normalization contained in EIGEN data, What Andre' does
         !RefGeo=RefGeo*Emass*G/Rearth/sqrt(2.0*pi)    
         RefGeo=Rearth*RefGeo    
         close(97)
      endif       


      do l=2, l_max_calc
         do m= -l, l
            if(l>=l_min) then
               do j=nr+1,1,-1
               geoidCO(l,m) =
     &            geoidCO(l,m) + anom(l,m,j)*GeoK(l,j)*laysize
               veloCO(l,m)  =
     &            veloCO(l,m)  + anom(l,m,j)*VelK(l,j)*laysize
               SurTopCO(l,m)=
     &            SurTopCO(l,m)+ anom(l,m,j)*SurK(l,j)*laysize
               CMBTopCO(l,m)=
     &            CMBTopCO(l,m)+ anom(l,m,j)*CMBK(l,j)*laysize
            !  Computation of the 3d fields
               if (visual3d ==1) then
                  geoidCO_3d(l,m,j) =
     &                  anom(l,m,j)*GeoK(l,j)*laysize
                  SurTopCO_3d(l,m,j)=
     &               anom(l,m,j)*SurK(l,j)*laysize
                  CMBTopCO_3d(l,m,j)=
     &               anom(l,m,j)*CMBK(l,j)*laysize
               endif
               enddo
!               Corr=Corr+geoidCO(l,m)*RefGeo(l,m)
!               SumTemp1= SumTemp1 + geoidCO(l,m)**2
!               SumTemp2= SumTemp2 + RefGeo(l,m)**2
               do i=1,grid
                  do j=1,grid2
                     if(m==0) then
                        sphharm=sqr(2*l+1)*0.5*sqpi*leg(l,m,i)
                     elseif(m<0.and.m>-l-1) then
                        sphharm=sqr(2*l+1)*sqpi/sqr(2)*facto(l+m)
     & /facto(l-m)*powone(m)*leg(l,-m,i)*
     & facto2(l+m)/facto2(l-m)*cos(m*phi(j))
                     elseif(m>0.and.m<l+1) then
                        sphharm=sqr(2*l+1)*sqpi/sqr(2)*facto(l-m)/
     & facto(l+m)*powone(m)*leg(l,m,i)*
     & facto2(l-m)/facto2(l+m)*sin(m*phi(j))
                     endif
                     geoid(i,j) =geoid(i,j) + geoidCO(l,m)*sphharm
                     velo (i,j) =velo (i,j) +  veloCO(l,m)*sphharm
                     SurTop(i,j)=SurTOP(i,j)+SurTOPCO(l,m)*sphharm
                     CMBTop(i,j)=CMBTop(i,j)+CMBTopCO(l,m)*sphharm
                     GeoFld(i,j)=GeoFld(i,j)+RefGeo  (l,m)*sphharm
                     if (geomax < abs(geoid(i,j) )) then
                         geomax=abs(geoid (i,j))
                     endif
                     if (surmax < abs(SurTop(i,j))) then 
                        surmax=abs(SurTop(i,j))
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo

      if (CorrFlag ==1) then
         Corr = Corr/sqrt(SumTemp1*SumTemp2)
      endif
     
!     Output of Spherical Harmonics
      call sphcout(geoidCO,'geoidCO')
      call sphcout(veloCO,'velCO')
      call sphcout(SurTopCO,'SurTopCo')
      call sphcout(CMBTopCO,'CMBTopCO')
      call sphcout(RefGeo,'RefGeoCO')

!     Output of Grid
      call fldout(SurTOP,       'SurTOP')
      call fldout(geoid ,       'Ge-oid')
      call fldout(velo  ,       'VelTOP')
      call fldout(CMBTop,       'CMBTOP')
      call fldout(GeoFld,       'RefGeo')     
      call fldout(SurTOP-geoid, 'DynTop')     


!     Write out the 3D-fields 
!      if (visual3d==1) then
!         call sphcout_3d(geoidCO_3d  ,'geoidCO_3d')
!         call sphcout_3d( SurTopCO_3d,'SurTopCo_3d')
!         call sphcout_3d( CMBTopCO_3d,'CMBTopCO_3d')
!      endif
      end subroutine fldcalc


! * function sqri
!     Needed for some of the calculations
      real function sqr(i)
      implicit none

      integer i

      sqr = sqrt(dble(i))

      return
      end function sqr

! *  function powone
!  To reduce the load of (-1)^m
      real function powone(m)
      implicit none

      integer m

      if (mod(m,2) .eq. 0) then
         powone = 1.0
      else
         powone = -1.0
      endif
      return
      end function powone
