      integer, parameter:: visual3d = 0

      integer,parameter:: lwmax=1000 
		  
!     misfit weightings (w1=dyntopo,w2=geiod;w3=CMB)
      real  , parameter:: w1=1.0
      real  , parameter:: w2=1.0 
      real  , parameter:: w3=0.0


!     Which type kappaindex==1 compressible else: incompressible
      integer, parameter:: kappaindex=1
!     Starting Depth for calculating the gravity based on
      real  , parameter:: top_depth = 200 
      real  , parameter:: bot_depth =-10

!     Boundary condition on top
      integer, parameter:: bdry=1

!     Boundary condition on bot
      integer, parameter:: bdry_bot=1

!     If Correlation with RefGeo should be calculated
      integer, parameter:: CorrFlag = 0

!     Background density
      real, parameter:: bgrho = 4.500e+03

!     Background Viscosity
      real, parameter:: visc0 = 1      

!     The changes of the Dynamic Topography at CMB surface
!           has to do with rho_c
      real, parameter:: rho_c =  9.90e3

!     The Dynamic topography on the surface 
      real, parameter:: rho_w = 0.0e3

!     Mathematical constants
      real, parameter:: G  = 6.67e-11
      real, parameter:: pi = 3.141592
      real, parameter:: sqpi = 1.0/sqrt(pi)
      
!     Geodetical Rearth
      real, parameter :: Rearth=6370000.0
      real, parameter :: Emass=5.974e24
   
!     Factor for conversion of anamolies ( Spherical Coefficients)
      real, parameter:: sqt=sqrt(4.0*pi)

!     Conversion of rad2grad or Vice Versa
      real  , parameter:: deg2rad = pi/180
      real  , parameter:: rad2deg = 180/pi


