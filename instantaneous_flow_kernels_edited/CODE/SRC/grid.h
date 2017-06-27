! Parameters for plotting data out
      real     , parameter :: lat1 = -90
      real     , parameter :: lat2 =  90
      real     , parameter :: lon1 =-179
      real     , parameter :: lon2 = 180
      integer  , parameter :: inc = 1
      integer  , parameter :: grid=
     &   (lat2 - lat1)/(inc)+1
      integer  , parameter :: grid2=
     &   (lon2 - lon1)/(inc)+1
      real     , parameter :: Vspac=inc
      real     , parameter :: Hspac=inc
