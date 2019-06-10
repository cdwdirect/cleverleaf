c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
c Description:   F77 routine to compute stable dt in 3d.
c
define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine stabledt3d(dx,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ngc0,ngc1,ngc2,
     &  advecspeed,uval,stabdt)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(FORTDIR/const.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL stabdt,dx(0:NDIM-1)
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     & ngc0,ngc1,ngc2
c
      REAL
     &  advecspeed(0:NDIM-1),
     &  uval(CELL3dVECG(ifirst,ilast,ngc))
c
      REAL maxspeed(0:NDIM-1)
c
      maxspeed(0)=zero
      maxspeed(1)=zero
      maxspeed(2)=zero

      maxspeed(0) = max(maxspeed(0), abs(advecspeed(0)))
      maxspeed(1) = max(maxspeed(1), abs(advecspeed(1)))
      maxspeed(2) = max(maxspeed(2), abs(advecspeed(2)))

c     Do the following with checks for zero
c      stabdt = min((dx(1)/maxspeed(1)),(dx(0)/maxspeed(0)))
c      stabdt = min((dx(2)/maxspeed(2)),stabdt)

      if ( maxspeed(0) .EQ. 0.0 ) then
         if( maxspeed(1) .EQ. 0.0 ) then
            stabdt = 1.0E9
         else 
            stabdt = dx(1)/maxspeed(1)
         endif
      elseif ( maxspeed(1) .EQ. 0.0 ) then
            stabdt = dx(0)/maxspeed(0) 
      else
         stabdt = min((dx(1)/maxspeed(1)),(dx(0)/maxspeed(0)))
      endif

      if (maxspeed(2) .NE. 0.0 ) then
         stabdt = min((dx(2)/maxspeed(2)),stabdt)
      endif

      return
      end
