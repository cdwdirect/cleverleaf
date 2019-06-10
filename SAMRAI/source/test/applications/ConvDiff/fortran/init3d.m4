c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
c Description:   F77 routine to initialize 3d convection diffusion equation.
c
define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine initsphere3d(dx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gcw0,gcw1,gcw2,
     &  var,
     &  i_var,o_var,
     &  center, radius,
     &  nequ)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
c***********************************************************************     
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer gcw0,gcw1,gcw2
      integer nequ
      REAL i_var(0:nequ-1),o_var(0:nequ-1)
      REAL radius,center(0:NDIM-1)
      REAL
     &     dx(0:NDIM-1),xlo(0:NDIM-1),xhi(0:NDIM-1)
c
c variables in 3d cell indexed         
      REAL
     &     var(CELL3dVECG(ifirst,ilast,gcw),0:nequ-1)
c
c***********************************************************************     
c
      integer ic0,ic1,ic2,ineq
      REAL xc(0:NDIM-1),x0,x1,x2

      do ic2=ifirst2,ilast2
        xc(2) = xlo(2)+dx(2)*(dble(ic2-ifirst2)+half)
        x2 = xc(2) - center(2)

        do ic1=ifirst1,ilast1
          xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
          x1 = xc(1)-center(1)

          do ic0=ifirst0,ilast0
            xc(0) = xlo(0)+dx(0)*(dble(ic0-ifirst0)+half)
            x0 = xc(0)-center(0)

            do ineq=0, nequ-1
              if ((x0**2+x1**2+x2**2).lt.radius**2) then
                 var(ic0,ic1,ic2,ineq) = i_var(ineq)
              else
                 var(ic0,ic1,ic2,ineq) = o_var(ineq)
              endif
            enddo

          enddo

        enddo

      enddo

      return
      end
