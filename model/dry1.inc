************************************************************************

      parameter ( nn=10000,ntp=20000,ny=410,nx=410)

************************************************************************

      implicit real*8(a-h,o-z)
      
      common/m/      grav, amax, epsh

      common/reals/  rain

      common/reals/  t_pond, start, finish

      common/reals/  tmax, dt, t, t_rain, beta, 
     &               dz

      common/reals/  dxdum, r8mV, r8etaV, r8alphaV,
     &               r8mB, r8etaB, r8alphaB

      common/array/  xc(0:ny,0:nx), yc(0:ny,0:nx), zc(0:ny,0:nx), 
     &               dzz(0:ny,0:nx,2), xnc(0:ny,0:nx), 
     &               vegc(0:ny,0:nx)

      common/array/  tt(0:ntp)
      
      common/array/  f(0:ny,0:nx,3,2), q(0:ny,0:nx,3), fdum(3), qs(3),
     &               xflux0(0:ny,0:nx), xflux1(0:ny,0:nx),
     &               yflux0(0:ny,0:nx), yflux1(0:ny,0:nx)
    
      common/array/  dh(0:ny,0:nx,2), du(0:ny,0:nx,2), dv(0:ny,0:nx,2)
      
      common/array/  itype(ny,nx,4), ipos(ny,nx,4), inum(ny,nx),
     &               nop(ny,nx,4), kbeg(ny), kend(ny)

      common/array/  sn(0:ny,0:nx,2), cn(0:ny,0:nx,2), 
     &               ds(0:ny,0:nx,2), dxi(0:ny,0:nx,2),
     &               deta(0:ny,0:nx,2), area(0:ny,0:nx), 
     &               zinflmap(0:ny,0:nx), zinflmap2(0:ny,0:nx)

      common/srctrm/ fix(ny,nx,3),
     &               sx(0:ny,0:nx), sy(0:ny,0:nx)

      common/ints/   ncol, np, ne, nbcell, ilim,
     &               nprt, it, itp, nt, iprt,
     &               nrow, isurf

      common/ints/   imodel, isveg

      common/array/  h(0:ny,0:nx), u(0:ny,0:nx), v(0:ny,0:nx),
     &               hn(nn), un(nn), vn(nn),
     &               hp(0:ny,0:nx), up(0:ny,0:nx), vp(0:ny,0:nx)

      common/array/  delta_theta(2), r8Ksat(2), H_i(2), Ao(2)
