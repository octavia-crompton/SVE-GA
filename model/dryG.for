C      program swe2d
C====================================================================
C     2D SVE solver coupled to Green Ampt
*********************************************************************
      include 'dry.inc'    

C     Input files
      open(2,file = 'input/coords.dat')
      open(3,file = 'input/boundary.dat')     
      open(4,file = 'input/params.dat') 
      open(5,file = 'input/veg.dat')
      open(11,file = 'input/nodes.dat')        
C       open(304,file = 'input/vanG.dat')

C     Output files: saved every nprt time steps (dt_p seconds)
      open(100,file = 'output/h.out')                
      open(101,file = 'output/time.out')        
      open(102,file = 'summary.txt')  
             
      open(104,file = 'output/hydro.out')
      open(110,file = 'output/fluxes1234.out')   ! boundary fluxes      
      open(111,file = 'output/dvol.out') ! SVE volume tracking              

C     Read input, set up grid, initial conditions.
      call cpu_time(start)
      call init
C       call initvg

C     write headings for fortran output files
      write(101,*)  "     time (s) |   CFL  "
      write(104,*)  "     time (s) |  hydrograph (m3/s)" 
C     heading for volume tracking file (dvol.out)
      write(111,* )  "time (s)      vol               zflux         ",
     &                "     zinfl "
C     heading for horizontal flux tracking file (fluxes1234.out)
      write(110,203)  "t", "flux1", "flux2", "flux3", "flux4", 
     &                "fluxin", "hydro"

      write(100, *) "   j   k    h              V               U ",
     &    "           zinflmap2        xflux0         yflux0",
     &    "         xflux1          yflux1"

C     Initialize boundary fluxes (m2/s)
      flux1 = 0.d0
      flux2 = 0.d0
      flux3 = 0.d0
      flux4 = 0.d0
      fluxin = 0.d0
      hydro = 0.d0
C     Initialize gridlevel fluxes (m2/s)      
      xflux0 = 0.0
      yflux0 = 0.0
      xflux1 = 0.0
      yflux1 = 0.0

C     Track infiltration volume
      zinfl = 0.d0  
      zinflmap2 = 0.d0
      zinflmap = 0.d0

      t_pond = 1e10        !  time of ponding  (placeholder value)
  

      dvol = 0.d0
      zflux = 0.d0

      vol = 0.D0
      do j=1,nrow 
        do k=kbeg(j),kend(j)
          vol = vol + h(j,k)*area(j,k)
        enddo
      enddo 
      
C     Write ICs to output files
C    fluxes1234.out
      write(110,202) t, flux1, flux2, flux3, flux4, fluxin, hydro 

      write(111,200) t, dvol, zflux, zinfl  ! dvol.out
      call myoutput
      
C    Begin time loop.
      do it = 0,nt-1
        t = t + dt
        if (t .gt. t_rain) rain = 0. 
        amax = 0.d0
    
C       Compute predictor.
        do j=1,nrow
          do k=kbeg(j),kend(j)
            call bconds(j,k,h,u,v)                              
            call predict(j,k)
          enddo
        enddo
        
C Loop over cells to compute fluxes.
        do j=1,nrow
        do k=kbeg(j),kend(j)
          call bconds(j,k,hp,up,vp)
          call fluxes(j-1,j,k,k,1)          ! horizontal faces.   
          call fluxes(j,j,k-1,k,2)          ! vertical faces.
          
          do l=1,inum(j,k)       
            if(ipos(j,k,l) .eq. 3) then
              call fluxes(j,j,k,k+1,2)       ! right boundaries. 
                flux3 =  flux3 + f(j,k+1,1,2)*ds(j,k+1,2)*dt
                hydro =  hydro + f(j,k,1,2)*ds(j,k,2)*dt
            
            elseif(ipos(j,k,l) .eq. 2) then
              call fluxes(j,j+1,k,k,1)         ! lower boundaries.
                flux2 = flux2  + f(j+1,k,1,1)*ds(j+1,k,1)*dt                
            
            elseif(ipos(j,k,l) .eq. 1) then         ! left boundary
              flux1 = flux1 + f(j,k,1,2)*ds(j,k,2)*dt
              fluxin = fluxin + f(j,k+1,1,2)*ds(j,k+1,2)*dt
            
            elseif(ipos(j,k,l) .eq. 4) then         ! top boundary
              flux4 = flux4 + f(j,k,1,1)*ds(j,k,1)*dt              
            endif
          enddo

C         Increment the fluxes at each grid cell  
          yflux0(j,k) = yflux0(j,k) + f(j,k,1,1)*ds(j,k,1)*dt        
          yflux1(j,k) = yflux1(j,k) + f(j+1,k,1,1)*ds(j+1,k,1)*dt   
          xflux0(j,k) = xflux0(j,k) + f(j,k,1,2)*ds(j,k,2)*dt   
          xflux1(j,k) = xflux1(j,k) + f(j,k+1,1,2)*ds(j,k+1,2)*dt       
        enddo
        enddo
      

C Compute corrector solution.
        do j=1,nrow
        do k=kbeg(j),kend(j)
          call source(j,k,hp(j,k),up(j,k),vp(j,k), 1)

C         Increment source (p-i) volume   
          zinfl = zinfl - dt*qs(1)*area(j,k) + dt*rain*area(j,k)

C         Compute infiltration for this grid cell and timestep. 
C          dzinfl, zinflmap :   units : m^3 
          dzinfl = - dt*qs(1)*area(j,k) + dt*rain*area(j,k) 
          zinflmap2(j,k) = zinflmap2(j,k) + dzinfl   
          zinflmap(j,k) = zinflmap(j,k) + dzinfl 
          
         do l=1,3
            q(j,k,l) = q(j,k,l) + dt/area(j,k)*(
     &        f(j,k,l,2)*ds(j,k,2) + f(j,k,l,1)*ds(j,k,1) -  
     &        f(j+1,k,l,1)*ds(j+1,k,1) - f(j,k+1,l,2)*ds(j,k+1,2)) 
     &        + dt*qs(l)
          enddo
       
        enddo
        enddo
C Store solution. 
C Solve continuity equation even in all cells, but 
C solve momentum equations in wet cells only.
        do j=1,nrow
          do k=kbeg(j),kend(j)
C Check for negative depth.
            if(q(j,k,1) .ge. 0.D0) then
              h(j,k) = q(j,k,1)
            else
              zinfl = zinfl + q(j,k,1)*area(j,k)*dt
              zinflmap2(j,k) = zinflmap2(j,k) + q(j,k,1)*area(j,k)*dt
              zinflmap(j,k) = zinflmap(j,k) + q(j,k,1)*area(j,k)*dt
              q(j,k,1) = 0.D0
              h(j,k) = 0.D0
            endif
C Neglect momentum in nearly dry cells.
           if(h(j,k) .lt. epsh) then  
              u(j,k) = 0.d0
              v(j,k) = 0.d0
              do l=2,3
                q(j,k,l) = 0.d0
              enddo
            elseif (h(j,k) .ge. epsh) then  
              u(j,k) = q(j,k,2)/h(j,k)
              v(j,k) = q(j,k,3)/h(j,k)
            endif
           enddo
       enddo 

        nprthydro = 1/dt 
        if (mod(it, nprthydro) .eq. nprthydro-1) then
          
          vol0 = vol
          vol = 0.D0
          do j=1,nrow 
            do k=kbeg(j),kend(j)
              vol = vol + h(j,k)*area(j,k)
            enddo
          enddo 
        
          dvol = vol - vol0
C         Sum of all lateral fluxes at boundaries
          zflux =   flux2 + flux3 - flux4  - flux1
C         Fluxes are positive  out of the domain
          write(111,200) t, dvol, zflux, zinfl ! dvol.out
          write(110,202) t, flux1, flux2, flux3, flux4, fluxin, hydro 

C           write(110,202) t, flux1, flux2, flux3, flux4  !fluxes1234.out    

          flux1 = 0.d0
          flux2 = 0.d0
          flux3 = 0.d0
          flux4 = 0.d0
          fluxin = 0.d0
          hydro = 0.d0
          zinfl = 0.d0
          
          call myhydro
          
        endif 

C      Below here only executed every nprt time steps        
        if (mod(it, nprt) .eq. nprt-1)  then         

C         Increment print count  
          itp = itp + 1     
          write(101, 204) t, amax*dt 

C         'output/h.out' 
          do j=1,nrow
            do k=kbeg(j),kend(j)
              write(100, 201) j, k, h(j,k), u(j,k), v(j,k),
     &                  zinflmap2(j,k), xflux0(j,k), yflux0(j,k),
     &                  xflux1(j,k), yflux1(j,k)
              zinflmap2(j,k) = 0.d0       
            enddo
          enddo
          write(100, 205) itp, t
  
          xflux0 = 0.0
          yflux0 = 0.0
          xflux1 = 0.0
          yflux1 = 0.0
          
        endif
         
         r8minvol = 0.1d0*epsh*dxdum**2*nrow*ncol
         if ((vol .le. r8minvol) .and. (t .gt. t_rain)) then
           call gracefulExit  
           write(102,*) 'exit: dry'
           return
        endif

        if(amax*dt.gt. 10.0d0) then
          call gracefulExit 
          write(*,*) 'AMAX too big', t, amax*dt
          write(102,*) 'exit: AMAX'
          return  ! leave early          
        
        endif
        if(t .gt. tmax) then
          call gracefulExit
          write(102,*) 'exit: t>tmax'
          return   ! leave early          
        endif

      enddo
      call gracefulExit 
      write(102,*) 'exit: finished'
      stop
      
 200  format(' ', f8.1, 5e19.8)  ! for writing to dvol.out
 201  format(' ', 2i4, 8e15.6) ! format h.out

 
 202  format(' ', 7e19.8)  ! for writing fluxes to fluxes1234.out
 203  format(' ', 7A19 ) ! format fluxes
 205  format(' ', i8, f9.2)
 
 204  format(' ', f10.2, f15.5)
      end

************************************************************************
      subroutine gracefulExit
      include 'dry.inc'
      
      call cpu_time(finish)
      write(102,*) '  ' 
      write(102,*) 'runtime:', (finish-start)
      write(102,*) 'ponding: ', t_pond
      write(102,*) 'final_time: ', t
      
      return 
      
 211  format(' ', A14, f12.2) 
 212  format(' ', A14, f12.3, A2)  
      end

************************************************************************
      subroutine myhydro
      include 'dry.inc'
      
      hydro = 0.d0
      do j=1,nrow
            hydro = hydro + f(j,ncol,1,2)*ds(j,ncol,2)  ! m3/s
      enddo
      write(104, 204)  t, hydro      
      
      return
 202  format(' ', i8, f9.2)
 204  format(' ', f15.2 , E20.10) 
      end
************************************************************************
      subroutine myoutput
      include 'dry.inc'
  
C     'output/time.out' 
      write(101, 204) t, amax*dt 

C     'output/h.out' 
      do j=1,nrow
        do k=kbeg(j),kend(j)
          write(100, 201) j, k, h(j,k), u(j,k), v(j,k),
     &                  zinflmap2(j,k), xflux0(j,k), yflux0(j,k),
     &                  xflux1(j,k), yflux1(j,k)
          zinflmap2(j,k) = 0.d0       
        enddo
      enddo
      write(100, 202) itp, t
      
      return
 202  format(' ', i8, f9.2)
 201  format(' ', 2i4, 8e15.6) ! format h.out
 204  format(' ', f10.2, f15.5)
      end
************************************************************************
      subroutine source(j,k,depth,udum,vdum,lstep)
      !   
      !   Input:
      !         j, k : grid cell
      !         depth  (h(j,k), hp(j,k))
      !         udum  (u(j,k), up(j,k))
      !         vdum  (v(j,k), vp(j,k)) 
      !   Modifies common variables: 
      !         qs(3):
      !              qs(1) = winflt = p-i  
      !
      include 'dry.inc'      
      
      real ( kind = 8 )   PI, depth, fm, feta, falpha
      real ( kind = 8 )   psi, delta_theta, F0, r8Ks, sorptivity, denom 

      
      F0 = abs(zinflmap(j,k)/area(j,k))

      ! determine whether the cell is vegetated or not
      ! define the roughness parameters (fm, falpha, feta)
      if (vegc(j,k) .gt. 0) then
        isveg = 1
        fm = r8mV
        falpha = r8alphaV
        feta = r8etaV
        d_theta = delta_theta(1)
        r8Ks = r8Ksat(1)
        psi = abs(H_i(1))
        sorptivity = Ao(1)
      else
        isveg = 2 
        fm = r8mB
        falpha = r8alphaB 
        feta = r8etaB  
        d_theta = delta_theta(2)
        r8Ks = r8Ksat(2)
        psi = abs(H_i(2))
        sorptivity = Ao(2)
      endif  

      if (imodel .eq. 2) then 
        denom = (sorptivity**2+4*r8Ks*F0)**0.5 - sorptivity 
        PI = r8Ks + r8Ks*sorptivity/denom

      elseif (imodel .eq. 1) then
        PI = r8Ks*(d_theta*psi + F0 )/F0
C       if (d_theta .lt. 1e-7) then ! See if this is necessary
C         PI = r8Ks
      endif
      
      PI = min(depth/dt, PI)

C     only update depth in corrector step.
      if (lstep .eq. 1)  then  

        if  (r8Ks .le. 1e-10) then   
           winflt = rain
C       Case 1: rain and no ponding
        elseif ((depth .le. 1e-8) .and. (rain .gt. 0.d0)) then  

          if (PI .lt. rain) then
            if (t_pond .gt. 1e5) then
                t_pond = t 
            endif
            winflt = rain - PI

          else    ! no ponding,  flux = rain
            winflt = 0.d0 
          endif  
C       Case 2: no ponding and no rain
        elseif ((depth.le.1e-9).and.(rain .le. 1e-10)) then  
          winflt = 0.d0
C       Case 3:  ponding  (rain or no rain)
        elseif (depth .gt. 0.d0) then  
          winflt = rain - PI
        endif
C     Predictor step:        
      else
        winfl = 0.d0
      endif 

      if (depth .gt. epsh) then 
        
        vmag = dsqrt(udum*udum + vdum*vdum) 

        if ((feta .eq. 0.5) .and. (fm .eq. 0.d0)) then  ! cylinder

          fricSx = (falpha)**(2.d0)*udum*vmag 
          fricSy = (falpha)**(2.d0)*vdum*vmag 

        elseif (feta .eq. 0.5) then  ! most schemes
          
           fricSx = (falpha/depth**fm)**(2.d0)*udum*vmag 
           fricSy = (falpha/depth**fm)**(2.d0)*vdum*vmag 

        elseif (feta .eq. 1) then  ! special case for poisseuille 
        
          fricSx = falpha*udum/depth**fm 
          fricSy = falpha*vdum/depth**fm 
          Rel = vmag*depth/1.e-6   

          if (Rel .gt. 500) then 
              ffact = 0.5  ! okay for smooth surfaces
              fricSx = ffact*udum*vmag/8./grav/depth
              fricSy = ffact*vdum*vmag/8./grav/depth      
          endif

        endif

         qs(1) = winflt
         qs(2) = 0.5D0*udum*winflt - grav*depth*fricSx -  
     &          grav*depth*sx(j,k)
         qs(3) = 0.5D0*vdum*winflt - grav*depth*fricSy - 
     &          grav*depth*sy(j,k)
       else
          
         qs(1) = winflt
         qs(2) = 0.d0
         qs(3) = 0.d0
       endif
  
      return
 209  format(' ', ' time of ponding is =',f8.1,'s ')           

      end
      
************************************************************************
      subroutine fluxes(jl,jr,kl,kr,i1)
      include 'dry.inc'
      !   
      !   Input:
      !             jl, jr, kl, kr : left/right grid cell indices
      !             i1 :  interface type (horizontal or vertical)
      
      !   Modifies f(jr,kr,1:3,1:2)
                  
C     MUSCL extrapolation at cell interface.
      hl = hp(jl,kl) + 0.5D0*dh(jl,kl,i1)
      ul = up(jl,kl) + 0.5D0*du(jl,kl,i1)
      vl = vp(jl,kl) + 0.5D0*dv(jl,kl,i1)
      hr = hp(jr,kr) - 0.5D0*dh(jr,kr,i1)
      ur = up(jr,kr) - 0.5D0*du(jr,kr,i1)
      vr = vp(jr,kr) - 0.5D0*dv(jr,kr,i1)
      snn = sn(jr,kr,i1)
      cnn = cn(jr,kr,i1)
      
      if(i1 .eq. 1) then
        dx =  deta(jr,kr,2)*area(jr,kr) 
        dy = - deta(jr,kr,1)*area(jr,kr)
      else
        dx = - dxi(jr,kr,2)*area(jr,kr) 
        dy =   dxi(jr,kr,1)*area(jr,kr)
      endif
C Needed for dry bed problems.
      if(hl .lt. 0.D0) hl = 0.D0
      if(hr .lt. 0.D0) hr = 0.D0
C Compute arithmatic averages for source terms.
      havg = 0.5D0*(hl + hr)
      uavg = 0.5D0*(ul + ur)
      vavg = 0.5D0*(vl + vr)
C Prevent leakage into cells with higher bed elevation.
      etal = hp(jl,kl) + zc(jl,kl)
      etar = hp(jr,kr) + zc(jr,kr)
C Fluxes and source terms.
      if(havg .le. 0.D0 ) then
        do i=1,3  
          f(jr,kr,i,i1) = 0.D0        
        enddo
      else   
        call solver(hl,hr,ul,ur,vl,vr,fdum,snn,cnn,dx,dy)
        do i=1,3
          f(jr,kr,i,i1) = fdum(i)
        enddo
      endif

      return
      end 
      
************************************************************************
      subroutine solver(hl,hr,ul,ur,vl,vr,ff,sndum,cndum,dx,dy)        
      !
      !   Input:
      !             hl,hr,ul,ur,vl,vr : left/right h,u,v
      !             i1 :  face numbers
    
      !   Output: 
      !             f(3)  --> fdum in fluxes

      include 'dry.inc'
      !implicit real*8(a-h,o-z)
      dimension ws(3), e(3,3), a(3), astar(3), da(3), ff(3), dum(3)
      !common/m/ grav, amax, epsh

C Compute Roe averages at cell face.
      duml  = dsqrt(hl)
      dumr  = dsqrt(hr)
      hhat  = duml*dumr
      uhat  = (duml*ul + dumr*ur)/(duml + dumr)
      vhat  = (duml*vl + dumr*vr)/(duml + dumr)
      chat  = dsqrt(0.5D0*grav*(hl + hr))
      uperp = uhat*cndum + vhat*sndum
C Compute eigenvalues.  Lambdahat
      a(1) = uperp - chat
      a(2) = uperp
      a(3) = uperp + chat
C Compute approximate wave strengths.
      dhdum   = hr - hl
      dudum   = ur - ul
      dvdum   = vr - vl
      dupar = -dudum*sndum + dvdum*cndum
      duperp=  dudum*cndum + dvdum*sndum
      !  deltaVhat (eqn 31 in 1999)
      ws(1) = 0.5D0*(dhdum - hhat*duperp/chat)
      ws(2) = hhat*dupar     
      ws(3) = 0.5D0*(dhdum + hhat*duperp/chat)
C Compute right eigenvectors.  Rhat
      e(1,1) = 1.D0
      e(2,1) = uhat - chat*cndum   
      e(3,1) = vhat - chat*sndum
      e(1,2) = 0.D0
      e(2,2) = -sndum
      e(3,2) =  cndum
      e(1,3) = 1.D0
      e(2,3) = uhat + chat*cndum
      e(3,3) = vhat + chat*sndum
C Entropy fix.
      dl = dsqrt(dx*dx + dy*dy)
      cl = dsqrt(grav*hl)
      cr = dsqrt(grav*hr)
      uperpl = ul*cndum + vl*sndum
      uperpr = ur*cndum + vr*sndum
      al1 = uperpl - cl
      al3 = uperpl + cl
      ar1 = uperpr - cr
      ar3 = uperpr + cr
      da(1) = dmax1(0.D0, 4.D0*(ar1 - al1))
      da(2) = 0.d0
      da(3) = dmax1(0.D0, 4.D0*(ar3 - al3))
      do i=1,3   
         if(dabs(a(i)) .lt. 0.5D0*da(i)) then
           astar(i) = a(i)*a(i)/da(i) + 0.25D0*da(i)  
          else
           astar(i) = dabs(a(i))
         endif
         if(astar(i)/dl .gt. amax) amax = astar(i)/dl
      enddo
C Compute flux increments.
      do i=1,3
        dum(i) = 0.D0
        do l=1,3
          dum(i) = dum(i) + (astar(l)*ws(l))*e(i,l) 
        enddo
      enddo
C Add flux to appropriate cell faces.
      ff(1) = 0.5D0*(f1(hl,uperpl) + f1(hr,uperpr) - dum(1))
      ff(2) = 0.5D0*(f2(hl,ul,uperpl,cndum) 
     &            + f2(hr,ur,uperpr,cndum) - dum(2))
      ff(3) = 0.5D0*(f3(hl,vl,uperpl,sndum)   
     &            + f3(hr,vr,uperpr,sndum) - dum(3))

      return
      end
      
************************************************************************
      subroutine predict(j,k)
      include 'dry.inc'
      
      !   Predictor step 
      !   Input:
      !             j, k : grid cell
      
      !   Modifies common variables:
      !             dh, du, dv  
      !             qs(3)  -  source term
      !             hp, up, vp  - predictor variables

      do kk=1,2                  ! loop over coord. directons.
        if(kk .eq. 1) then
          jr = j + 1
          jl = j - 1
          kr = k
          kl = k
        else
          jr = j
          jl = j
          kr = k + 1
          kl = k - 1
        endif
C Compute gradients only in wet cells.
        if(h(j,k) .ge. epsh) then
C Limit afree surface elevation to reduce dissipation..
          dh1 = h(j,k) + zc(j,k) - h(jl,kl) - zc(jl,kl)
          dh2 = h(jr,kr) + zc(jr,kr) - h(j,k) - zc(j,k)
C Needed to minimize dissipation at wet/dry interfaces.
          if(h(jl,kl) .lt. epsh) dh1 = 2.d0*dh1 
          if(h(jr,kr) .lt. epsh) dh2 = 2.d0*dh2
          call limitr(ilim,beta,dh1,dh2,dhh)
          dh(j,k,kk) = dhh - dzz(j,k,kk)
C U velocity.
          du1 = u(j,k) - u(jl,kl)
          du2 = u(jr,kr) - u(j,k)
          call limitr(ilim,beta,du1,du2,duu)
          du(j,k,kk) = duu
C V velocity.
          dv1 = v(j,k) - v(jl,kl)
          dv2 = v(jr,kr) - v(j,k)
          call limitr(ilim,beta,dv1,dv2,dvv) 
          dv(j,k,kk) = dvv
        else
          dh(j,k,kk) = 0.D0
          du(j,k,kk) = 0.D0
          dv(j,k,kk) = 0.D0
        endif
      enddo
C Generalized velocities.
      uxi = u(j,k)*dxi(j,k,1) + v(j,k)*dxi(j,k,2)
      ueta = u(j,k)*deta(j,k,1) + v(j,k)*deta(j,k,2)
C Predictor.
!     Call source before updating the predictor valuess
      call source(j, k, h(j,k), u(j,k), v(j,k), 0)
      if(h(j,k) .ge. epsh**(0.75d0)) then           
        qs(2) = qs(2)/h(j,k)
        qs(3) = qs(3)/h(j,k)
      else
        qs(2) = 0.d0
        qs(3) = 0.d0
      endif
      hp(j,k) = h(j,k) - 0.5D0*dt*(
     &    uxi*dh(j,k,1) + h(j,k)*(dxi(j,k,1)*du(j,k,1) + 
     &                            dxi(j,k,2)*dv(j,k,1)) +   
     &    ueta*dh(j,k,2) + h(j,k)*(deta(j,k,1)*du(j,k,2) +
     &                             deta(j,k,2)*dv(j,k,2)) + qs(1))   
      up(j,k) = u(j,k) - 0.5D0*dt*(      
     &    grav*dxi(j,k,1)*dh(j,k,1) + uxi*du(j,k,1) +
     &    grav*deta(j,k,1)*dh(j,k,2) + ueta*du(j,k,2) + qs(2))
      vp(j,k) = v(j,k) - 0.5D0*dt*(
     &    grav*dxi(j,k,2)*dh(j,k,1) + uxi*dv(j,k,1) +
     &    grav*deta(j,k,2)*dh(j,k,2) + ueta*dv(j,k,2) + qs(3))
C Correct any negative depths.
      if(hp(j,k) .lt. 0.d0) then
        hp(j,k) = 0.d0
        dh(j,k,1) = 0.d0
        dh(j,k,2) = 0.d0
      endif
C Neglect momentum in nearly dry cells.
      if(hp(j,k) .le. epsh) then
        up(j,k) = 0.D0
        vp(j,k) = 0.D0
        do i=1,2
          du(j,k,i) = 0.d0
          dv(j,k,i) = 0.d0
        enddo
      endif

      return
      end
      
************************************************************************
      subroutine bconds(j,k,hdum,udum,vdum)
      include 'dry.inc'

      dimension hdum(0:ny,0:nx), udum(0:ny,0:nx), vdum(0:ny,0:nx)
C     Loop over all boundary faces in the cell. 
      
      do i=1, inum(j,k)
        if(ipos(j,k,i) .eq. 1) then ! front face.         
            jj = j
            kk = k-1
            jl = j
            kl = k 
            j2 = j
            k2 = k+1
            io = 2 
        elseif(ipos(j,k,i) .eq. 2) then ! right face.
            jj = j+1
            kk = k 
            jl = j+1 
            kl = k  
            j2 = j-1
            k2 = k
            io = 1
        elseif(ipos(j,k,i) .eq. 3) then ! back face.
            jj = j
            kk = k+1
            jl = j
            kl = k+1
            j2 = j
            k2 = k-1 
            io = 2 
        elseif(ipos(j,k,i) .eq. 4) then ! left face.
            jj = j-1
            kk = k 
            jl = j
            kl = k 
            j2 = j+1
            k2 = k  
            io = 1  
        endif
C   Open boundary.
        if(itype(j,k,i) .eq. 0)then
          
          dh(jj,kk,io) = dh(j,k,io)
          du(jj,kk,io) = du(j,k,io)
          dv(jj,kk,io) = dv(j,k,io)
          
          hdum(jj,kk) = 2.D0*hdum(j,k) - hdum(j2,k2)
          udum(jj,kk) = 2.D0*udum(j,k) - udum(j2,k2)
          vdum(jj,kk) = 2.D0*vdum(j,k) - vdum(j2,k2)
                   
          
C Wall boundary.
       elseif(itype(j,k,i) .eq. 1)then
          dh(jj,kk,io) = dh(j,k,io)
          du(jj,kk,io) = du(j,k,io)
          dv(jj,kk,io) = dv(j,k,io)
          
          hdum(jj,kk) = 2.D0*hdum(j,k) - hdum(j2,k2)

          udum(jj,kk) = udum(j,k)*(sn(jl,kl,io)*sn(jl,kl,io) - 
     &                                 cn(jl,kl,io)*cn(jl,kl,io)) -
     &                   2.D0*vdum(j,k)*sn(jl,kl,io)*cn(jl,kl,io)
          vdum(jj,kk) = vdum(j,k)*(cn(jl,kl,io)*cn(jl,kl,io) - 
     &                      sn(jl,kl,io)*sn(jl,kl,io)) -
     &                      2.D0*udum(j,k)*sn(jl,kl,io)*cn(jl,kl,io)
        
C Specified depth and velocity (supercritical).
     
       elseif(itype(j,k,i) .eq. 2)then
          
          dh(jj,kk,io) = 0.D0
          du(jj,kk,io) = 0.D0
          dv(jj,kk,io) = 0.D0
          hdum(jj,kk) = fix(j,k,1) 
          udum(jj,kk) = fix(j,k,2)
          vdum(jj,kk) = fix(j,k,3)
          if(isurf .eq. 1) hdum(jj,kk) = hdum(jj,kk) - zc(j,k)
      
C Specified flow rate (subcritical).
        elseif(itype(j,k,i) .eq. 4)then 
          du(jj,kk,io) =  0.d0
          dv(jj,kk,io) =  0.d0 
          if(hdum(j,k) .ge. epsh) then
            hdum(jj,kk) = hdum(j,k) 
            dh(jj,kk,io) =  0.d0
           elseif(hdum(j,k)*hdum(j2,k2) .ge. epsh) then
            hdum(jj,kk) = 2.D0*hdum(j,k) - hdum(j2,k2)
            dh(jj,kk,io) = dh(j,k,io)
          endif
          if(hdum(jj,kk) .ge. epsh) then
            udum(jj,kk) = fix(j,k,2)/hdum(jj,kk)
            vdum(jj,kk) = fix(j,k,3)/hdum(jj,kk)
          else
            udum(jj,kk) = 0.D0
            vdum(jj,kk) = 0.D0
          endif        
        endif
        
        if ((itype(j,k,i) .eq. 4) .and. (t .gt. t_rain)) then 
          itype(j,k,i) = 1
        endif  

        if (hdum(jj,kk)  .lt. 0.D0) hdum(jj,kk) = 0
        
        if(hdum(jj,kk) .lt. epsh) then
          udum(jj,kk) = 0.D0
          vdum(jj,kk) = 0.D0
        endif
        
      enddo

      return
      end
      
***********************************************************************
      subroutine grid
      include 'dry.inc'
      !     Read grid data from  'coords.dat', 'veg.dat' and 'nodes.dat'
      !       file 2 = 'input/coords.dat'
      !       file 5 = 'input/veg.dat'
      !       file 11 = 'input/nodes.dat'

      real ( kind = 8 ) x(nn), y(nn), zz(nn)
      real ( kind = 8 ) veg(nn)

      area = 0.d0
      read(2,*) np, ne
      if(np .gt. nn) then
        write(*,*) np
        write(*,*) 'ERROR:  parameter nn in file dry.inc is too small'
        stop
      endif
      do i=1,np
        read(2,*) x(i), y(i), zz(i)
      enddo
      do i=1,np
        read(5,*) veg(i)
      enddo 
      do j=1,nrow
        do k=kbeg(j),kend(j)
          read(11,*) (nop(j,k,i), i=1,4)
        enddo
      enddo
      
C Compute grid metrics.
      do j=1,nrow
      do k=kbeg(j),kend(j)
        n1 = nop(j,k,1)
        n2 = nop(j,k,2)
        n3 = nop(j,k,3)
        n4 = nop(j,k,4)
        xc(j,k) = 0.25D0*(x(n1) + x(n2) + x(n3) + x(n4))
        yc(j,k) = 0.25D0*(y(n1) + y(n2) + y(n3) + y(n4))
        vegc(j,k) = veg(n1) 
        if (vegc(j,k) .gt. 0) then
          xnc(j,k) = r8alphaV
        else
          xnc(j,k) = r8alphaB
        endif        
        dxdxi = 0.5D0*(-x(n1) + x(n2) + x(n3) - x(n4))
        dxdeta = 0.5D0*(-x(n1) - x(n2) + x(n3) + x(n4))
        dydxi = 0.5D0*(-y(n1) + y(n2) + y(n3) - y(n4))
        dydeta = 0.5D0*(-y(n1) - y(n2) + y(n3) + y(n4))
        area(j,k) = dxdxi*dydeta - dxdeta*dydxi
        if(area(j,k) .le. 0.D0)then
           write(*,*) 'area error in cell ',j,k
           stop
        endif
        dxi(j,k,1) = dydeta/area(j,k)
        deta(j,k,1) = -dydxi/area(j,k)

        dxi(j,k,2) = -dxdeta/area(j,k)
        deta(j,k,2) = dxdxi/area(j,k)
        sx(j,k)=((zz(n2)-zz(n4))*(y(n3)-y(n1))-(zz(n3)-zz(n1))
     &     *(y(n2)-y(n4)))/(2.D0*area(j,k))
        sy(j,k)=((zz(n3)-zz(n1))*(x(n2)-x(n4))-(zz(n2)-zz(n4))
     &     *(x(n3)-x(n1)))/(2.D0*area(j,k))
        zc(j,k) = 0.25D0*(zz(n1) + zz(n2) + zz(n3) + zz(n4))
        dzz(j,k,1) = sx(j,k)*dxdxi + sy(j,k)*dydxi
        dzz(j,k,2) = sx(j,k)*dxdeta + sy(j,k)*dydeta
      enddo
      enddo
      dxdum = dxdxi          
C Compute cell face angles.
      do j=1,nrow
      do k=kbeg(j),kend(j)
       ddx = x(nop(j,k,2)) - x(nop(j,k,1))
       ddy = y(nop(j,k,2)) - y(nop(j,k,1))
       ds(j,k,2) = dsqrt(ddx*ddx + ddy*ddy)
       sn(j,k,2) = ddx/ds(j,k,2)            ! Horizontal face.
       cn(j,k,2) = -ddy/ds(j,k,2)          ! Horizontal face.
       ddx = x(nop(j,k,4)) - x(nop(j,k,1))
       ddy = y(nop(j,k,4)) - y(nop(j,k,1))
       ds(j,k,1) = dsqrt(ddx*ddx + ddy*ddy)
       sn(j,k,1) = -ddx/ds(j,k,1)            ! Vertical face.
       cn(j,k,1) = ddy/ds(j,k,1)
       do i=1,inum(j,k)
        if(ipos(j,k,i) .eq. 3) then
           ddx = x(nop(j,k,3)) - x(nop(j,k,4))
           ddy = y(nop(j,k,3)) - y(nop(j,k,4))
           ds(j,k+1,2) = dsqrt(ddx*ddx + ddy*ddy)
           sn(j,k+1,2) = ddx/ds(j,k+1,2)     ! Top (boundary) faces.
           cn(j,k+1,2) = -ddy/ds(j,k+1,2)
         elseif(ipos(j,k,i) .eq. 2) then
           ddx = x(nop(j,k,3)) - x(nop(j,k,2))
           ddy = y(nop(j,k,3)) - y(nop(j,k,2))
           ds(j+1,k,1) = dsqrt(ddx*ddx + ddy*ddy)
           sn(j+1,k,1) = -ddx/ds(j+1,k,1)     ! Right (boundary) faces.
           cn(j+1,k,1) = ddy/ds(j+1,k,1)
        endif
       enddo
      enddo
      enddo

C Set some things in ghost cells.
      do j=1,nrow
      do k=kbeg(j),kend(j) 
        do i=1, inum(j,k)
          call findbc(i,j,k,jj,kk,j2,k2)
          area(jj,kk) = area(j,k)
          sx(jj,kk) = sx(j,k)
          sy(jj,kk) = sy(j,k)
          dxi(jj,kk,1) =dxi(j,k,1)
          deta(jj,kk,1) = deta(j,k,1)
          dxi(jj,kk,2) = dxi(j,k,2)
          deta(jj,kk,2) = deta(j,k,2)
          xc(jj,kk) = 2.D0*xc(j,k) - xc(j2,k2)
          yc(jj,kk) = 2.D0*yc(j,k) - yc(j2,k2)
          zc(jj,kk) = 2.D0*zc(j,k) - zc(j2,k2)
          xnc(jj,kk) = 2.D0*xnc(j,k) - xnc(j2,k2)
          vegc(jj,kk) = vegc(j,k)
        enddo
      enddo
      enddo


      return
      
 210  format(' ', A14, f12.2)
      end 
      
************************************************************************
      subroutine findbc(i,j,k,jj,kk,j2,k2)
      include 'dry.inc'
      if(ipos(j,k,i) .eq. 1) then  
         jj = j
         kk = k-1
         j2 = j
         k2 = k+1
      elseif(ipos(j,k,i) .eq. 2) then
         jj = j+1
         kk = k 
         j2 = j-1
         k2 = k
       elseif(ipos(j,k,i) .eq. 3) then 
         jj = j
         kk = k+1
         j2 = j
         k2 = k-1
       elseif(ipos(j,k,i) .eq. 4) then 
         jj = j-1
         kk = k
         j2 = j+1
         k2 = k 
      endif
      return
      end   
************************************************************************
      subroutine limitr(i,beta,dq1,dq2,dq)
C               
      implicit real*8(a-h,o-z)
C Lax-Wendroff.
      if(i .eq. 1)then
        dq = dq2 
C Beam-Warming.
      elseif(i .eq. 2)then
        dq = dq1 
C Fromm
      elseif(i .eq. 3)then
        dq = 0.5D0*(dq1 + dq2) 
C Double minmod.
      elseif(i .eq. 4)then
        a = 0.5D0*(dq1 + dq2)
        b = 2.D0*dq1
        c = 2.D0*dq2
        if(a*b .gt. 0.D0 .and. b*c .gt. 0.D0) then
           dq = fmin1(a, b, c)  
         else
           dq = 0.D0
        endif
C Beta Family.
      else
        if(dq1*dq2 .le. 0.D0) then
           dq = 0.D0
         else
           dq = fmin2(fmax2(dq1, dq2), beta*fmin2(dq1, dq2))
        endif
      endif

      return     
      end      
************************************************************************
      subroutine init
             
      include 'dry.inc' 
      ! Read  file 'input/params.dat'.
      read(4,'(a72)') dum
      read(4,*) grav, dt
      read(4,'(a72)') dum
      read(4,*) tmax, t_rain
      read(4,'(a72)') dum
      read(4,*) rain, nt
      read(4,'(a72)') dum      
      read(4,*) epsh, beta
      read(4,'(a72)') dum
      read(4,*) nprt 
      read(4,'(a72)') dum
      read(4,*) h0, u0, v0
      read(4,'(a72)') dum
      read(4,*) r8mV,r8etaV, r8alphaV  
      read(4,'(a72)') dum
      read(4,*) r8mB, r8etaB, r8alphaB

      read(4,'(a72)') dum
      read(4,*) imodel 

      read(4,'(a72)') dum
      read(4,'(a72)') dum
      read(4,*)  r8Ksat(1), H_i(1), delta_theta(1), Ao(1)

      read(4,'(a72)') dum
      read(4,'(a72)') dum
      read(4,*) r8Ksat(2), H_i(2), delta_theta(2), Ao(2)
 

      ! Read  file 'input/boundary.dat'.
      read(3,'(a72)') dum
      read(3,*) nbcell
      read(3,'(a72)') dum
      do ii=1,nbcell
         read(3,*)j,k,inum(j,k),
     &   (itype(j,k,i),i=1,inum(j,k)),(ipos(j,k,i),i=1,inum(j,k))
      enddo
      read(3,'(a72)') dum
      read(3,*) nrow
      read(3,'(a72)') dum
      read(3,*) ncol
      read(3,'(a72)') dum
      
      do j=1,nrow       
         read(3,*) idum, kbeg(j), kend(j)
      enddo
      kbeg(nrow+1) = kbeg(nrow)
      kend(nrow+1) = kend(nrow)   
      read(3,'(a72)') dum
      read(3,*) ndir

      read(3,'(a72)') dum
      do i=1,ndir
        read(3,*) j,k,fix(j,k,1),fix(j,k,2),fix(j,k,3)
C         write(*,*) j,k,fix(j,k,1),fix(j,k,2),fix(j,k,3)
      enddo

      isurf = 2    ! 2 to set flow depth
      ilim =  5   

C Set up grid.
      call grid
      
C Set initial conditions.    
      t = 0.D0
      do j=1,nrow
        do k=kbeg(j),kend(j)
          
          n1 = nop(j,k,1)
          n2 = nop(j,k,2)
          n3 = nop(j,k,3)
          n4 = nop(j,k,4)
          
          if(isurf .eq. 1) then
            h(j,k) = h0 - zc(j,k)
           else
            h(j,k) = h0
          endif
          
          u(j,k) = u0
          v(j,k) = v0
        enddo
      enddo     

      do j=1,nrow
      do k=kbeg(j),kend(j)
C       if (h(j,k) .lt. 0.D0) h(j,k) = 0.D0
        q(j,k,1) = h(j,k)
        q(j,k,2) = h(j,k)*u(j,k)
        q(j,k,3) = h(j,k)*v(j,k)
        
        do i=1,inum(j,k)
C For fixed flux BCs, must initialize depth.
          if(itype(j,k,i) .eq. 4) then
            
            call findbc(i,j,k,jj,kk,j2,k2)
            if(ipos(j,k,i) .eq. 1 .or. ipos(j,k,i) .eq. 3) then
              qflux = fix(j,k,2)*cn(j,k,2) + fix(j,k,3)*sn(j,k,2)

              dx = -dxi(j,k,2)*area(j,k)
              dy = dxi(j,k,1)*area(j,k)  
              dss = dsqrt(dx*dx + dy*dy)

              if(h(j,k) .lt. epsh) then
                 if(qflux*dzz(j,k,2).lt.0.D0) then
                   qflux = dabs(qflux)
                   hnorm=(qflux*xnc(j,k)/dsqrt(dabs(dzz(j,k,2)/dss)))
     &                      **(3.D0/5.D0)     

C                    write(*,211) ' normal depth = ', hnorm,
C      &                              'm specified'
                   write(102,*)  'hnorm: ', i,j,k,ipos(j,k,i), hnorm     
                   h(jj,kk) = hnorm 
                   hp(jj,kk) = hnorm     
                else
C                   write(*,*)  qflux, dzz(j,k,2)
                  write(*,*) 'adverse slope or zero Manning n in cell'
                  write(*,*) 'flux boundary ', j,k,ipos(j,k,i)
                  !read(*,*) h(jj,kk)                  
                  h(jj,kk) = hnorm                         
                  hp(jj,kk) = hnorm
                endif
              endif
              qflux = dabs(qflux)
              
            endif
          endif
        enddo
        
      enddo
      enddo
      
!     Initialize common time variables 
!     Make this every second instead
!       read(303,*) nt
!         do i=0,nt
!           read(303,*) tt(i), rain(i)
!         enddo
!
      
      
      if (ntp .lt. nt) then 
        write(*,*) 'nt too small'
        write(*,*) 'ntp = ', ntp
        write(*,*) 'nt = ', nt                
        stop
      endif
                   


      return
 210  format(' ', A14, f12.3)
 211  format(' ', A, e12.3, A) 
 212  format(' ', A14, f12.3, A2)  
      end
      
************************************************************************
      real*8 function f1(h,up)
      implicit real*8(a-h,o-z)
      f1 = h*up
      return
      end

************************************************************************
      real*8 function f2(h,u,up,cn)
      implicit real*8(a-h,o-z)
      common/m/ grav, amax, epsh
      f2 = h*u*up + 0.5D0*grav*h*h*cn
      return
      end

************************************************************************
      real*8 function f3(h,v,up,sn)
      implicit real*8(a-h,o-z)
      common/m/ grav, amax, epsh
      f3 = h*v*up + 0.5D0*grav*h*h*sn
      return
      end
************************************************************************
      real*8 function fmin1(a,b,c)
      implicit real*8(a-h,o-z)
      if (a .lt. 0.D0) then
         fmin1 = -1.D0*dmin1(dabs(a),dabs(b),dabs(c)) 
        else
         fmin1 = dmin1(a,b,c) 
      endif 
      return
      end
************************************************************************
      real*8 function fmin2(a,b)
      implicit real*8(a-h,o-z)
      if (a .lt. 0.D0) then
          fmin2 = -dmin1(dabs(a),dabs(b))
         else
          fmin2 = dmin1(a,b)
      endif 
      return
      end
************************************************************************
      real*8 function fmax2(a,b)
      implicit real*8(a-h,o-z)
      if (a .lt. 0.D0) then
          fmax2 = -dmax1(dabs(a),dabs(b))
         else
          fmax2 = dmax1(a,b) 
      endif
      return
      end

************************************************************************
      subroutine fdiag(aa, bb)
    ! Make a diagonal 2D matrix (hh) from 1D array (dh)
    ! 
    ! Input:
    !             dh  (real,)
    !             gh (real, kind = 8)    -  dimension
    ! Output:
    !             gtheta  (real, kind= 8) - volumetric moisture content
    !             gK  (real, kind= 8) - hydraulic conductivity
      include 'dry.inc'
      
      real(kind=8), dimension(nz) :: aa
      real(kind=8), dimension(nz, nz) :: bb
      
      bb = 0.0
      do k=1,nz
        bb(k,k) = aa(k)
      enddo

      return
      end

