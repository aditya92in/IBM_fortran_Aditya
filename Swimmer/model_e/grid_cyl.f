      program grid_cyl
c
      implicit none
c
      integer imax, jmax
c
      parameter (imax=60)
      parameter (jmax=60)
c
      real      xh(0:imax+1)
      real      dx(0:imax+1)
      real      yh(0:jmax+1)
      real      dy(0:jmax+1)
c
c
c  z-direction:
c
c
c   0------0.5------1.---------2.2---------3.4
c   H               H                       H
c
c y-direction:
c
c
c   0------------1.5---------2.5-----------4.
c                  H          H
c
      real dxstart
      parameter (dxstart=0.1)
      integer imina, imaxa, i, j
      integer ibegin, iend
      real    xbegin, xend
c
      xh = 0.
      dx = 0.
      yh = 0.
      dy = 0.
c
      xbegin =  1. 
      xend   =  0.0
      ibegin = 6
      iend   = 0
      imina  = 0
      imaxa  = jmax + 1
      call grstr(xbegin,xend,ibegin,iend,
     & dxstart,imina,imaxa,yh)
c
c
      xbegin =  1.0
      xend   =  3.5
      ibegin = 6
      iend   = 33
      imina  = 0
      imaxa  = jmax + 1
      call grstr(xbegin,xend,ibegin,iend,
     & dxstart,imina,imaxa,yh)
c
c
      xbegin = 6.0
      xend   = 3.5
      ibegin = 50
      iend   = 33
      imina  = 0
      imaxa  = jmax + 1
      call grstr(xbegin,xend,ibegin,iend,
     & dxstart,imina,imaxa,yh)
c
c
      xbegin = 6.0
      xend   = 8.0
      ibegin = 50
      iend   = 60
      imina  = 0
      imaxa  = jmax + 1
      call grstr(xbegin,xend,ibegin,iend,
     & dxstart,imina,imaxa,yh)
c
      yh(jmax+1) = yh(jmax) + (yh(jmax) - yh(jmax-1))
c
      do j=1,jmax+1
         dy(j) = yh(j) - yh(j-1)
      enddo
c
      open(unit=10,file='gridcube_j_60_005.dat')
      dy(0) = dy(1)
      do j=0,jmax+1
         write(6,1000)j,yh(j),dy(j)
         write(10,*)j,yh(j),dy(j)
      enddo
      close(10)
c
1000  format(i4,2(f16.8))
c
      xbegin = 1.0
      xend   = 0.0
      ibegin = 6
      iend   = 0
      imina  = 0
      imaxa  = imax + 1
      call grstr(xbegin,xend,ibegin,iend,
     & dxstart,imina,imaxa,xh)
c
c
      xbegin = 1.0
      xend   = 4.0
      ibegin = 6
      iend   = 29
      imina  = 0
      imaxa  = imax + 1
      call grstr(xbegin,xend,ibegin,iend,
     & dxstart,imina,imaxa,xh)
c
c
      xbegin = 7.0
      xend   = 4.0
      ibegin = 52
      iend   = 29
      imina  = 0
      imaxa  = imax + 1
      call grstr(xbegin,xend,ibegin,iend,
     & dxstart,imina,imaxa,xh)
c
c
      xbegin = 7.0
      xend   = 16.0
      ibegin = 52
      iend   = 60
      imina  = 0
      imaxa  = imax + 1
      call grstr(xbegin,xend,ibegin,iend,
     & dxstart,imina,imaxa,xh)
c
c
      xh(imax+1) = xh(imax) + (xh(imax) - xh(imax-1))
c
      do i=1,imax+1
         dx(i) = xh(i) - xh(i-1)
      enddo
c
      dx(0) = dx(imax)
      open(unit=10,file='gridcube_i_60_005.dat')
      do i=0,imax+1
         write(6,1000)i,xh(i),dx(i)
         write(10,*)i,xh(i),dx(i)
      enddo
      close(10)

      stop
      end

      subroutine grstr(xbegin,xend,ibegin,iend,
     & dxstar,imin,imax,xarray)
c
c
c purpose: give a geometric grid distribution in one direction
c          for a piece of calculation domain which extends from
c          xbegin to xend (cell boundaries), where cell boundaries
c          xbegin has index i and cell boundary xend has index iend
c          
c output:  array xarray will be filled from i = ibegin 
c                                      to   i = iend
c          with a grid space distribution ranging from
c                                      x(ibegin) = xbegin
c
c                                                 to
c
c                                      x(iend)   = xend
c         
c  ibegin and iend do not have to be in ascending order, but
c  sign(ibegin - iend) must be equal to sign(xbegin - xend)
c
c
c
           implicit none
           integer i,ibegin,iend
           integer imax,imin
c
c eps = tolerance in iteration process for determining the
c       stretch parameter in the geometric distribution
c
c CAUTION: this parameter may have to be changed on using 
c          double or single precision
c          
c          
           real eps
c
c single precision:
c
           parameter (eps = 1.e-14)
c
c double precision:
c          parameter (eps = 1.e-12)
c
           integer iturna, npoint,  niter
           real xbegin, xend, xarray(imin:imax)
           real dxarray(imin:imax)
           real xlengt, strech, dxstar
c
c take care, that xend can be SMALLER then xleft:
c
	   if(iend-ibegin .lt.0)then
              iturna = -1
           else
	      iturna =  1
           endif
c
c determine length to be covered
c
           xlengt = abs(xend - xbegin)  - abs(dxstar)
           if(abs(xend - xbegin).lt.abs(dxstar)+1.e-6)stop 'CRAZY GRID'
           write(6,*)'xlengt, dxstar',xlengt, dxstar
c
c determine number of grid points in the length to be covered
c
           npoint = abs(iend - ibegin   ) - 1
           write(6,*)'npoint = ', npoint
c
c this is the number of cells covered
c
c use the geometric series summation formula to obtain the stretch
c ratio between the successive cells
c
c
c  the distribution will be:
c  
c  dxstar  | strech   | stretch**2  |  strech**3   |... | stretch**n |
c
c          ^                                                         ^
c          |                                                         |
c
c      boundary with xbegin                               boundary end
c
c                                        n
c  xlengt = dxstar  * stretch*   stretch    - 1
c                                ---------------
c                                stretch    - 1
c
c
c                                        n
c  xlengt             stretch*   stretch    - 1
c  ------- =                     ---------------
c  dxstar                        stretch    - 1
c
c    
c use bisection method for finding the root to  the equation
c
      call strsub(xlengt, dxstar, npoint, strech, niter)
      write(6,*)'strech strsub   = ',strech
c                                        n
c
           write(6,*)'niter, strech, xlengt, npoint', niter, strech,
     &     xlengt, npoint, ibegin, iend
c
c strech determined, calculate grid distribution
c
      xarray(ibegin) = xbegin
      xarray(ibegin+1*iturna) = xbegin + iturna * dxstar
c
c
      do i=2,npoint+1
         xarray(ibegin + iturna*i) = 
     &                               xarray(ibegin + iturna*i -iturna) 
     &                             + iturna*dxstar*(strech**(i-1))
      enddo
      write(6,*)'i, ibegin ,iend,xbegin,xend,iturna,dxstar,strech',
     &i, ibegin ,iend,xbegin,xend,iturna,dxstar,strech
c
C     write(6,*)'imax, ib, iend',imax, ib, iend
C     write(6,*)'xarray(iend) = ', xarray(iend)
C     xarray(iend+1) = xarray(iend) + (xarray(iend) - xarray(iend-1))
C     write(6,*)'xarray(iend+1) = ', xarray(iend+1)
C     do i=imin+1,iend! +1
C        dxarray(i) = xarray(i)-xarray(i-1)
C     enddo
C        dxarray(imin) = dxarray(imin+1)
      return
      end
c
      subroutine strsub(xlengt, dxstar, n, strech, niter)
c
c solve the equation for the stretch factor strech in an iterative way
c using the bisection method. Note, that we use two different formulas
c for the  cases of a grid with increasing grid cells (grid cells become
c smaller going away from the first grid cell) and the case of 
c increasing grid cells (grid cells become larger going away from the
c first grid cell)
c
      implicit none
c
      real    xlengt, dxstar
      integer n
      integer niter
      real  xleft, xmid, xrigh
      real  fleft, fmid, frigh
      real  eps,  strech
      real    A
c
c A = ratio length line segment : initial cell length
c if A > n then the grid cells must increase
c if A < n then the grid cells must decrease
c
      A = xlengt / dxstar
c
      eps = 1.e-14
c
c number of necessary iterations
c
      niter = 0
c
c check for an equidistant grid,resulting in strech = 1.
c
      if(abs(A - n).le.eps)then
         strech = 1.
         niter  = 1
         return
      endif
c
c increasing grid cells
c
      if(A.gt.n)then
c
c limit the interval where the solution must be found
c
c strech must be > 1 and smaller than some (big) number
c
      xleft  = 1.000001
      xrigh  = 10.
c
      fleft  = xleft - (1. + A*(xleft - 1.)/xleft)**(1./n)
      frigh  = xrigh - (1. + A*(xrigh - 1.)/xrigh)**(1./n)
c
      write(6,*) 'A gt n', fleft, frigh, A, n
c
      if(abs(sign(1.,fleft) - sign(1.,frigh)).gt.1.e-3)then
100   continue
         niter = niter + 1
         xmid = 0.5*(xleft + xrigh)
         fmid = xmid  - (1. + A*(xmid  - 1.)/xmid )**(1./n)
         if(abs(sign(1.,fmid) - sign(1.,fleft)).lt.1.e-3)then
            xleft = xmid
         else
            xrigh = xmid
         endif
         write(6,*)niter, xmid
         if(abs((xleft - xrigh)).lt.eps)then
            write(6,*)'niter, ROOT = ', niter, xmid
            goto 400
         endif
         goto 100
      endif
c
c decreasing grid cells
c
        elseif(A.le.n)then
c
c
c limit the interval where the solution must be found
c
c strech must be < 1 and < 0
c
      xleft  = 0.
      xrigh  = 0.99999
      fleft = xleft**(n+1) - xleft - A*(xleft-1)
      frigh = xrigh**(n+1) - xrigh - A*(xrigh-1)
c
      write(6,*) 'A lt n', fleft, frigh, A, n
c
      if(abs(sign(1.,fleft) - sign(1.,frigh)).gt.1.e-3)then
200   continue
         niter = niter + 1
         xmid = 0.5*(xleft + xrigh)
         fmid = xmid **(n+1) - xmid  - A*(xmid -1)
         if(abs(sign(1.,fmid) - sign(1.,fleft)).lt.1.e-3)then
            xleft = xmid
         else
            xrigh = xmid
         endif
         write(6,*)niter, xmid
         if(abs((xleft - xrigh)).lt.eps)then
            write(6,*)'niter, ROOT = ', niter, xmid
            goto 400
         endif
         goto 200
      endif
c    
      endif
400   continue
      strech = xmid
      write(6,*)'strech strsub routine = ',strech
      return
      end
