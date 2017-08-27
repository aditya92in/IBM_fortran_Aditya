      subroutine reagr2(lx,ly)
c
c version: DOES NOT USE ghost points for HO scheme
c          ALSO not for normal velocities
c
      implicit none
c
      include 'param.txt'
      include 'common.txt'
      real    lx,ly,lz
      integer idummy
      integer iequid
      real    factor
      real    dtmp
      real    zstart, dzstart, rlength
c
c     dx(0)  dx(1)                dx(imax)   dx(imax+1)
c  |       |       |             |        |             |
c  xh(0) xh(1)                    xh(imax)    xh(i1)
c
c periodic:  dx(0) = dx(imax)
c            dx(imax+1) = dx(1)
c            this is accomplished by making dx(1) = dx(imax)
c            and using the standard boundary treatment:
c
c            dx(0) = dx(1) etc.
c
c
c iadd = 1: NO ghost points used for HO scheme
c      = 0:    ghost points used for HO scheme
      integer iadd
      integer iadd2
      iadd = 0
      iadd2= 0
c
      xh  = 0.
      dx  = 0.
      dxh = 0.

      yh  = 0.
      dy  = 0.
      dyh = 0.

      dtmp = 15./imax
      xh(0)= -5
      do i=1,i1
         xh(i) = xh(i-1) + dtmp
      enddo
c
      dtmp = 1./jmax
      yh(0)= 0.
      do j=1,j1
         yh(j) = yh(j-1) + dtmp
      enddo
!     open(unit=10,file='gridcube_i_60_005.dat')
!     do i=0,imax
!        read(10,*)idummy,xh(i)
!        write(6,*)'i, xh(i) ',i, xh(i)
!     enddo
!     close(10)
!     open(unit=10,file='gridcube_j_60_005.dat')
!     do j=0,jmax
!        read(10,*)idummy,yh(j)
!        write(6,*)'j, yh(j) ',j, yh(j)
!     enddo
!     close(10)
c ATT: periodic!!!
      yh(j1) = yh(jmax) + yh(1) - yh(0)
      xh(i1) = xh(imax) + xh(1) - xh(0)
c
      do i=imax+2,imax+ih
         xh(i) = xh(i-1) + (xh(i-1)-xh(i-2))
      enddo
c
      do i=-1,1-ih,-1
         xh(i) = xh(i+1) - (xh(i+2) - xh(i+1))
      enddo

      rlength = yh(jmax) - yh(0)

      do j=1-jh,-1
         yh(j) = yh(j + jmax) - rlength
      enddo
      do j=jmax+2,jmax+jh
         yh(j) = yh(j-jmax) + rlength
      enddo

c form all other quantities from zh, yh, xh ..
c
      do i=2-ih,imax+ih
         dx(i) = xh(i) - xh(i-1)
      enddo
c
CATT      dx(0) = dx(1)
CNONPERIODIC      dx(1-ih) = dx(1-ih+imax)
      dx(1-ih) = dx(1-ih+1)
c
      do i=1-ih,imax+ih
         xp(i) = xh(i) - 0.5*dx(i)
      enddo
c
      do i=1-ih,imax+ih-1
         dxh(i) = xp(i+1) - xp(i)
      enddo
c
      dxh(imax+ih)   = dxh(ih)
c
      do i=1-ih,imax+ih
         dxi(i)  = 1./dx(i)
         dxhi(i) = 1./dxh(i)
      enddo
      write(6,*)'dx :',dx(1) , dx(imax+1)
C      if(abs(dx(1) - dx(imax+1)).gt.1.e-8)STOP 'dx !!!'
C      if(abs(dx(ibl1-1) - dx(ibl1)).gt.1.e-8)then
C      write(6,*)dx(ibl1-1),dx(ibl1)
C      STOP 'dx1 !!!'
C      endif
C      if(abs(dx(ibl2+1) - dx(ibl2)).gt.1.e-8)STOP 'dx2 !!!'
C      if(abs(dxh(ibl2) - dx(ibl2)).gt.1.e-8)then
C      write(6,*) 'dx3 !!!'
C      write(6,*)ibl2, dxh(ibl2),dx(ibl2)
C      stop
C      endif
C      if(abs(dxh(ibl1-1) - dx(ibl1)).gt.1.e-8)then
C      write(6,*) 'dx4 !!!'
C      write(6,*) ibl1, dxh(ibl1-1) , dx(ibl1), dx(ibl1-1)
C      stop
C      endif
C      if(abs(dxh(0) - dxh(imax)).gt.1.e-8)STOP 'dx5 !!!'
c
c y --
c
      do j=2-jh,jmax+jh
         dy(j) = yh(j) - yh(j-1)
      enddo
c
CATTT      dy(0) = dy(1)
      dy(1-jh) = dy(1-jh+jmax)
c
      do j=1-jh,jmax+jh
         yp(j) = yh(j) - 0.5*dy(j)
      enddo
c
      do j=1-jh,jmax+jh-1
         dyh(j) = yp(j+1) - yp(j)
      enddo
c
      dyh(jmax+jh)   = dyh(jh)
c
      do j=1-jh,jmax+jh
         dyi(j)  = 1./dy(j)
         dyhi(j) = 1./dyh(j)
      enddo
         
      do i=1-ih,imax+ih
         write(6,1000)'i, dx, dxi, dxh, dxhi, ',
     & i, dx(i), dxi(i), dxh(i), dxhi(i)
      enddo
      do j=1-jh,jmax+jh
         write(6,1000)'j, dy, dyi, dyh, dyhi ',
     & j, dy(j), dyi(j), dyh(j), dyhi(j)
      enddo
      write(6,*)
      do i=1-ih,imax+ih
         write(6,*)'i, xp(i) ', i, xp(i)
      enddo
      write(6,*)
      do j=1-jh,jmax+jh
         write(6,*)'j, yp(j) ', j, yp(j)
      enddo
c
      Rp = yp
      Rv = yh
c
      Rpi = 1/Rp
      do j=1-jh,jmax+jh
         if (j .ne. 0) Rvi(j) = 1/Rv(j)
      enddo
      Rvi(0) = 0.
c
1000  format(a,i4,5(f16.6))
c
c make additional help quantities for higher order Veldman
c discretisation
c
c dxi, dyi and dzi for larger volume
c
         write(6,*)'X4 quantities'
      do i=1,imax
         dx4W (i)  = xp(i)-xp(i-3)
         dx4E (i)  = xp(i+3)-xp(i)
         dxh4W(i)  = xh(i)-xh(i-3)
         dxh4E(i)  = xh(i+3)-xh(i)

         dxi4W (i) = 1./dx4W (i)
         dxi4E (i) = 1./dx4E (i)
         dxhi4W(i) = 1./dxh4W(i)
         dxhi4E(i) = 1./dxh4E(i)

         write(6,5000)dxi4W (i),dxi4E (i),dxhi4W(i),dxhi4E(i)
      enddo
5000  format (5(f16.8))
c
      do j=1,jmax
         dy4S (j)  = yp(j)-yp(j-3)
         dy4N (j)  = yp(j+3)-yp(j)
         dyh4S(j)  = yh(j)-yh(j-3)
         dyh4N(j)  = yh(j+3)-yh(j)

         dyi4S (j) = 1./dy4S (j)
         dyi4N (j) = 1./dy4N (j)
         dyhi4S(j) = 1./dyh4S(j)
         dyhi4N(j) = 1./dyh4N(j)

      enddo
c
      ordaru = 0.
      ordarv = 0.
c
      return
      ENd
