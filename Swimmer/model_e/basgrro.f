      subroutine basgrro(lx,ly,lz,Re)
c
c version: DOES NOT USE ghost points for HO scheme
c          ALSO not for normal velocities
c
      implicit none
c
      include 'param.txt'
      include 'common.txt'
      include 'block.inc'
      real    lx,ly,lz,Re
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
      iadd = 1
      iadd2= 1
c
      xh  = 0.
      dx  = 0.
      dxh = 0.

      yh  = 0.
      dy  = 0.
      dyh = 0.

      zh  = 0.
      dz  = 0.
      dzh = 0.
c
      dtmp = 1./kmax
      zh(0)= 0.
      do k=1,k1
         zh(k) = zh(k-1) + dtmp
      enddo
      rlength = 1.
      zstart  = 0.
      call basgr(zp, zh, dz, dzh, dzi, dzhi, kmax,kh,rlength,zstart,
     &                 dzstart)
c
C     close(33)
      zh(k1) = 2*zh(kmax) - zh(kmax-1)
c
      dtmp = lx/imax
      xh(0)= 0.
      do i=1,i1
         xh(i) = xh(i-1) + dtmp
      enddo
       xh(i1) = xh(imax) + (xh(imax)-xh(imax-1))
c
      dtmp = ly/jmax
      yh(0)= 0.
      do j=1,j1
         yh(j) = yh(j-1) + dtmp
      enddo
c
      yh(j1) = yh(jmax) + yh(1) - yh(0)
c
c ATT: periodic!!!
      xh(i1) = xh(imax) + xh(1) - xh(0)
c
C     close(33)
c
c periodic continuation for x and y
c
      rlength = xh(imax) - xh(0)
c
      do i=1-ih,-1
         xh(i) = xh(i + imax) - rlength
      enddo
      do i=imax+2,imax+ih
         xh(i) = xh(i-imax) + rlength
      enddo
c
      rlength = yh(jmax) - yh(0)

      do j=1-jh,-1
         yh(j) = yh(j + jmax) - rlength
      enddo
      do j=jmax+2,jmax+jh
         yh(j) = yh(j-jmax) + rlength
      enddo

c
c equidistant continuation for z
c
      dtmp = zh(1) - zh(0)
      do k=-1,1-kh,-1
         zh(k) = zh(k+1) - dtmp
      enddo

      dtmp = zh(kmax) - zh(kmax-1)
      do k=kmax+1,kmax+kh
         zh(k) = zh(k-1) + dtmp
      enddo

c
c form all other quantities from zh, yh, xh ..
c
      do i=2-ih,imax+ih
         dx(i) = xh(i) - xh(i-1)
      enddo
c
CATT      dx(0) = dx(1)
      dx(1-ih) = dx(1-ih+imax)
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
c y --
c
      do j=2-jh,jmax+jh
         dy(j) = yh(j) - yh(j-1)
      enddo
c
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
c
      do k=2-kh,kmax+kh
         dz(k) = zh(k) - zh(k-1)
      enddo
c
      dz(1-kh) = dz(kh)
c
      do k=1-kh,kmax+kh
         zp(k) = zh(k) - 0.5*dz(k)
      enddo
c
      do k=1-kh,kmax+kh-1
         dzh(k) = zp(k+1) - zp(k)
      enddo
c
      dzh(kmax+kh)   = dzh(kmax+kh-1)
c
      do k=1-kh,kmax+1
         dzi(k)  = 1./dz(k)
         dzhi(k) = 1./dzh(k)
      enddo
c
c interpolation factors 
c
c
      do i=0,imax
         fx(i)  = dx(i+1)/(dx(i) + dx(i+1))
      enddo
         fx(i1) = fx(imax)
c
      do j=0,jmax
         fy(j)  = dy(j+1)/(dy(j) + dy(j+1))
      enddo
         fy(j1) = fy(jmax)
c
c
      do k=0,kmax
         fz(k)  = dz(k+1)/(dz(k) + dz(k+1))
      enddo
         fz(k1) = fz(kmax)
c
c ATT adapt grid for vonvegence of solver !!!!!
c
c
c periodicity:
c
!     dx(1)  = dx(imax)
!     dxi(1) = dxi(imax)
c
c block walls:
c
!     dx(ibl1-1) = dx(ibl1)
!     dx(ibl2+1) = dx(ibl2)
c
c
!        dx(0)  = dx(1)
!        dxi(0) = dxi(1)

!        dx(i1)  = dx(imax)
!        dxi(i1) = dxi(imax)
c
c
c periodicity:
c
!     dy(1) = dy(jmax)
!     dyi(1) = dyi(jmax)
c
c
c block walls:
c
!     dy(jbl1-1) = dy(jbl1)
!     dy(jbl2+1) = dy(jbl2)
c
c
!        dy(0)  = dy(1)
!        dyi(0) = dyi(1)

!        dy(j1)  = dy(jmax)
!        dyi(j1) = dyi(jmax)
c
c block walls:
c
!     dz(kbl1-1) = dz(kbl1)
!     dz(kbl2+1) = dz(kbl2)
c
c
!        dz(0)  = dz(1)
!        dzi(0) = dzi(1)
c
!        dz(k1)  = dz(kmax)
!        dzi(k1) = dzi(kmax)
c
         
      do i=1-ih,imax+ih
         write(6,1000)'i, dx, dxi, dxh, dxhi, fx ',
     & i, dx(i), dxi(i), dxh(i), dxhi(i), fx(i)
      enddo
      do j=1-jh,jmax+jh
         write(6,1000)'j, dy, dyi, dyh, dyhi, fy ',
     & j, dy(j), dyi(j), dyh(j), dyhi(j), fy(j)
      enddo
      do k=1-kh,kmax+kh
         write(6,2000)'k, dz, dzi, dzh, dzhi, fz, z+ ',
     & k, dz(k), dzi(k), dzh(k), dzhi(k), fz(k), zh(k)*Re
      enddo
      write(6,*)
      do i=1-ih,imax+ih
         write(6,*)'i, xp(i) ', i, xp(i)
      enddo
      write(6,*)
      do j=1-jh,jmax+jh
         write(6,*)'j, yp(j) ', j, yp(j)
      enddo
      write(6,*)
      do k=1-kh,kmax+kh
         write(6,*)'k, zp(k) ', k, zp(k)
      enddo
c
1000  format(a,i4,6(f16.6))
2000  format(a,i4,7(f16.6))
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

         dx4   (i) = dx(i-1)  + dx(i)  + dx(i+1)
         dxh4  (i) = dxh(i-1) + dxh(i) + dxh(i+1)
         dxi4  (i) = 1./dx4   (i)
         dxhi4 (i) = 1./dxh4  (i)
         write(6,5000)dxi4W (i),dxi4E (i),dxhi4W(i),dxhi4E(i),
     &             dxi4  (i),dxhi4 (i)
      enddo
5000  format (6(f16.8))
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

         dy4   (j) = dy(j-1)  + dy(j)  + dy(j+1)
         dyh4  (j) = dyh(j-1) + dyh(j) + dyh(j+1)
         dyi4  (j) = 1./dy4   (j)
         dyhi4 (j) = 1./dyh4  (j)
      enddo
c
      do k=1,kmax
         dz4B (k)  = zp(k)-zp(k-3)
         dz4T (k)  = zp(k+3)-zp(k)
         dzh4B(k)  = zh(k)-zh(k-3)
         dzh4T(k)  = zh(k+3)-zh(k)

         dzi4B (k) = 1./dz4B (k)
         dzi4T (k) = 1./dz4T (k)
         dzhi4B(k) = 1./dzh4B(k)
         dzhi4T(k) = 1./dzh4T(k)

         dz4   (k) = dz(k-1)  + dz(k)  + dz(k+1)
         dzh4  (k) = dzh(k-1) + dzh(k) + dzh(k+1)
         dzi4  (k) = 1./dz4   (k)
         dzhi4 (k) = 1./dzh4  (k)
      enddo
c
c interpolation factors for larger volume, for staggered and
c non-staggered variables 
c
c NOTES:
c
c because we shift only over one i-, j- or k-index each time we go to
c another point, we need two kinds of interpolation factors, namely
c fx4W and fx4E, because otherwise we have to jump over index ranges
c with size 4. Moreover, we need fx4W and fxh4W, for interpolation 
c in staggered- or non-staggered direction
c
c there is another kind of interpolation factor too: we need it to
c interpolate V, W velocities to the right position of a U cel face.
c In the 2nd order case this kind of term is absent.
c
         write(6,*)'X4 quantities'
      do i=1,imax
         fx4W(i)  = ( 0.5*dx (i) + dx (i-1) )/
     &             ( 0.5*dx (i) + dx (i-1) + dx (i-2) + 0.5*dx (i-3))
         fx4E(i)  = ( dx (i+2) + 0.5*dx (i+3) )/
     &             ( 0.5*dx (i) + dx (i+1) + dx (i+2) + 0.5*dx (i+3))
         fxh4W(i) = ( 0.5*dxh(i) + dxh(i-1) )/
     &             ( 0.5*dxh(i) + dxh(i-1) + dxh(i-2) + 0.5*dxh(i-3))
         fxh4E(i) = ( dxh(i+2) + 0.5*dxh(i+3) )/
     &             ( 0.5*dxh(i) + dxh(i+1) + dxh(i+2) + 0.5*dxh(i+3))
         fxe4(i)  = ( 0.5*dxh(i) + dxh(i+1) )/
     &              ( dxh(i) + dxh(i+1) + dx (i-1) )
         write(6,4000)fx4W(i), fx4E(i), fxh4W(i), fxh4E(i), fxe4(i)
      enddo
4000  format (5(f16.8))
c
      do j=1,jmax
         fy4S(j)  = ( 0.5*dy (j) + dy (j-1) )/
     &             ( 0.5*dy (j) + dy (j-1) + dy (j-2) + 0.5*dy (j-3))
         fy4N(j)  = ( dy (j+2) + 0.5*dy (j+3) )/
     &             ( 0.5*dy (j) + dy (j+1) + dy (j+2) + 0.5*dy (j+3))
         fyh4S(j) = ( 0.5*dyh(j) + dyh(j-1) )/
     &             ( 0.5*dyh(j) + dyh(j-1) + dyh(j-2) + 0.5*dyh(j-3))
         fyh4N(j) = ( dyh(j+2) + 0.5*dyh(j+3) )/
     &             ( 0.5*dyh(j) + dyh(j+1) + dyh(j+2) + 0.5*dyh(j+3))
         fye4(j)  = ( 0.5*dyh(j) + dyh(j+1) )/
     &              ( dyh(j) + dyh(j+1) + dy (j-1) )
      enddo
c
      do k=1,kmax
         fz4B(k)  = ( 0.5*dz (k) + dz (k-1) )/
     &             ( 0.5*dz (k) + dz (k-1) + dz (k-2) + 0.5*dz (k-3))
         fz4T(k)  = ( dz (k+2) + 0.5*dz (k+3) )/
     &             ( 0.5*dz (k) + dz (k+1) + dz (k+2) + 0.5*dz (k+3))
         fzh4B(k) = ( 0.5*dzh(k) + dzh(k-1) )/
     &             ( 0.5*dzh(k) + dzh(k-1) + dzh(k-2) + 0.5*dzh(k-3))
         fzh4T(k) = ( dzh(k+2) + 0.5*dzh(k+3) )/
     &             ( 0.5*dzh(k) + dzh(k+1) + dzh(k+2) + 0.5*dzh(k+3))
         fze4(k)  = ( 0.5*dzh(k) + dzh(k+1) )/
     &              ( dzh(k) + dzh(k+1) + dz (k-1) )
      enddo
c
c ATT set coefficients equal to 0.5!!
c
      fxe4 = 0.5
      fye4 = 0.5
      fze4 = 0.5
c
      ordaru = 1.
      ordarv = 1.
      ordarw = 1.
c
c
      do k=1-kh,2+iadd2
         do j=1-jh,jmax+jh
         do i=1-ih,imax+ih
            ordarw(i,j,k) = 0.
         enddo
         enddo
      enddo
c
      do k=1-kh,2+iadd
         do j=1-jh,jmax+jh
         do i=1-ih,imax+ih
            ordaru(i,j,k) = 0.
            ordarv(i,j,k) = 0.
         enddo
         enddo
      enddo
c
      do k=kmax-2-iadd2,kmax+kh
         do j=1-jh,jmax+jh
         do i=1-ih,imax+ih
            ordarw(i,j,k) = 0.
         enddo
         enddo
      enddo
c
      do k=kmax-1-iadd,kmax+kh
         do j=1-jh,jmax+jh
         do i=1-ih,imax+ih
            ordaru(i,j,k) = 0.
            ordarv(i,j,k) = 0.
         enddo
         enddo
      enddo
c
C     ordaru = 0.
C     ordarv = 0.
C     ordarw = 0.
c
      return
      ENd
      subroutine basgr(z, zh, dz, dzh,dzi,dzhi,kmax,kh,rlength,zstart,
     &                 dzstart)
c
c
      implicit none
c
      integer k,  kmax, kh
c
      real z(1-kh:kmax+kh)
      real zh(1-kh:kmax+kh)
      real dz(1-kh:kmax+kh)
      real dzh(1-kh:kmax+kh)
      real dzi(1-kh:kmax+kh)
      real dzhi(1-kh:kmax+kh)
      real rlength
      real zstart
      real dzstart ! dz start for other grid
c
c help variables
c
      real pi, expk, zhelp, aconst, atanh, ztmp, power
      parameter (power = 1.3) ! 1.1
c
      integer myid, ierr
c
      pi = 4.*atan(1.)
c
C     kmax = k1 - 1
c
c
c first distribution zh
c
      do k=0, kmax/2
        zh(k) = (1.*k)**power
      enddo
      ztmp = zh(kmax/2)
      do k=0, kmax/2
        zh(k) = rlength*0.5*zh(k)/ztmp
      enddo
      do k= kmax/2+1, kmax
         zh(k) = rlength-zh(kmax-k)
      enddo
c
c
      zh(kmax+1) = zh(kmax) + (zh(kmax) - zh(kmax-1))
c
c add zstart
c
      do k=0,kmax+1
c        zh(k) = zh(k) + zstart
      enddo
c
      do k=1,kmax+1
         z(k) = zh(k-1) + 0.5*(zh(k) - zh(k-1))
      enddo
c
      z(0) =  z(1) -2.*(z(1) - zh(0))
c
c
      do k=1,kmax+1
         dz(k) = zh(k) - zh(k-1)
      enddo
c
c
      dz(0) = dz(1)
c
c
      do k=0,kmax
         dzh(k) = z(k+1) - z(k)
      enddo
c
c
      dzh(kmax+1) = dzh(kmax)
c
c
c write
c
c
c
      write(6,*)' k, zh, z, dz, dzh'
      do k=0,kmax+1
         write(6,100)k, zh(k), z(k), dz(k), dzh(k)
      enddo
      do k=0,kmax+1
         dzi(k) = 1./dz(k)
         dzhi(k)= 1./dzh(k)
      enddo
c
c
100   format(i5,4(x,f16.10))
c
      dzstart = dz(1)
c
      return
      end
