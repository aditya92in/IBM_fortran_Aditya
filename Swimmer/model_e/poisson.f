c solbox: all solid walls
c
cMathieu added variables wi,wj,xrt,yrt,a,b,c,bb,d,vfti,vftj
c in header to prevent use of stack variables
c 28 Oct 1996
      subroutine poisson(p,imax,jmax,ih,jh,
     &                      dxi,dyi,dxhi,dyhi,Rp,Rv,Rpi,Rvi,iperiox,iperioy)
      implicit none
      integer ih,jh
      integer imax,jmax
      real work(2*imax*jmax)
      real dxi(1-ih:imax+ih), dyi (1-jh:jmax+jh)
      real dxhi(1-ih:imax+ih),dyhi(1-jh:jmax+jh)
      real Rp(1-jh:jmax+jh),Rv(1-jh:jmax+jh),Rpi(1-jh:jmax+jh),Rvi(1-jh:jmax+jh)
      integer iperiox,iperioy
      integer i,j
      integer ier
      real     p(1-ih:imax+ih,1-jh:jmax+jh)
c  LOCAL VARIABLES
      real    pi,z,d(imax,jmax)
      real bin(imax)
      real ax(imax),bx(imax),cx(imax)
      real ay(jmax),by(jmax),cy(jmax)
      real  y(imax,jmax)

      pi =    4.*atan(1.)
C****************************
      doi=1,imax
        ax(i) = 1.*dxi(i)*dxhi(i-1)
        cx(i) = 1.*dxi(i)*dxhi(i)
        bx(i) = -ax(i) - cx(i)
      enddo
c
c BC: p given in virtual pint OUTSIDE domain
c
c Neumann
c
      if(iperiox.eq.0)then
      bx(1)    = bx(1) + ax(1)
      ax(1)    = 0.
      else
c periodic
c     bx(1)    = bx(1)
c     ax(1)    = 0.
      endif
c
c Dirichlet
c
      if(iperiox.eq.0)then
      bx(imax) = bx(imax) - cx(imax)
      cx(imax) = 0.
      else
c periodic
      bx(imax) = bx(imax) 
c     cx(imax) = 0.
      endif
c
c Neumann
c
c     bx(imax) = bx(imax) + cx(imax)
c     cx(imax) = 0.
c
c fill coefficients in j-direction
c
      do j=1,jmax
         ay(j) = Rpi(j)*Rv(j-1)*dyi(j)*dyhi(j-1)
         cy(j) = Rpi(j)*Rv(j)*dyi(j)*dyhi(j)
         by(j) = -ay(j) - cy(j)
      enddo
c
c BC: Neumann
c
      if(iperioy.eq.0)then
      by(1)    = by(1) + ay(1)
      ay(1)    = 0.
      by(jmax) = by(jmax) + cy(jmax)
      cy(jmax) = 0.
      endif
c periodic
c     by(1)    = by(1) 
c     ay(1)    = 0.
c     by(jmax) = by(jmax) 
c     cy(jmax) = 0.
c
c
      do i=1,imax
         bin(i) = bx(i) 
      enddo
c
c
      do j=1,jmax
      do i=1,imax
         y(i,j) = P(i,j) 
      enddo
      enddo
      CALL BLKTRI(0,1-iperioy,jmax,ay,by,cy,1-iperiox,imax,ax,bin,cx,imax,y
c
c                   ^ 0 for perioic BC
c
     $     ,ier,work)
c     write(6,*)'ier = ', ier
c
      CALL BLKTRI(1,1-iperioy,jmax,ay,by,cy,1-iperiox,imax,ax,bin,cx,imax,y
     $     ,ier,work)
c
c     write(6,*)'ier = ', ier
c
      do j=1,jmax
      do i=1,imax
         P(i,j) = y(i,j) 
      enddo
      enddo
c
      return
      end
