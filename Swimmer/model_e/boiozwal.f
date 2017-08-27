      subroutine boiozwal(U,V,W,iafcor)
c
c wall at z-planes, io in x-direction, wall in y-direction
c
      include 'param.txt'
      include 'block.inc'
      include 'common.txt'
      integer  iafcor
      real U(1-ih:imax+ih,1-jh:jmax+jh,1-kh:kmax+kh)
      real V(1-ih:imax+ih,1-jh:jmax+jh,1-kh:kmax+kh)
      real W(1-ih:imax+ih,1-jh:jmax+jh,1-kh:kmax+kh)
c
      real Umax, Height
      Umax = 1.5
      Height = yh(jmax)-yh(0)
c
      if(iafcor.eq.0)then
      do k=0,k1
        do j=0,j1
          U(0 ,j,k)=1.!4.*Umax*yp(j)*(Height-yp(j))/Height**2
          V(0 ,j,k)=-V(1,j,k)*0.
          W(0 ,j,k)=-W(1,j,k)*0.
          U(i1,j,k)=U(imax,j,k)
          V(i1,j,k)=V(imax,j,k)
          W(i1,j,k)=W(imax,j,k)
C         U(imax,j,k)=U(imax-1,j,k)
C         V(i1,j,k)=V(imax,j,k)
C         W(i1,j,k)=W(imax,j,k)
         enddo
       enddo
c
       elseif(iafcor.eq.1)then
c
c U at end of domain NOT reset!!
c
c i-direction:
c
      do k=0,k1
        do j=0,j1
          U(0 ,j,k)= 1.!4.*Umax*yp(j)*(Height-yp(j))/Height**2
          V(0 ,j,k)=-V(1,j,k)*0.
          W(0 ,j,k)=-W(1,j,k)*0.
C         U(i1,j,k)=U(imax,j,k)
C         V(i1,j,k)=V(imax,j,k)
C         W(i1,j,k)=W(imax,j,k)
c         U(imax,j,k)=U(imax-1,j,k)
          V(i1,j,k)=V(imax,j,k)
          W(i1,j,k)=W(imax,j,k)
         enddo
       enddo
       endif
c
c j-direction:
c
       do k=0,k1
         do i=0,i1
           U(i,0,k )= U(i,1,k)
           V(i,0,k )= 0.
           W(i,0,k )= W(i,1,k)

           U(i,j1,k)= U(i,jmax   ,k)
           V(i,j1,k)= V(i,jmax   ,k)
           V(i,jmax,k)=0.
           W(i,j1,k)= W(i,jmax   ,k)
         enddo
       enddo
c
c k-direction:
c
      do i=1-ih,imax+ih
        do j=1-jh,jmax+jh
         W(i,j,0    )   =      0. 
         W(i,j,kmax )   =      0.
         V(i,j,0    )   =   V(i,j,1   )
         V(i,j,k1   )   =   V(i,j,kmax)
         U(i,j,0    )   =   U(i,j,1   )
         U(i,j,k1   )   =   U(i,j,kmax)
        enddo
      enddo
c
      return
      end
