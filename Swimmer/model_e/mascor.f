      subroutine mascor(U,V,W)
c
c correct overall mass inflow-outflow
c
      include 'param.txt'
      include 'block.inc'
      include 'common.txt'
      integer  iafcor
      real U(0:i1,0:j1,0:k1),V(0:i1,0:j1,0:k1),
     $     W(0:i1,0:j1,0:k1)
c
      real summas, sumU
      real in, out
c
      summas = 0.
c
      do i=1,imax
        do j=1,jmax
c          summas = summas +( -W(i,j,0 )  + W(i,j,imax))*dx(i)*dy(j)
           summas = summas +( -W(i,j,0 )  + W(i,j,kmax))*dx(i)*dy(j)
        enddo
      enddo
c
       do k=1,kmax
         do i=1,imax
           summas = summas +( -V(i,0,k ) + V(i,jmax,k))*dx(i)*dz(k)
         enddo
       enddo
c
C      write(6,*)'summas inbetween = ', summas
       sumU = 0.
       in   = 0.
       out  = 0.
c
      in  = 0.
      out = 0.
      do k=1,kmax
        do j=1,jmax
           summas = summas +( -U(0 ,j,k) + U(imax,j,k))*dy(j)*dz(k)
           sumU   = sumU   + U(imax,j,k)*dy(j)*dz(k)
            out   =  out   + U(imax,j,k)*dy(j)*dz(k)
             in   =   in   + U(0,j,k)*dy(j)*dz(k)
         enddo
       enddo
c
C      write(6,*)'summas, in, out = ', summas, in, out
c
c correct mass flow at outflow boundary so that total mass inflow-outflow = 0.
c
c summas is too big by an amount summas. By multiplying U, V, W at the
c outflow boundary by a correction factor corfac, it will be exactly 0.
c
c  SUM (1.-corfac)*U*dy*dz       = summas -->
c
c       1. - corfac              = summas / (SUM U*dy*dz)
c
c            corfac              = 1 - summas / (SUM U*dy*dz)
c
      corfac = 1. - summas / sumU
c
c     write(6,*)' summas, sumU, corfac', summas, sumU, corfac
c
c perform correction
c
      do k=0,kmax+1
        do j=0,jmax+1
           U(imax,j,k) = U(imax,j,k) * corfac
c          V(i1,j,k)   = V(i1,j,k)   * corfac
c          W(i1,j,k)   = W(i1,j,k)   * corfac
         enddo
       enddo
c      
c correct mass flow at outflow boundary so that total mass inflow-outflow = 0.
c      
      return
      end
