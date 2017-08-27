
      subroutine corgen(PL,dtin)
c
c general correc routine: ALL interior points inlcuded, PL
c is assumed to have the right BC
c
      implicit none
      include 'param.txt'
      include 'common.txt'
      real PL(1-ih:imax+ih,1-jh:jmax+jh)
      real dtin
      do j=1,jmax
      do i=1,imax
       dudt(i,j)= 
     1  dudt(i,j) - dtin*dxhi(i) * ( PL(i+1,j)-PL(i,j) )
      enddo
      enddo
      do j=1,jmax
      do i=1,imax
       dvdt(i,j)= 
     1  dvdt(i,j) - dtin*dyhi(j) * ( PL(i,j+1)-PL(i,j) )
      enddo
      enddo

         do j=1-jh,jmax+jh
            do i=1-ih,imax+ih
            Uold(i,j)=Unew(i,j)
            Unew(i,j)=dUdt(i,j)
            Vold(i,j)=Vnew(i,j)
            Vnew(i,j)=dVdt(i,j)
            enddo
        enddo
      return
      end
