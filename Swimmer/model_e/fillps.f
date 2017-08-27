
      subroutine fillps(pvar,dtin)
      implicit none
      include 'param.txt'
      include 'common.txt'
c     real pvar(imax, jmax, kmax)
      real pvar(1-ih:imax+ih,1-jh:jmax+jh)
      real   dti, dtin
      real   sum
      dti =1./dtin
         do j=1,jmax
            do i=1,imax
               pvar(i,j) =  dti*  (
     2  Rpi(j)*(Rv(j)*dVdt(i,j)-Rv(j-1)*dVdt(i,j-1) ) * dyi(j) +
     3  (dUdt(i,j)-dUdt(i-1,j) ) * dxi(i)     )
CVOL               pvar(i,j) =  dti*  (
CVOL     1  (dUdt(i,j)-dUdt(i-1,j) ) * dy(j) +
CVOL     2  (dVdt(i,j)-dVdt(i,j-1) ) * dx(i)  )
            enddo
         enddo
c 

      return
      end
