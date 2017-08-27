
      subroutine stepscal
c
c
      implicit none
      include 'param.txt'
      include 'common.txt'

         do j=1-jh,jmax+jh
            do i=1-ih,imax+ih
            Cold(i,j)=Cnew(i,j)
            Cnew(i,j)=dCdt(i,j)
            enddo
        enddo
      return
      end
