      subroutine boucha(U,V,iperiox,iperioy)
      include 'param.txt'
      include 'common.txt'
      real U(1-ih:imax+ih,1-jh:jmax+jh)
      real V(1-ih:imax+ih,1-jh:jmax+jh)
      integer iperiox,iperioy

c
      do i=1,ih
        do j=1-jh,jmax+jh
          U(1-i ,j)  =U(imax+1-i,j)
          V(1-i ,j)  =V(imax+1-i,j)
          U(imax+i,j)=U(i,j)
          V(imax+i,j)=V(i,j)
         enddo
       enddo
c
c
      if(iperioy.eq.0)then
      do j=1,jh
         do i=1-ih,imax+ih
Cperio           U(i,1-j )  =U(i,jmax+1-j)
Cperio           V(i,1-j )  =V(i,jmax+1-j)
Cperio           U(i,jmax+j)=U(i,j   )
Cperio           V(i,jmax+j)=V(i,j   )
           U(i,1-j )  =-U(i,j)
           V(i,1-j )  =0.
           U(i,jmax+j)=-U(i,jmax+1-j )
           V(i,jmax+j-1)=0.
         enddo
       enddo
       else
      do j=1,jh
         do i=1-ih,imax+ih
           U(i,1-j )  =U(i,jmax+1-j)
           V(i,1-j )  =V(i,jmax+1-j)
           U(i,jmax+j)=U(i,j   )
           V(i,jmax+j)=V(i,j   )
         enddo
       enddo
       endif
c
       return
       end
