      subroutine bouio(U,V)
      include 'param.txt'
      include 'common.txt'
      real U(1-ih:imax+ih,1-jh:jmax+jh)
      real V(1-ih:imax+ih,1-jh:jmax+jh)

c
      do i=1,ih
        do j=1-jh,jmax+jh
          U(0,j)     = 1.
!!!!! ATT          if(yp(j).lt.0.)U(0,j)     = 0.
          V(0,j)     = 0.
          U(imax+1,j)=U(imax,j)
          V(imax+1,j)=V(imax,j)
         enddo
       enddo
c
c
      do j=1,jh
         do i=1-ih,imax+ih
Cperio           U(i,1-j )  =U(i,jmax+1-j)
Cperio           V(i,1-j )  =V(i,jmax+1-j)
Cperio           U(i,jmax+j)=U(i,j   )
Cperio           V(i,jmax+j)=V(i,j   )
           U(i,1-j )  =  U(i,j)
           V(i,1-j )  =0.
           U(i,jmax+j)= -U(i,jmax+1-j )
           V(i,jmax+j-1)=0.
         enddo
       enddo
c
       return
       end
