      subroutine bouioc(C,U,V)
      include 'param.txt'
      include 'common.txt'
      real C(1-ih:imax+ih,1-jh:jmax+jh)
      real U(1-ih:imax+ih,1-jh:jmax+jh)
      real V(1-ih:imax+ih,1-jh:jmax+jh)

c
      do i=1,ih
        do j=1-jh,jmax+jh
          C(0,j)     = 1.
          if(yp(j).lt.0.)C(0,j)     = 0.
          C(imax+1,j)=C(imax,j)
         enddo
       enddo
c
c
      do j=1,jh
         do i=1-ih,imax+ih
Cperio           C(i,1-j )  =C(i,jmax+1-j)
Cperio           C(i,jmax+j)=C(i,j   )
           C(i,1-j )  = C(i,j)
           C(i,jmax+j)= C(i,jmax+1-j )
         enddo
       enddo
c
       return
       end
