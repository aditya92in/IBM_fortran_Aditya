      subroutine boupio(pvar)
      include 'param.txt'
      include 'common.txt'
      real pvar(1-ih:imax+ih,1-jh:jmax+jh)
      do i=1,ih
        do j=1-jh,jmax+jh
          pvar(1-i ,j)=pvar(i,j)
          pvar(imax+i,j)=-pvar(imax+1-i,j)
         enddo
       enddo
      do j=1,jh
         do i=1-ih,imax+ih
Cperio           pvar(i,1-j )=pvar(i,jmax+1-j)
Cperio           pvar(i,jmax+j)=pvar(i,j   )
           pvar(i,1-j )=pvar(i,j)
           pvar(i,jmax+j)=pvar(i,jmax+1-j   )
         enddo
       enddo
c
c
       return
       end
