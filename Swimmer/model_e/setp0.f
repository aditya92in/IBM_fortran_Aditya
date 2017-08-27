      subroutine setp0(pvar)
      include 'param.txt'
      include 'common.txt'
      real pvar(1-ih:imax+ih,1-jh:jmax+jh)
      real p0
c
      p0 = p(1,jmax/2)
c
      do j=1,jmax
      do i=1,imax
         pvar(i,j) = pvar(i,j) - p0
      enddo
      enddo
c
       return
       end
