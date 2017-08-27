      subroutine tec
c
c output in tecplot format
c
      implicit none
c
      include 'param.txt'
      include 'common.txt'
c
      integer iindex, jindex
      real uin, vin
       
c
c
      open(unit=10,file='tec2D_z.dat')
c
              write (10,*) 'title= "domain 2"'
        write (10,*) 'variables = "x","y"'
        write (10,*)             '"U","V","P","C"'
        write (10,*) 'zone t= "zone-title", f=point'
c
        write(10,*)'i = ', imax, ' , j = ', jmax 
c
        do j=1,jmax
        do i=1,imax
           uin = 0.5*(Unew(i,j) + Unew(i-1,j))
           vin = 0.5*(Vnew(i,j) + Vnew(i,j-1))
           write(10,2000)xp(i),yp(j),uin,vin,p(i,j),Cnew(i,j)
        enddo
        enddo
        close(10)
c
c
2000    format (6e14.6)
3000    format (6e14.6)
c
c gnuplot output
c
c
      open(unit=10,file='gnu1D_jvar')
c
      iindex = imax/2
      do j=1,jmax
         write(10,5000)j,yp(j),yh(j),unew(iindex,j),vnew(iindex,j),p(iindex,j)
      enddo
c
5000  format(i6,5(e16.8,1x))
      close(10)
c
c
      open(unit=10,file='gnu1D_ivar')
      jindex = jmax/2
      do i=1,imax
         write(10,1000)i,xp(i),xh(i),Unew(i,jindex),
     &                  Vnew(i,jindex),p(i,jindex)
      enddo
1000  format(i4, 5(e16.8,x))
      close(10)
c
c
      open(unit=10,file='gnu2D.dat')
c
        do j=1,jmax
        do i=1,imax
           uin = 0.5*(Unew(i,j) + Unew(i-1,j))
           vin = 0.5*(Vnew(i,j) + Vnew(i,j-1))
           write(10,2000)xp(i),yp(j),uin,vin,p(i,j),Cnew(i,j)
        enddo
        write(10,'(a)')
        enddo
        close(10)
c
c

        return
        end


