      subroutine tec_cropped(ib,ie,jb,je)
c
c output in tecplot format
c
      implicit none
c
      include 'param.txt'
      include 'common.txt'
      integer  ib,ie,jb,je
c
      integer iindex, jindex
      real uin, vin
       
c
c
      open(unit=10,file='tec2D_z_cropped.dat')
c
              write (10,*) 'title= "domain 2"'
        write (10,*) 'variables = "x","y"'
        write (10,*)             '"U","V","P"'
        write (10,*) 'zone t= "zone-title", f=point'
c
        write(10,*)'i = ', ie-ib+1, ' , j = ', je-jb+1 
c
        do j=jb,je
        do i=ib,ie
           uin = 0.5*(Unew(i,j) + Unew(i-1,j))
           vin = 0.5*(Vnew(i,j) + Vnew(i,j-1))
           write(10,2000)xp(i),yp(j),uin,vin,p(i,j)
        enddo
        enddo
        close(10)
c
c
2000    format (5e14.6)
c
        return
        end
