

      subroutine chkdiv
      implicit none
      include 'param.txt'
      include 'common.txt'
      real div   ,divmax,divtot
      real divflu,divmaf,divtof
      integer iindex, jindex
      divmax = -9999.9999
      divtot = 0.
      divmaf = -9999.9999
      divtof = 0.
      iindex = -20
      jindex =-20
         do j=1,jmax
            do i=1,imax
            div = (
     2   Rpi(j)*( Rv(j)*Vnew(i,j)-Rv(j-1)*Vnew(i,j-1) ) * dyi(j) +
     3   ( Unew(i,j)-Unew(i-1,j) ) * dxi(i) )
            div = div*Rp(j)*dy(j)*dx(i)
            divflu = (
     2   ( Rv(j)*Vnew(i,j)-Rv(j-1)*Vnew(i,j-1) ) * dx(i)   +
     3   Rp(j)*( Unew(i,j)-Unew(i-1,j) ) * dy(j)   )
c
c
           divtot = divtot + div*dx(i)*Rp(j)*dy(j)
           divmax = max ( abs(div),divmax)
           divtof = divtof + divflu
           if (abs(divflu) .gt. divmaf)then
              iindex = i
              jindex = j
           endif
           divmaf = max ( abs(divflu),divmaf)
!          write(6,*)'TS ',i,j,abs(divflu),abs(div)*Rp(j)*dx(i)*dy(j)
!          write(6,*)'TS ',i,j,Rpi(j)
           enddo
         enddo
      write(6,111) divtot,divmax
111   FORMAT('Mass loss/gain : Tot =',e13.6,' Max = ',e13.6)
      write(6,222) divtof,divmaf
222   FORMAT('Flux mass loss/gain : Tot =',e13.6,' Max = ',e13.6)
c     if (divtot .gt. 10e-4) stop 'Fatal Error **DIVERGENCE TO LARGE**'
      write(6,*)'iindex, jindex = ', iindex, jindex
      return
      end
