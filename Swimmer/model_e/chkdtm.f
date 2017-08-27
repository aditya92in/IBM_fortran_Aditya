      subroutine chkdtm
c
c also write indices grid cell with maximum limitation on time step
c and the velocity components here
c
      implicit none
      include 'param.txt'
      include 'common.txt'
      real    Courant,dtemp,dtmax, dtadv, dtdif, dtadvt, dtdift
      real dtmin
      integer itime
      data itime /0/
      save itime
      integer iindex,jindex
      real    umax, vmax
      Courant  = .1
      Courant = Courant !+ 0.05*real(istep)/real(nstep)
      dtmax   = 1.e-14
      dtadv   = 1.e-14
      dtdif   = 1.e-14
c
      iindex = 0
      jindex = 0
      Umax = 999.999
      Vmax = 999.999
c
         do j=1,jmax
            do i=1,imax
            dtemp =
     1      abs(Unew(i,j) * dxi(i)) +
     2      abs(Vnew(i,j) * dyi(j)) +
C ATT
c    4      2. * visc(i,j) * (dxi*dxi + dyi*dyi )
     4   2.*visc(i,j) * (dxi(i)*dxi(i) + dyi(j)*dyi(j) )



            dtadvt =      
     1      abs(Unew(i,j) * dxi(i)) +
     2      abs(Vnew(i,j) * dyi(j)) 
            dtdift =
     4   visc(i,j) * (dxi(i)*dxi(i) + dyi(j)*dyi(j) )
            if(dtadvt.gt.dtadv)then
               iindex = i
               jindex = j
               umax = Unew(i,j)
               vmax = Vnew(i,j)
            endif
            dtadv = max(dtadv, dtadvt)
            dtdif = max(dtdif, dtdift)
            dtmax = max ( dtmax , dtemp )
            enddo
         enddo
      dt = Courant / dtadv
      write(6,*)'dt adv = ',dt
      dt = Courant / dtdif
      write(6,*)'dt dif = ',dt
      write(6,1000)'indeces etc:',
     & iindex, jindex,  Umax, Vmax
1000   format(a,2i6,2(f16.8,1x))
      
      dt = Courant / dtmax
C     dtmin = 0.0004 + 0.0002 *(real(itime)/nstep)
C     dtmin = 0.0001
C     dt = min(dtmin, Courant / dtmax)
      itime = itime + 1
      write(6,*)'dt = ',dt
      return
      end
