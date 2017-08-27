      subroutine basgr(z, zh, dz, dzh, dzi, dzhi, k1)
c
c
      implicit none
      include 'mpif.h'
c
      integer k, k1, kmax
c
      real z(0:k1)
      real zh(0:k1)
      real dz(0:k1)
      real dzh(0:k1)
      real dzi(0:k1)
      real dzhi(0:k1)
c
c help variables
c
      real pi, expk, zhelp, aconst, atanh, ztmp, power
      parameter (power = 1.1)
c
      integer myid, ierr
c
      pi = 4.*atan(1.)
c
      kmax = k1 - 1
c
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
c
c first distribution zh
c
      do k=0, kmax/2
        zh(k) = (1.*k)**power
      enddo
      ztmp = zh(kmax/2)
      do k=0, kmax/2
        zh(k) = 0.5*zh(k)/ztmp
      enddo
      do k= kmax/2+1, kmax
         zh(k) = 1.-zh(kmax-k)
      enddo
c
c
      zh(k1) = zh(kmax) + (zh(kmax) - zh(kmax-1))
c
c
      do k=1,k1
         z(k) = zh(k-1) + 0.5*(zh(k) - zh(k-1))
      enddo
c
      z(0) =  z(1) -2.*(z(1) - zh(0))
c
c
      do k=1,k1
         dz(k) = zh(k) - zh(k-1)
      enddo
c
c
      dz(0) = dz(1)
c
c
      do k=0,kmax
         dzh(k) = z(k+1) - z(k)
      enddo
c
c
      dzh(k1) = dzh(kmax)
c
c
c write
c
c
      if(myid.eq.0)then
c
      write(6,*)' k, zh, z, dz, dzh'
      do k=0,k1
         write(6,100)k, zh(k), z(k), dz(k), dzh(k)
      enddo
      do k=0,k1
         dzi(k) = 1./dz(k)
         dzhi(k)= 1./dzh(k)
      enddo
c
      endif
c
100   format(i5,4(x,f16.10))
c
      return
      end
