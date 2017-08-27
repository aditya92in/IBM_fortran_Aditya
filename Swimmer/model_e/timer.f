      subroutine timer(time)
c
      implicit none
c
      double precision time
      real   timereal
c
      call CPU_TIME(timereal)
c 
      time = timereal
c
C     time = time + 1.
      return
      end
