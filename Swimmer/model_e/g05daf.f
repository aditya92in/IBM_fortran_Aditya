      double precision function g05daf(t1,t2)
c
c uses f90 random number generator
c
      implicit none
c
      double precision t1, t2
      real ran
c
c     call srandom(1)
c
      call random_number(ran )
c
      g05daf = t1 + (t2-t1)*ran
c
c
      return
      end

