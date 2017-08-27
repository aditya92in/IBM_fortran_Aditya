      program ruk
      double precision g05daf
      external g05cbf,g05daf
      integer i
      call g05cbf(i)
      x =  g05daf(-1.d0,1.d0)
      write(6,*)x
      x =  g05daf(-1.d0,1.d0)
      write(6,*)x
      x =  g05daf(-1.d0,1.d0)
      write(6,*)x
      x =  g05daf(-1.d0,1.d0)
      write(6,*)x
      x =  g05daf(-1.d0,1.d0)
      write(6,*)x
      stop
      end
