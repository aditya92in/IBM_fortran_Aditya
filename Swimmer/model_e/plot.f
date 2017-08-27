      program plot
      npt = 100
      dx = 1./8
      r1d4dx=1./(4.*dx) 
      r1d2dx=1./(2.*dx) 
      xt = 0.5
      pi = 4. * atan(1.)
      sum = 0.
      open(unit=10,file='pl')
      do i=1,npt
      x = (i-1.)/(npt-1.)
      y = r1d4dx*(1.+cos(pi*(xt-x)*r1d2dx))
      if(abs((xt-x)).gt.2.*dx)y=0.
      write(10,*)x,y
      sum = sum + y/(npt-1.)
      enddo
      close(10)
      write(6,*)'integral = ', sum
      stop
      end
