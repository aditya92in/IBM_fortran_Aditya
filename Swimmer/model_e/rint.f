      real function rint(x1,x2,r)
      implicit none
      real x1, x2, r
      real upper, lower
      if(r**2 .le. x2**2)stop 'KUT2'
      if(r**2 .le. x1**2)stop 'KUT1'
      upper = 0.5*x2*sqrt(r**2-x2**2) + 0.5*r**2*asin(x2/r)
      lower = 0.5*x1*sqrt(r**2-x1**2) + 0.5*r**2*asin(x1/r)
      rint = upper - lower
      return
      end
