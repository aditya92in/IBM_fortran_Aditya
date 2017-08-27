      subroutine mea2(icase)
c
c  read average, adapt average, write average
c
      implicit none
c
      include 'param.txt'
      include 'common.txt'
c
      integer icase
c
      real Uaver(0:i1,0:j1)
      real Vaver(0:i1,0:j1)
      real Waver(0:i1,0:j1)
      real Uavsq(0:i1,0:j1)
      real Vavsq(0:i1,0:j1)
      real Wavsq(0:i1,0:j1)
      real UV   (0:i1,0:j1)
c
      save Uaver
      save Vaver
      save Waver
      save Uavsq
      save Vavsq
      save Wavsq
      save UV
c
      integer naver
      save    naver
      real    one_ze, const1, const2
      real    uwkp, uwkp4
      real    uvjp, uvjp4
      real    Etot, vol
      real    inbloc
      integer ip, jp, kp
c
c    
      integer iunit
c
      parameter(iunit=15)
c
      if(icase.eq.0)then
c
c initialise: read number of samples already taken, and old averages
c
      open(unit=iunit,file='data_time_averages.txt',form='unformatted')
      open(unit=27,file='data_Etot')
100     read(27,*,end=200)
        goto 100
200     continue
        backspace 27
c
c    
      naver = 0
 
      Uaver     = 0.
      Vaver     = 0.
      Waver     = 0.
 
      Uavsq     = 0.
      Vavsq     = 0.
      Wavsq     = 0.
      UV        = 0.
c
c     read(iunit)naver
C
c     read(iunit)Uaver
c     read(iunit)Vaver
c     read(iunit)Waver
c
c     read(iunit)Uavsq
c     read(iunit)Vavsq
c     read(iunit)Wavsq
c     read(iunit)UV
c
      close(iunit)
c
      write(6,*)'read means: naver  = ', naver
      endif
c
      if(icase.eq.1)then
c
c adapt time-average
c
      do j=1,jmax
      do i=1,imax
c
      ip = i + 1
      jp = j + 1
c
      Uaver(i,j) = (1./(naver+1))*( naver*Uaver(i,j) 
     &             + Unew(i,j) )
      Vaver(i,j) = (1./(naver+1))*( naver*Vaver(i,j) 
     &             + Vnew(i,j) )
c
      Uavsq(i,j) = (1./(naver+1))*( naver*Uavsq(i,j) 
     &             + Unew(i,j)**2 )
      Vavsq(i,j) = (1./(naver+1))*( naver*Vavsq(i,j) 
     &             + Vnew(i,j)**2 )

          uvjp  =  0.25*( Unew(i,jp)+Unew(i,j) )*
     &             ( Vnew(ip,j)+Vnew(i,j)  )
          uvjp4 = 0.25*(Unew(i,j)    +Unew(i,j+3))
     &                *(Vnew(i-1,j+1)+Vnew(i+2,j+1))
        one_ze = ordaru(i,j)

        const1 = one_ze*  9./8. + (1.-one_ze)*1.
        const2 = one_ze*(-1./8.)


      UV   (i,j) = (1./(naver+1))*( naver*UV   (i,j) 
     & + const1 * uvjp + const2 * uvjp4)
c
      enddo
      enddo
c
      naver  = naver + 1
c
c
c calculate total energy
c
      Etot = 0.
      vol  = 0.
      do j=1,jmax
      do i=1,imax

         Etot = Etot + 
     &   0.5*
     &   (
     &     Uavsq(i,j)    + Uavsq(i-1,j)   
     &   - Uaver(i,j)**2 - Uaver(i-1,j)**2
     &   + Vavsq(i,j)    + Vavsq(i,j-1)   
     &   - Vaver(i,j)**2 - Vaver(i,j-1)**2
     &   )*dx(i)*dy(j)
         vol = vol + inbloc * dx(i)*dy(j)
      enddo
      enddo

c      write(27,1000)time, Etot, Etot/vol
c1000  format(3(f16.8,1x))

c
      endif

c
c
      if(icase.eq.2)then
c
      open(unit=iunit,file='data_time_averages.txt',form='unformatted')
c
      write(iunit)naver
c
      write(iunit)Uaver
      write(iunit)Vaver
      write(iunit)Waver
c
      write(iunit)Uavsq
      write(iunit)Vavsq
      write(iunit)Wavsq
      write(iunit)UV
c
      close(iunit)
c
      close(27)
c
      endif
      return
      end
