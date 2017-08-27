      subroutine adamsles(string)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real   dold(1-ih:imax+ih,1-jh:jmax+jh)
      real   dnew(1-ih:imax+ih,1-jh:jmax+jh)
      real   pr  (1-ih:imax+ih,1-jh:jmax+jh)
      save   pr
      real   cof1,cof2
      character*5 string
*     Adams-Bashfort  cof1=1.5  , cof2 = -0.5
*     Euler-Forward   cof1=  1  , cof2 =    0
      integer ifirst
      data ifirst /0/
      save ifirst
      if (string .eq. 'euler') then
      cof1=1.
      cof2=0.
      endif
      if (string .eq. 'adams') then
      cof1= 1.5
      cof2=-0.5
      endif

      call momxles(dold,uold,vold,
     &          1,imax,1,jmax,visc,pold)
      call momxles(dnew,unew,vnew,
     &          1,imax,1,jmax,visc,p)
        do j=1,jmax
          do i=1,imax
          dudt(i,j)=Unew(i,j) +
     1    dt * ( cof1 * dnew(i,j) + cof2 * dold(i,j)  )
     1    - dt * dxhi(i) * ( p(i+1,j)-p(i,j) )
          enddo
        enddo

      call momyles(dold,uold,vold,
     &          1,imax,1,jmax,visc,pold)
      call momyles(dnew,unew,vnew,
     &          1,imax,1,jmax,visc,p)
        do j=1,jmax
          do i=1,imax
          dvdt(i,j)=Vnew(i,j)+
     1    dt * ( cof1 * dnew(i,j) + cof2 * dold(i,j) )
     1    - dt*dyhi(j) * ( p(i,j+1)-p(i,j) )
          enddo
        enddo
      return
      end
