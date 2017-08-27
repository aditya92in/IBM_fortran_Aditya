      subroutine fkdat
c   
      implicit none
c
      include 'param.txt'
      include 'common.txt'
c
      real G05DAF,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, Yplus
      realCplus, kappa,Wvelo,ran
      parameter (kappa=0.4)
      parameter (Cplus = 5.3)
      integer inbloc
C DGF Benchmark
      real Umax, Height
      Umax = 1.5
      Height = yh(jmax)-yh(0)
c
c
      call G05CBF()
c
c *** Disturbances: order of magnitude friction velocity ( = 20 U*)
c
      tmp1 = - 1.5 ! 0. + 0.1
      tmp2 =   1.5 ! 0. + 0.1
      tmp3 = - 1.5 ! 0. + 0.1
      tmp4 =   1.5 ! 0. + 0.1
!     tmp1 = 0.
!     tmp2 = 0.
!     tmp3 = 0.
!     tmp4 = 0.
      ran  = 0.
 
        do j=1,jmax
          do i=1,imax
c
         Wvelo = 1.
         inbloc = 1
C        if(( (i.ge.ibl1-1) .and. (i.le. ibl2)).and.
C    &      ( (j.ge.jbl1  ) .and. (j.le. jbl2)).and.
C    &                                        )inbloc = 0
C     Wvelo = 1.*inbloc + Wvelo * 0.
C     Wvelo =  4.*Umax*yp(j)*(Height-yp(j))/Height**2 ! DFG
c
      Wvelo = 1.
!     ran = G05DAF( tmp5 , tmp6 )
      Unew(i,j) = Wvelo + ran
!!!!!! ATT         if(yp(j).lt.0.)Unew(0,j)     = 0.
!     ran = G05DAF( tmp5 , tmp6 )
      Uold(i,j) = Wvelo + ran
!!!!!! ATT          if(yp(j).lt.1.)Uold(0,j)     = 0.

!     ran = G05DAF( tmp3 , tmp4 )
      Vnew(i,j) = ran
!     ran = G05DAF( tmp3 , tmp4 )
      Vold(i,j) = ran
      Cold(i,j) = 0.
      Cnew(i,j) = 0.
c
      enddo
      enddo
c
c     unew(imax/2,jmax/2) = unew(imax/2,jmax/2)+0.02
c     uold(imax/2,jmax/2) = uold(imax/2,jmax/2)+0.02
c
      return
      end


