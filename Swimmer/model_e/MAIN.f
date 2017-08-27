c  
C    3D channel model for the direct numerical simulation
C    of turbulent flows.
C   
C    Bendiks Jan Boersma
C    adapted mathieu pourquie
C
C 
C     parameter list:
C     Uold(..),Vold(...),Wold(...) Velocities at the oldest
C                                  timestep
C  
c version: in-out BC, 1 pressure-correction per time step
c
 




      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'timing.inc'
      real     deltap(1-ih:imax+ih,1-jh:jmax+jh)
      real     psum
C     real     pr   (1-ih:imax+ih,1-jh:jmax+jh)
      real     pnorm, cof
      real    lx,ly,lz,Re
      integer iter, maxit
      real    parlen, zan
!      integer istart, ipt
      integer  ipt
      integer iperiox,iperioy
      
      real    zstart, rlength, dzstart
      read(5,*)istart
      read(5,*) Ufree
      read(5,*) cs  
      read(5,*) lx
      read(5,*) ly
      read(5,*) lz
      read(5,*) Re
      read(5,*) nstep
      iperiox = 0
      iperioy = 0
********
      pold = 0.
      unew = 0.
      vnew = 0.
      uold = 0.
      vold = 0.
      dudt = 0.
      dvdt = 0.

c     call basgrro(lx,ly,Re)
      call reagr2(lx,ly)
      zstart  = 0.2
      rlength = 100.
      Rei = 1.   / Re
      unew = 0.
      vnew = 0.
      uold = 0.
      vold = 0.
      dudt = 0.
      dvdt = 0.
      if(istart.eq.0)then
      call fkdat
      p      = 0.
      pold   = 0.
C     call makepr(pr,1)
      else
c     call loadd(0)
C     call loadd_peskin(part,part0,upart,imax,jmax,ih,jh,npts,0)
c      open(unit=10,file='time')
c      read(10,*)time
c      close(10)

C     call makepr(pr,0)
      endif
C      pold = p
C      call addist

      cof = 1.
      call bouio(Uold,Vold)
      call bouio(Unew,Vnew)
      call bouioc(Cnew,Unew,Vnew)
      call bouioc(Cold,Uold,Vold)
c
c overwrite for wall law
c
c
         deltap = 0.
c
c
c initialise/update  means
c See the routine to check it does what you want!
c
      write(6,*)'HIER'
c     call mea2(0)
c
c
c
c initialise force function
c
      open(unit=27,file='liftdrag',access='append')
      write(6,*)'HIER2'
      open(unit=33,file='conv')
      visc = Rei
      write(6,*)'Rei, visc = ', Rei
      call timer(CPUprog)
      CPUforce = 0.
c
c MAIN loop starts here
c
c
c perturb velocity, remove after adding aobstacle
c
       unew(imax/2,jmax/2) = 0.

      do istep = 1,100
      call chkdtm
      time = time + dt
      write(6,*) 'time, dt = ', time,dt
      call adamsb('adams')
      
      iter = iter + 1
      call bouio   (dUdt,dVdt)
      call bouioc(dCdt,dUdt,dVdt)
      call fillps   (deltap,dt)

5000  format(3(i4),f16.10)
      
      call poisson(deltap,imax,jmax,ih,jh,
     &                      dxi,dyi,dxhi,dyhi,Rp,Rv,Rpi,Rvi,iperiox,iperioy)
         call setp0(deltap)
         call boupio(deltap)
         if((istart.eq.0).and.(istep.eq.1))then
         pold = deltap
         else
         pold = p
         endif
         p = p + deltap
         call setp0(p)
         call setp0(pold)
         call boupio(p)
         call boupio(pold)
CFORC        call boupob(p)
      call corgen(deltap,dt)
      call stepscal
     
c
c     call copy
      call bouio   (Unew,Vnew)
      call bouioc(Cnew,Unew,Vnew)
c
      call chkdiv ! WAS call chkdiv ! WAS chkblo
c
c update means
         if((istep/10)*10 .eq. istep)then
          call mea2(1)
          end if
c
       end do  ! END MAIN loop
c      call timer(t1)
c      CPUprog = t1 - CPUprog
c      write(6,*)'CPU program = ',CPUprog
c      close(33)
c
c       call loadd(1)
c
c write means
c
              call mea2(2)
      
              call tec
             open(unit=10,file='time')
             write(10,*)time
             close(10)
c      close(17)

c      close(31)
      
      stop
      end
