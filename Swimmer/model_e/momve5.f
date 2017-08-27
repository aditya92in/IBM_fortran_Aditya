      subroutine momx4(du,u,v,
     &                ib,ie,jb,je,
     $                vi, iorder,pvar)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer ib,ie,jb,je
      integer ip,im,jp,jm
      integer iorder
      real du(1-ih:imax+ih,1-jh:jmax+jh),
     $     u   (1-ih:imax+ih,1-jh:jmax+jh),
     $     v   (1-ih:imax+ih,1-jh:jmax+jh),
     $     pvar(1-ih:imax+ih,1-jh:jmax+jh)
      real vi
c
c
      real uuip ,uuim ,uvjp ,uvjm ,uwkp ,uwkm
      real du1ip,du1im,du1jp,du1jm,du1kp,du1km
c
      real pip , pim
c
c
      
        do j=jb,je
          do i=ib,ie
*
          ip = i + 1
          jp = j + 1
          im = i - 1
          jm = j - 1

          uuip  = 0.25 * ( U(ip,j)+U(i,j) )*( U(ip,j)+U(i,j)  )
          uuim  = 0.25 * ( U(im,j)+U(i,j) )*( U(im,j)+U(i,j)  )
          uvjp  =  0.25*( U(i,jp)+U(i,j) )*
     &             ( V(ip,j)+V(i,j)  )
          uvjm  =  0.25*( U(i,jm)+U(i,j) )*
     &             ( V(ip,jm)+V(i,jm))


          du1ip = dxi(ip)  * (U(ip,j)-U(i ,j) )
          du1im = dxi(i)   * (U(i ,j)-U(im,j) )
          du1jp = dyhi(j)  * (U(i,jp)-U(i, j) )
          du1jm = dyhi(jm) * (U(i,j )-U(i,jm) )

c
          pip   = 0.!pvar(i+1,j)
          pim   = 0.!pvar(i  ,j)
c

C  ATT         du(i,j) =2.+
C        du(i,j) =8./50.+
C        du(i,j) =8./100. + 
         du(i,j) =0.00 + 
     1  dxhi(i)*( -uuip + vi*du1ip+uuim - vi*du1im -pip + pim )+
     2  dyi(j)*Rpi(j)*(  Rv(j)*(-uvjp + vi*du1jp)-Rv(j-1)*(-uvjm + vi*du1jm) )
            enddo
         enddo
      return
      end

      subroutinemomy4(dv,u,v,
     &                ib,ie,jb,je,
     $                vi, iorder,pvar)

      implicit none
      include 'param.txt'
      include 'common.txt'
      integer ib,ie,jb,je
      integer im,ip,jm,jp
      integer iorder
      real dv(1-ih:imax+ih,1-jh:jmax+jh),
     $     u   (1-ih:imax+ih,1-jh:jmax+jh),
     $     v   (1-ih:imax+ih,1-jh:jmax+jh),
     $     w   (1-ih:imax+ih,1-jh:jmax+jh),
     $     pvar(1-ih:imax+ih,1-jh:jmax+jh)
      real vi
c
      real uvip ,uvim ,vvjp ,vvjm ,wvkp ,wvkm
      real dv1ip,dv1im,dv1jp,dv1jm,dv1kp,dv1km
c
      real pjp , pjm
c
        do j=jb,je
          do i=ib,ie
*
          ip = i + 1
          jp = j + 1
          im = i - 1
          jm = j - 1
          uvip  =  (U(i ,j)+U(i ,jp))*
     &                   0.25*( V(i,j)+V(ip,j) )
          uvim  =  (U(im,j)+U(im,jp))*
     &                   0.25*( V(i,j)+V(im,j) )
          vvjp  = 0.25 * (V(i,j )+V(i,jp) )*
     &                   ( V(i,j)+V(i,jp) )
          vvjm  = 0.25 * (V(i,j )+V(i,jm) )*
     &                   ( V(i,j)+V(i,jm) )

c
c
          dv1ip = dxhi(i)  * ( V(ip,j) - V(i ,j) )
          dv1im = dxhi(im) * ( V(i ,j) - V(im,j) )
          dv1jp = dyi(jp)  * ( V(i,jp) - V(i,j ) )
          dv1jm = dyi(j)   * ( V(i,j ) - V(i,jm) )

C
          pjp   = 0.!pvar(i,j+1)
          pjm   = 0.!pvar(i,j  )
c
          dv(i,j)= 
     1 dxi(i)*((-uvip + vi    *dv1ip)-(-uvim + vi     *dv1im))+
     2 dyhi(j)*Rvi(j)*((Rp(j+1)*(-vvjp +vi    *dv1jp))-(Rp(j)*(-vvjm+vi       *dv1jm))) + dyhi(j)*(-pjp + pjm)
     3  -v(i,j)*Rvi(j)**2
                   
          enddo
        enddo
      return
      end
