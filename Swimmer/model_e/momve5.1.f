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
      real uuip4 ,uuim4 ,uvjp4 ,uvjm4 ,uwkp4 ,uwkm4
      real du1ip4,du1im4,du1jp4,du1jm4,du1kp4,du1km4
      real const1, const2, one_ze
      real pip , pim
      real pip4, pim4
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

          uuip4 = 0.25*(U(i,j)    +U(i+3,j))**2
          uuim4 = 0.25*(U(i-3,j)  +U(i,j)  )**2
          uvjp4 = 0.25*(U(i,j)    +U(i,j+3))
     &                *(V(i-1,j+1)+V(i+2,j+1))
          uvjm4 = 0.25*(U(i,j-3)  +U(i,j)  )
     &                *(V(i-1,j-2)+V(i+2,j-2))

          du1ip = dxi(ip)  * (U(ip,j)-U(i ,j) )
          du1im = dxi(i)   * (U(i ,j)-U(im,j) )
          du1jp = dyhi(j)  * (U(i,jp)-U(i, j) )
          du1jm = dyhi(jm) * (U(i,j )-U(i,jm) )

          du1ip4= dxhi4E(i)   * (U(i+3,j)-U(i ,j) )
          du1im4= dxhi4W(i)   * (U(i ,j) -U(i-3,j) )
          du1jp4= dyi4N (j)   * (U(i,j+3)-U(i, j) )
          du1jm4= dyi4S (j)   * (U(i,j ) -U(i,j-3) )
c
          pip   = pvar(i+1,j)
          pim   = pvar(i  ,j)
          pip4  = pvar(i+2,j)
          pim4  = pvar(i-1,j)
c
        one_ze = ordaru(i,j)
c
        const1 = one_ze*  9./8. + (1.-one_ze)*1.
        const2 = one_ze*(-1./8.) 
c

C  ATT         du(i,j) =2.+
C        du(i,j) =8./50.+
C        du(i,j) =8./100. + 
         du(i,j) =0.00 + 
     1  const1 *(
     1  dxhi(i)*( -uuip + vi*du1ip+uuim - vi*du1im -pip + pim )+
     2  dyi(j)*(  -uvjp + vi*du1jp+uvjm - vi*du1jm )
     7          )
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
      real const1, const2, one_ze
c
      real uvip4 ,uvim4 ,vvjp4 ,vvjm4 ,wvkp4 ,wvkm4
      real dv1ip4,dv1im4,dv1jp4,dv1jm4,dv1kp4,dv1km4
      real pjp , pjm
      real pjp4, pjm4
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
          uvip4 =  (U(i+1 ,j-1)+U(i+1 ,j+2))*
     &                   0.25*( V(i,j)+V(i+3,j) )
          uvim4 =  (U(i-2,j-1)+U(i-2,j+2))*
     &                   0.25*( V(i,j)+V(i-3,j) )
          vvjp4 =  0.25*(V(i,j )  +V(i,j+3) )**2
          vvjm4 =  0.25*(V(i,j-3 )+V(i,j)   )**2

          dv1ip = dxhi(i)  * ( V(ip,j) - V(i ,j) )
          dv1im = dxhi(im) * ( V(i ,j) - V(im,j) )
          dv1jp = dyi(jp)  * ( V(i,jp) - V(i,j ) )
          dv1jm = dyi(j)   * ( V(i,j ) - V(i,jm) )

C
C
          dv1ip4= dxi4E (i)   * ( V(i+3,j  )  - V(i  ,j) )
          dv1im4= dxi4W (i)   * ( V(i  ,j  )  - V(i-3,j) )
          dv1jp4= dyhi4N(j)   * ( V(i,j+3  )  - V(i  ,j  ) )
          dv1jm4= dyhi4S(j)   * ( V(i,  j  )  - V(i  ,j-3) )
c
          pjp   = pvar(i,j+1)
          pjm   = pvar(i,j  )
          pjp4  = pvar(i,j+2)
          pjm4  = pvar(i,j-1)
c
        one_ze = ordarv(i,j)
c
        const1 = one_ze*  9./8. + (1.-one_ze)*1.
        const2 = one_ze*(-1./8.) 
c
          dv(i,j)= 
     1 const1 *    (
     1 dxi(i)*((-uvip + vi    *dv1ip)-(-uvim + vi     *dv1im))+
     2 dyhi(j)*((-vvjp +vi    *dv1jp)-(-vvjm+vi       *dv1jm)-pjp + pjm)
     1             )
                   
          enddo
        enddo
      return
      end
