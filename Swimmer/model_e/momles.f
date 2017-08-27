      subroutine momxles(du,u,v,
     &               ib,ie,jb,je,
     $                vi,pvar)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer ib,ie,jb,je
      integer ip,im,jp,jm
      real du(1-ih:imax+ih,1-jh:jmax+jh),
     $     u   (1-ih:imax+ih,1-jh:jmax+jh),
     $     v   (1-ih:imax+ih,1-jh:jmax+jh),
     $     pvar(1-ih:imax+ih,1-jh:jmax+jh),
     $     vi(1-ih:imax+ih,1-jh:jmax+jh)
      real uuip ,uuim ,uvjp ,uvjm , vyp,vym
      real du1ip,du1im,du1jp,du1jm
c
      real pip , pim
c
      real fy(1-jh:jmax+jh)
c
      fy = 0.5
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


          du1ip = 2.*dxi(ip)  * (U(ip,j)    -  U(i ,j) )
          du1im = 2.*dxi(i)   * (U(i ,j)    -  U(im,j) )
          du1jp = dyhi(j)     * (U(i,jp)    -  U(i, j) ) +
     &            dxhi(i)     * (V(ip,j)    -  V(i,j)  )
          du1jm = dyhi(jm)    * (U(i,j )    -  U(i,jm) ) +
     &            dxhi(i)     * (V(ip,jm)   -  V(i,jm) )
c
      vyp =
     &   0.5*(
     &    fy(j)     *(vi(i,j)   + vi(ip,j))    +
     &   (1.-fy(j)) *(vi(ip,jp) + vi(i,jp))
     &       )
      vym =
     &   0.5*(
     &    fy(jm)    * (vi(ip,jm)+ vi(i,jm))    +
     &   (1.-fy(jm))* (vi(i,j)  + vi(ip,j))   
     &       )
c
c
          pip   = 0.!pvar(i+1,j)
          pim   = 0.!pvar(i  ,j)
c
         du(i,j) = 2.+
     1  dxhi(i)*( -uuip + vi(ip,j)*du1ip+uuim - vi(i,j)*du1im -pip + pim )+
     2  dyi(j)*( -uvjp + vyp*du1jp+uvjm - vym*du1jm )
            enddo
         enddo
      return
      end

      subroutinemomyles(dv,u,v,
     &                ib,ie,jb,je,
     $                vi,pvar)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer ib,ie,jb,je
      integer im,ip,jm,jp
      real dv(1-ih:imax+ih,1-jh:jmax+jh),
     $     u   (1-ih:imax+ih,1-jh:jmax+jh),
     $     v   (1-ih:imax+ih,1-jh:jmax+jh),
     $     w   (1-ih:imax+ih,1-jh:jmax+jh),
     $     pvar(1-ih:imax+ih,1-jh:jmax+jh),
     $     vi(1-ih:imax+ih,1-jh:jmax+jh)
      real uvip ,uvim ,vvjp ,vvjm 
      real dv1ip,dv1im,dv1jp,dv1jm
      real vxp,vxm
c
      real pjp , pjm
      real fx(1-ih:imax+ih)
c
      fx = 0.5
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

          dv1ip = dxhi(i)     * ( V(ip,j)   - V(i ,j)  ) + 
     &            dyhi(j)     * ( U(i,jp)   - U(i,j)   )
          dv1im = dxhi(im)    * ( V(i ,j)   - V(im,j)  ) +
     &            dyhi(j)     * ( U(im,jp)  - U(im,j)  )
          dv1jp = 2.*dyi(jp)  * ( V(i,jp)   - V(i,j )  )
          dv1jm = 2.*dyi(j)   * ( V(i,j )   - V(i,jm)  )
      vxp =
     &   0.5*(
     &  fx(i)     * (vi(i,j)  + vi(i,jp) )     +
     & (1.-fx(i)) * (vi(ip,j) + vi(ip,jp))
     &       )
      vxm =
     &   0.5*(
     & fx(im)     * (vi(im,j) + vi(im,jp))       +
     & (1.-fx(im))* (vi(i,jp) + vi(i,j)) 
     &       )
C
          pjp   = 0.!pvar(i,j+1)
          pjm   = 0.!pvar(i,j  ) 

          dv(i,j)= 
     1 dxi(i)*((-uvip + vxp         *dv1ip)-(-uvim + vxm    *dv1im))+
     2 dyhi(j)*((-vvjp +vi(i,jp)*dv1jp)-(-vvjm+vi(i,j)*dv1jm)-pjp + pjm)
          enddo
        enddo
      return
      end
