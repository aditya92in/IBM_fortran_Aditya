      subroutine adscal(dc,c,u,v,
     &                ib,ie,jb,je,
     $                vi)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer ib,ie,jb,je
      integer ip,im,jp,jm
      integer iorder
      real dc(1-ih:imax+ih,1-jh:jmax+jh),
     $     c   (1-ih:imax+ih,1-jh:jmax+jh),
     $     u   (1-ih:imax+ih,1-jh:jmax+jh),
     $     v   (1-ih:imax+ih,1-jh:jmax+jh)
      real vi
c
c
      real ucip ,ucim ,vcjp ,vcjm 
      real dc1ip,dc1im,dc1jp,dc1jm
c
c
      
        do j=jb,je
          do i=ib,ie
*
          ip = i + 1
          jp = j + 1
          im = i - 1
          jm = j - 1

          ucip  = 0.5 * U(ip,j)*( C(ip,j)+C(i,j)  )
          ucim  = 0.5 * U(im,j)*( C(im,j)+C(i,j)  )
          vcjp  =  0.5*V(i,jp)*
     &             ( C(i,jp)+C(i,j)  )
          vcjm  =  0.5*V(i,jm)*
     &             ( C(i,jm)+C(i,j))


          dc1ip = dxi(ip) * (C(ip,j)-C(i ,j) )
          dc1im = dxi(i)  * (C(i ,j)-C(im,j) )
          dc1jp = dyi(j)  * (C(i,jp)-C(i, j) )
          dc1jm = dyi(jm) * (C(i,j )-C(i,jm) )


         dc(i,j) =
     1  dxi(i)*(  -ucip + vi*dc1ip+ucim - vi*dc1im )+
     2  dyi(j)*(  -vcjp + vi*dc1jp+vcjm - vi*dc1jm )
            enddo
         enddo
      return
      end
