subroutine second_correction(u_star,v_star,phi_x,phi_y,u_solid,v_solid,u_d_star,v_d_star,imax,jmax)

integer imax,jmax,i,j;
!real A_u(1:imax+2,1:jmax+1,count_star),D_u(1:imax+2,1:jmax+1,count_star),u_old(1:imax+2,1:jmax+1)
!real A_v(1:imax+1,1:jmax+2,count_star),D_v(1:imax+1,1:jmax+2,count_star),v_old(1:imax+1,1:jmax+2)
real u_star(1:imax+2,1:jmax+1),v_star(1:imax+1,1:jmax+2)
real dt,dx,dy,nu
real u_d_star(1:imax+2,1:jmax+1),v_d_star(1:imax+1,1:jmax+2)
real phi_x(1:imax+2,1:jmax+1),phi_y(1:imax+1,1:jmax+2)
real u_solid(1:imax+2,1:jmax+1),v_solid(1:imax+1,1:jmax+2)

!print*,'entered into second correction';
!function [u_d_star,v_d_star]=second_correction(u_star,v_star,phi_x,phi_y,u_s,v_s1)

u_d_star=u_star;
v_d_star=v_star;
!print*,'u_d_star',size(u_d_star,dim=1),size(u_d_star,dim=2)
!print*,'u_star',size(u_star,dim=1),size(u_star,dim=2)
!print*,'u_solid',size(u_solid,dim=1),size(u_solid,dim=2)
!print*,'phi_x',size(phi_x,dim=1),size(phi_x,dim=2)
!print*,imax,jmax;
!print*,'equalizing done';

!do j=1,jmax+1
! do i=1,imax+2
!print*,'entered for loop';
!print*,i,j,u_star(i,j)
!u_d_star(i,j)=u_star(i,j)-(phi_x(i,j)*u_star(i,j))+(u_solid(i,j)*phi_x(i,j));
 !end do 
 !end do
!print*,'u_d_star computed';
u_d_star=u_star-(phi_x*u_star)+(u_solid*phi_x);
v_d_star=v_star-(phi_y*v_star)+(v_solid*phi_y);
!u_d_star(2:end-1,2:end-1)=u_star(2:end-1,2:end-1)-(pnew(2:end-1,3:end-1)-pnew(2:end-1,2:end-2)).*(dt/dx);
            
!v_d_star(2:end-1,2:end-1)=v_star(2:end-1,2:end-1)-(pnew(3:end-1,2:end-1)-pnew(2:end-2,2:end-1)).*(dt/dy);



end subroutine
