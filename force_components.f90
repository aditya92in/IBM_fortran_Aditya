subroutine force_components(imax,jmax,u_d_star,v_d_star,u_star,v_star,dt,dx,dy,rho,F_x,F_y)

!function [F_x,F_y]=force_components(u_d_star,v_d_star,u_star,v_star,dt,dx,dy,rho)
integer imax,jmax,count_star,i,j;
!real A_u(1:imax+2,1:jmax+1,count_star),D_u(1:imax+2,1:jmax+1,count_star),u_old(1:imax+2,1:jmax+1)
!real A_v(1:imax+1,1:jmax+2,count_star),D_v(1:imax+1,1:jmax+2,count_star),v_old(1:imax+1,1:jmax+2)
real u_star(1:imax+2,1:jmax+1),v_star(1:imax+1,1:jmax+2)
real u_d_star(1:imax+2,1:jmax+1),v_d_star(1:imax+1,1:jmax+2)
real dt,dx,dy,temp_u,temp_v,F_x,F_y;
!real pguess(1:imax+2,1:jmax+2)
temp_u=0.0;
temp_v=0.0;

do j=1,jmax+1
 do i=1,imax+2
   temp_u=temp_u + (u_d_star(i,j)-u_star(i,j))*(1/dt);
    end do 
end do 

do j=1,jmax+2
 do i=1,imax+1
!f_x=(u_d_star-u_star)*(1/dt);
temp_v=temp_v+ (v_d_star(i,j)-v_star(i,j))*(1/dt);
end do
end do

F_x=-(temp_u)*rho*dx*dy;
F_y=-(temp_v)*rho*dx*dy;

end subroutine
