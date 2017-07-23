subroutine get_rhs(u_star,v_star,dx,dy,dt,imax,jmax,rhs)
integer imax,jmax;
!real A_u(1:imax+2,1:jmax+1,count_star),D_u(1:imax+2,1:jmax+1,count_star),u_old(1:imax+2,1:jmax+1)
!real A_v(1:imax+1,1:jmax+2,count_star),D_v(1:imax+1,1:jmax+2,count_star),v_old(1:imax+1,1:jmax+2)
real u_star(1:imax+2,1:jmax+1),v_star(1:imax+1,1:jmax+2)
real rhs(1:imax+2,1:jmax+2)
real dt,dx,dy,nu
!real pguess(1:imax+2,1:jmax+2)
m_u=size(u_star,dim=1);
n_u=size(u_star,dim=2);
m_v=size(v_star,dim=1);
n_v=size(v_star,dim=2);
!function rhs=get_rhs(u_star,v_star,dx,dy,dt)
rhs=(u_star(2:m_u-1,2:n_u)-u_star(2:m_u-1,1:n_u-1))*(1/(dt*dx))+ (v_star(2:m_v,2:n_v-1)-v_star(1:m_v-1,2:n_v-1))*(1/(dt*dy));
end
