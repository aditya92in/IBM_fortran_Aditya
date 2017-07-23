subroutine divergence(u_new,v_new,dx,dy,imax,jmax,div)

integer imax,jmax;
!real A_u(1:imax+2,1:jmax+1,count_star),D_u(1:imax+2,1:jmax+1,count_star),u_old(1:imax+2,1:jmax+1)
!real A_v(1:imax+1,1:jmax+2,count_star),D_v(1:imax+1,1:jmax+2,count_star),v_old(1:imax+1,1:jmax+2)
real u_new(1:imax+2,1:jmax+1),v_new(1:imax+1,1:jmax+2)
real div(1:imax+2,1:jmax+2);
real dx,dy
!real pguess(1:imax+2,1:jmax+2)
m_u=size(u_new,dim=1);
n_u=size(u_new,dim=2);
m_v=size(v_new,dim=1);
n_v=size(v_new,dim=2);


!function div=divergence(unew,vnew,dx,dy)

div=(u_new(2:m_u-1,2:n_u)-u_new(2:m_u-1,1:n_u-1))*(1/(dx))+ (v_new(2:m_v,2:n_v-1)-v_new(1:m_v-1,2:n_v-1))*(1/(dy));
end subroutine
