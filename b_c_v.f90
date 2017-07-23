subroutine b_c_v(v_old,V_w,V_e,V_s,V_n,imax,jmax)

integer imax,jmax
!real A_u(1:imax+2,1:jmax+1,count_star),D_u(1:imax+2,1:jmax+1,count_star),u_old(1:imax+2,1:jmax+1)
!real A_v(1:imax+1,1:jmax+2,count_star),D_v(1:imax+1,1:jmax+2,count_star),v_old(1:imax+1,1:jmax+2)
real v_old(1:imax+1,1:jmax+2)
real V_w,V_e,V_s,V_n;
!v_star(1:jmax+1,1:imax+2)
!real dt,dx,dy,nu
!real pguess(1:imax+2,1:jmax+2)
m_v=size(v_old,dim=1);
n_v=size(v_old,dim=2);

v_old(2:m_v-1,1)=2*V_w-v_old(2:m_v-1,2); !V_w=0;
v_old(2:m_v-1,n_v)=v_old(2:m_v-1,n_v-1); !dv_e/dx=0;
v_old(1,2:n_v-1)=V_s; !V_s=0; no penetration south
v_old(m_v,2:n_v-1)=V_n; !V_n=0; no penetration north
!end
end
