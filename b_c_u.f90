subroutine b_c_u(u_old,U_w,U_e,U_s,U_n,imax,jmax)

integer imax,jmax
!real A_u(1:imax+2,1:jmax+1,count_star),D_u(1:imax+2,1:jmax+1,count_star),u_old(1:imax+2,1:jmax+1)
!real A_v(1:imax+1,1:jmax+2,count_star),D_v(1:imax+1,1:jmax+2,count_star),v_old(1:imax+1,1:jmax+2)
real u_old(1:imax+2,1:jmax+1)
real U_w,U_e,U_s,U_n;
!v_star(1:jmax+1,1:imax+2)
!real dt,dx,dy,nu
!real pguess(1:imax+2,1:jmax+2)
m_u=size(u_old,dim=1);
n_u=size(u_old,dim=2);
!m_v=size(v_old,dim=1);
!n_v=size(v_old,dim=2);
!function uold=b_c_u()
u_old(2:m_u-1,1)=U_w; !U_west
!%uold(2,2:end-1)=2*U_w-uold(1,2:end-1); %test-free slip
!%uold(2,2:end-1)=2*U_s-uold(1,2:end-1); % U_s=0 no slip
u_old(1,2:n_u-1)=u_old(2,2:n_u-1); !free slip
!%uold(end-1,2:end-1)=2*U_w-uold(end,2:end-1); %test -free slip
!%uold(end-1,2:end-1)=2*U_s-uold(end,2:end-1); %U_n=0 no slip
u_old(m_u,2:n_u-1)=u_old(m_u-1,2:n_u-1); !free slip
u_old(2:m_u-1,n_u)=u_old(2:m_u-1,n_u-1); !du_e/dx=0;
!%uold(2:end-1,end)=uold(2:end-1,end-1) %du_e/dy=0; outlet bc
end
