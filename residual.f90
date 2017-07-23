subroutine residual(u_old,u_new,v_old,v_new,imax,jmax,res_x,res_y)

integer imax,jmax;
!real A_u(1:imax+2,1:jmax+1,count_star),D_u(1:imax+2,1:jmax+1,count_star),
real u_old(1:imax+2,1:jmax+1),u_new(1:imax+2,1:jmax+1)
!real A_v(1:imax+1,1:jmax+2,count_star),D_v(1:imax+1,1:jmax+2,count_star),v_old(1:imax+1,1:jmax+2)
real v_new(1:imax+1,1:jmax+2),v_old(1:imax+1,1:jmax+2)
real res_x,res_y
integer m_u,n_u,m_v,n_v!m_p,n_p;
real pguess(1:imax+2,1:jmax+2)


m_u=size(u_old,dim=1);
n_u=size(u_old,dim=2);
m_v=size(v_old,dim=1);
n_v=size(v_old,dim=2);
!m_p=size(pguess,dim=1);
!n_p=size(pguess,dim=2);

res_x=norm2(reshape(u_new(2:m_u-1,2:n_u-1),(/1,imax*(jmax-1)/)) &
-reshape(u_old(2:m_u-1,2:n_u-1),(/1,imax*(jmax-1)/)));
res_y=norm2(reshape(v_new(2:m_v-1,2:n_v-1),(/1,(imax-1)*jmax/)) &
-reshape(v_old(2:m_v-1,2:n_v-1),(/1,(imax-1)*jmax /)));

!function [e_x,e_y]=res(uold,unew,vold,vnew)
!e_x=norm(unew(2:end-1,2:end-1)-uold(2:end-1,2:end-1));
!e_y=norm(vnew(2:end-1,2:end-1)-vold(2:end-1,2:end-1));
end subroutine
