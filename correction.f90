subroutine correction(u_d_star,v_d_star,pguess,dt,dx,dy,imax,jmax,u_new,v_new,res_u,res_v)
!function [unew,vnew,err_corr_u,err_corr_v]= correction(u_d_star,v_d_star,pnew,dt,dx,dy)
integer imax,jmax
!real A_u(1:imax+2,1:jmax+1,count_star),D_u(1:imax+2,1:jmax+1,count_star),u_old(1:imax+2,1:jmax+1)
!real A_v(1:imax+1,1:jmax+2,count_star),D_v(1:imax+1,1:jmax+2,count_star),v_old(1:imax+1,1:jmax+2)
real u_d_star(1:imax+2,1:jmax+1),v_d_star(1:imax+1,1:jmax+2)
real u_new(1:imax+2,1:jmax+1),v_new(1:imax+1,1:jmax+2)
!real
real dt,dx,dy
real pguess(1:imax+2,1:jmax+2)
m_u=size(u_d_star,dim=1);
n_u=size(u_d_star,dim=2);
m_v=size(v_d_star,dim=1);
n_v=size(v_d_star,dim=2);
m_p=size(pguess,dim=1);
n_p=size(pguess,dim=2);
    !unew=u_star;
    u_new=u_d_star;
    v_new=v_d_star;
     !vnew=v_star;
     u_new(2:m_u-1,2:n_u-1)=u_d_star(2:m_u-1,2:n_u-1) &
-(pguess(2:m_p-1,3:n_p-1)-pguess(2:m_p-1,2:n_p-2))*(dt/dx);
     v_new(2:m_v-1,2:n_v-1)=v_d_star(2:m_v-1,2:n_v-1) &
-(pguess(3:m_p-1,2:n_p-1)-pguess(2:m_p-2,2:n_p-1))*(dt/dy);

!print*,'enter res_u'
res_u=norm2(reshape(u_new(2:m_u-1,2:n_u-1),(/1,imax*(jmax-1)/)) &
-reshape(u_d_star(2:m_u-1,2:n_u-1),(/1,imax*(jmax-1)/)));
!print*,'exit res_u';
!print*,'enter res_v';
res_v=norm2(reshape(v_new(2:m_v-1,2:n_v-1),(/1,(imax-1)*jmax/)) &
-reshape(v_d_star(2:m_v-1,2:n_v-1),(/1,(imax-1)*jmax/)));
!print*,'exit res_v';
!unew(2:end-1,2:end-1)=u_d_star(2:end-1,2:end-1)+fnew_x(2:end-1,2:end-1).*dt;
!vnew(2:end-1,2:end-1)=v_d_star(2:end-1,2:end-1)+fnew_y(2:end-1,2:end-1).*dt;


!disp([' error in unew and u** ',num2str(err_corr_u),'error in vnew and v**',num2str(err_corr_v)]);
!pause(1);
!boundary conditions
! unew(1,:)=-unew(2,:);
! unew(end,:)=-unew(end-1,:);
 
! %boundary conditions
!  vnew(1,:)=0;          
!  vnew(end,:)=0;
!  vnew(:,1)=-vnew(:,2);
!  vnew(:,end)=vnew(:,end-1);
end subroutine
