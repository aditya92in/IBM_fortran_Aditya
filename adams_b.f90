subroutine adams_b(imax,jmax,A_u,D_u,A_v,D_v,u_old,v_old,dt,dx,dy,nu,count_star,count_star_temp,pguess,u_star,v_star)

integer imax,jmax,count_star,count_star_temp;
real A_u(1:imax+2,1:jmax+1,count_star),D_u(1:imax+2,1:jmax+1,count_star),u_old(1:imax+2,1:jmax+1)
real A_v(1:imax+1,1:jmax+2,count_star),D_v(1:imax+1,1:jmax+2,count_star),v_old(1:imax+1,1:jmax+2)
real u_star(1:imax+2,1:jmax+1),v_star(1:imax+1,1:jmax+2)
real dt,dx,dy,nu
real pguess(1:imax+2,1:jmax+2)
m_u=size(u_old,dim=1);
n_u=size(u_old,dim=2);
m_v=size(v_old,dim=1);
n_v=size(v_old,dim=2);
m_p=size(pguess,dim=1);
n_p=size(pguess,dim=2);
if(count_star==1) then
!!!!!!!!!!!!!u-star
u_star(2:m_u-1,2:n_u-1)=u_old(2:m_u-1,2:n_u-1)- &
(dt/(dx*dy))*(A_u(2:m_u-1,2:n_u-1,count_star_temp))+ &
(nu*dt/(dx*dy))*(D_u(2:m_u-1,2:n_u-1,count_star_temp)) &
-(pguess(2:m_p-1,3:n_p-1)-pguess(2:m_p-1,2:n_p-2))*(dt/dx);
!!!!!!!!! v-star
v_star(2:m_v-1,2:n_v-1)=v_old(2:m_v-1,2:n_v-1) &
-(dt/(dx*dy))*(A_v(2:m_v-1,2:n_v-1,count_star_temp))&
+(nu*dt/(dx*dy))*(D_v(2:m_v-1,2:n_v-1,count_star_temp)) &
-(pguess(3:m_p-1,2:n_p-1)-pguess(2:m_p-2,2:n_p-1))*(dt/dy);
else if(count_star==2)  then
!!!!!!!!!!!u-star
u_star(2:m_u-1,2:n_u-1)=u_old(2:m_u-1,2:n_u-1)- &
(dt/(dx*dy))*(1.5*A_u(2:m_u-1,2:n_u-1,count_star_temp) &
-0.5*A_u(2:m_u-1,2:n_u-1,count_star_temp-1))+ &
(nu*dt/(dx*dy))*(1.5*D_u(2:m_u-1,2:n_u-1,count_star_temp) &
-0.5*(D_u(2:m_u-1,2:n_u-1,count_star_temp-1))) &
 -(pguess(2:m_p-1,3:n_p-1)-pguess(2:m_p-1,2:n_p-2))*(dt/dx);
!!!!!! v-star
v_star(2:m_v-1,2:n_v-1)=v_old(2:m_v-1,2:n_v-1) &
-(dt/(dx*dy))*(1.5*A_v(2:m_v-1,2:n_v-1,count_star_temp) &
-0.5*A_v(2:m_v-1,2:n_v-1,count_star_temp-1)) &
+(nu*dt/(dx*dy))*(1.5*D_v(2:m_v-1,2:n_v-1,count_star_temp) &
-0.5*D_v(2:m_v-1,2:n_v-1,count_star_temp-1)) &
-(pguess(3:m_p-1,2:n_p-1)-pguess(2:m_p-2,2:n_p-1))*(dt/dy);

    else
!!!!!!!!!u-star
   u_star(2:m_u-1,2:n_u-1)=u_old(2:m_u-1,2:n_u-1)- &
(dt/(dx*dy))*((23.0/12.0)*A_u(2:m_u-1,2:n_u-1,count_star_temp)- &
(4.0/3.0)*A_u(2:m_u-1,2:n_u-1,count_star_temp-1)+ &
(5.0/12.0)*A_u(2:m_u-1,2:n_u-1,count_star_temp-2))+ &
(nu*dt/(dx*dy))*((23.0/12.0)*D_u(2:m_u-1,2:n_u-1,count_star_temp) &
-(4.0/3.0)*(D_u(2:m_u-1,2:n_u-1,count_star_temp-1)) &
+(5.0/12.0)*D_u(2:m_u-1,2:n_u-1,count_star_temp-2)) &
-(pguess(2:m_p-1,3:n_p-1)-pguess(2:m_p-1,2:n_p-2))*(dt/dx);  
   !!!!!!!!! v-star                      
   v_star(2:m_v-1,2:n_v-1)=v_old(2:m_v-1,2:n_v-1)- &
(dt/(dx*dy))*((23.0/12.0)*A_v(2:m_v-1,2:n_v-1,count_star_temp)- &
(4.0/3.0)*A_v(2:m_v-1,2:n_v-1,count_star_temp-1)+ &
(5.0/12.0)*A_v(2:m_v-1,2:n_v-1,count_star_temp-2))&
            +(nu*dt/(dx*dy))*((23.0/12.0)*D_v(2:m_v-1,2:n_v-1,count_star_temp)&
-(4.0/3.0)*(D_v(2:m_v-1,2:n_v-1,count_star_temp-1)) &
+(5.0/12.0)*D_v(2:m_v-1,2:n_v-1,count_star_temp-2)) &
-(pguess(3:m_p-1,2:n_p-1)-pguess(2:m_p-2,2:n_p-1))*(dt/dy);
        
end if

end subroutine

