program main
use mod_partial_volume_fractions
!use mod_momentum
integer ilow,ihigh,jlow,jhigh,n,imax,jmax,m,i,j,count_star,n_time_steps,count_star_temp;
parameter(imax=576,jmax=480,n_time_steps=20000,count_star=3)
real L_x,L_y;
real A_u(1:imax+2,1:jmax+1,count_star),D_u(1:imax+2,1:jmax+1,count_star),u_old(1:imax+2,1:jmax+1)
real A_v(1:imax+1,1:jmax+2,count_star),D_v(1:imax+1,1:jmax+2,count_star),v_old(1:imax+1,1:jmax+2)
real u_star(1:imax+2,1:jmax+1),v_star(1:imax+1,1:jmax+2)
real u_new(1:imax+2,1:jmax+1),v_new(1:imax+1,1:jmax+2)
real u_d_star(1:imax+2,1:jmax+1),v_d_star(1:imax+1,1:jmax+2)
real F_x,F_y,c_d,c_l,rho,res_p,res_u,res_v,res_x,res_y;
real U_s,U_e,U_w,U_n;
real V_s,V_e,V_w,V_n;
real u_solid(1:imax+2,1:jmax+1),v_solid(1:imax+1,1:jmax+2)
!real c_d_arr(1:n_time_steps),c_l_arr(1:n_time_steps)
real div(1:imax+2,1:jmax+2),rhs(1:imax+2,1:jmax+2)
!real dt,dx,dy,nu,m,rho;
!real x_c,y_c
real pguess(1:imax+2,1:jmax+2),pcorr(1:imax+2,1:jmax+2)
!real x(1:imax+2),y(1:jmax+2)
real,allocatable,dimension(:)::x,y,x_general,y_general;
real, allocatable, dimension(:,:)::phi,phi_x,phi_y; 
real dx,dy,nu,dt
real x_c,y_c,h,r_x,r_y,area_tot,d_general; !circular IB
!real x_c,y_c,wi,he,area_tot;      !square IB, ENABLE if you want a square and comment out circular 
real, parameter :: pi = 3.1415927
!--------------------------PART -1 ----------------------------------------------
!imax=1152;!no of INTERIOR grid points in the Y direction
!jmax=960; !no of INTERIOR grid points in the X direction
r_x=0.5; ! Horizontal radius of the IB
r_y=0.5 ! Vertical radius of the IB  (Useful if we change the IB to an ellipse -(case of deforming cylinder)
h=0.01; ! step sizing for the number of points DEFINING the IB 
n=int(2*pi/(h)); !number of points DEFINING the IB
L_x=30; !length of the domain in the X direction
L_y=36; !length of the domain in the Y direction
x_c=L_x/2; !x coordinate of centroid of the IB
y_c=L_y/2; !y coordinate of centroid of the IB
dx=L_x/(jmax);
dy=L_y/(imax);
dt=0.01;
rho=1;nu=0.01;
d_general=2*r_x;
!wi=10; !WIDTH of the square IB 
!he=1; !HEIGHT of the square IB
count_star_temp=1;!not used here 
U_s=0.0;U_n=0.0;U_w=1.0;U_e=0.0;
V_s=0.0;V_n=0.0;V_w=0.0;V_e=0.0;
u_solid(:,:)=0.0;v_solid(:,:)=0.0;
u_d_star(:,:)=0.0;v_d_star(:,:)=0.0;
u_new(:,:)=U_w; v_new(:,:)=0.0;
A_u(:,:,:)=0.0;
D_u(:,:,:)=0.0;
D_v(:,:,:)=0.0;
D_u(:,:,:)=0.0;
allocate (phi(1:imax+2,1:jmax+2),x(1:jmax+2),y(1:imax+2),x_general(1:n),y_general(1:n));
m_phi=size(phi,dim=1);
n_phi=size(phi,dim=2);

!phi = volume fractions defined at the cell centres 
!phi (could also be defined at the velocity nodes)
! The module used to compute phi in the latter case remains same.
x(:)=0.0;y(:)=0.0;phi(:,:)=0.0;
pguess(:,:)=0.0;pcorr(:,:)=0.0;
x(1)=-dx/2; y(1)=-dy/2;
!------------------------------x and y (cell centers values)----------------------------------- 
do m=2,jmax+2
x(m)=x(m-1)+dx; !x- array
end do
do m=2,imax+2
y(m)=y(m-1)+dy; !y-array
end do
!-------------------------------x_general and y_general for a cirlce-------------------------------
do i=1,n
x_general(i)=x_c+r_x*cos(i*h);
y_general(i)=y_c+r_y*sin(i*h);
end do

!-----------UNCOMMENT IF you want a SQUARE IB----------------------
!x_general=[x_c-(wi/2), x_c+(wi/2), x_c+(wi/2), x_c-(wi/2)];
!y_general=[y_c-(he/2), y_c-(he/2), y_c+(he/2), y_c+(he/2)];
!n=4;
!note, vertices of the square should be defined in anti-clockwise dirn

!-------------------------Finding Bounds for the Bounding Box Circle /Ellipse----------------------------------------
ilow=maxloc(y,dim=1,mask=(y<y_c-r_y));ihigh=minloc(y,dim=1,mask=(y>y_c+r_y));
jlow=maxloc(x,dim=1,mask=(x<x_c-r_x));jhigh=minloc(x,dim=1,mask=(x>x_c+r_x));
!-------------------------Finding Bounds Square---------------------------
!---------UNCOMMENT if you want a square IB
!ilow=maxloc(y,dim=1,mask=(y<y_c-(he/2)));ihigh=minloc(y,dim=1,mask=(y>y_c+(he/2)));
!jlow=maxloc(x,dim=1,mask=(x<x_c-(wi/2)));jhigh=minloc(x,dim=1,mask=(x>x_c+(wi/2)));

!--------------------------Setting Phi------------------------------------------
print*,'calling routine for computing  partial volume fractions';
call set_phi(x,y,ilow,ihigh,jlow,jhigh,dx,dy,phi,x_general,y_general,n);
print*,'end of computing the volume fractions';
do j=1,jmax+2
 do i=1,imax+2
   area_tot=area_tot+phi(i,j)*dx*dy;
end do
end do
print*,'area_tot=',area_tot; !check for area 

allocate(phi_x(1:imax+2,1:jmax+1),phi_y(1:imax+1,1:jmax+2))
phi_x=0.5*(phi(1:imax+2,1:jmax+1)+phi(1:imax+2,2:jmax+2));
phi_y=0.5*(phi(1:imax+1,1:jmax+2)+phi(2:imax+2,1:jmax+2));



!end program test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!--Main Loop!!!!!!!!!!!!!!!!!
open(unit=5,file='C_d_C_l.txt')
open(unit=1,file='residuals.txt')
do i=1,20
  u_old=u_new;
  v_old=v_new;
!print*,'enter bcu';
!%**************apply boundary conditions for uold
!print*,'U_W',U_w;
call b_c_u(u_old,U_w,U_e,U_s,U_n,imax,jmax);
!print*,'exit bcu';
!%***************apply boundary conditions for vold
!print*,'enter bcv';
call b_c_v(v_old,V_w,V_e,V_s,V_n,imax,jmax);
!print*,'exit bcv';
print*,'enter momentum';
u_old(:,:)=0.0;
v_old(:,:)=0.0;
call momentum (u_old,v_old,A_u,D_u,A_v,D_v,dx,dy,imax,jmax,count_star,count_star_temp)
print*,'exit momentum';
print*,'enter adams_b';
call adams_b(imax,jmax,A_u,D_u,A_v,D_v,u_old,v_old,dt,dx,dy,nu,count_star,count_star_temp,pguess,u_star,v_star)
print*,'u_star(50,50)',u_star(50,50);
print*,'exit adams_b';
!print*,'enter bcu_ustar';
call b_c_u(u_star,U_w,U_e,U_s,U_n,imax,jmax);
!print*,'exit bcu_ustar';
!print*,'enter bcv vstar';
call b_c_v(v_star,V_w,V_e,V_s,V_n,imax,jmax);
!print*,'enter bcv vstar';
!print*,'enter secondcorr';
call second_correction(u_star,v_star,phi_x,phi_y,u_solid,v_solid,u_d_star,v_d_star,imax,jmax)
print*,'u_d_star',u_d_star(50,50);
!print*,'exit secondcorr';
!print*,'enter forcecomp';
call force_components(imax,jmax,u_d_star,v_d_star,u_star,v_star,dt,dx,dy,rho,F_x,F_y)
!print*,'exit forcecomp';
!print*,'enter integralquant';
print*,'F_y=',F_y;
call integral_quant(F_x,F_y,d_general,U_w,rho,c_d,c_l)
!print*,'exit integralqunat';
print*,'C_d',c_d,'C_l',c_l;
write(9,*)c_d,c_l;
!print*,'enter get_rhs';
call get_rhs(u_star,v_star,dx,dy,dt,imax,jmax,rhs)
print*,'rhs(50,50)',rhs(50,50);
!print*,'exit get_rhs';
!print*,'enter solve_poisson';
call solve_poisson(pcorr,imax,jmax,dx,dy,res_p,rhs)
!print*,'exit solve_poisson';
!print*,'enter correction';
call correction(u_d_star,v_d_star,pcorr,dt,dx,dy,imax,jmax,u_new,v_new,res_u,res_v)
!print*,'exit correction';
!print*,'enter divergence';
call divergence(u_new,v_new,dx,dy,imax,jmax,div)
!print*,'exit divergence';
!print*,'enter residual';
call residual(u_old,u_new,v_old,v_new,imax,jmax,res_x,res_y)
!print*,'exit residual';
print*,'res_x=',res_x,'res_y=',res_y;
write(10,*) res_x,res_y;

!for PLOTTING Only
![up,vp]=inter_p(unew,vnew);

!CHECKS

!******************check the divergence and courant number 

!*********************compute the resdiuals
![res_x,res_y]=res(uold,unew,vold,vnew);
!%toc;

!psi=STREAMFUNCTION(imax,jmax,dx,dy,unew,vnew);
!vort=VORTICITY(imax,jmax,dx,dy,unew,vnew,U_w);
!vort(vort<0 & vort >-0.1)=0;
!vort(vort>0 & vort <0.1)=0;%*****************post processing
 pguess=pguess+pcorr;


if(count_star_temp>=3) then
        A_u(:,:,count_star_temp-2)=A_u(:,:,count_star_temp-1);
        A_v(:,:,count_star_temp-2)=A_v(:,:,count_star_temp-1);
        D_u(:,:,count_star_temp-2)=D_u(:,:,count_star_temp-1);
        D_v(:,:,count_star_temp-2)=D_v(:,:,count_star_temp-1);
        A_u(:,:,count_star_temp-1)=A_u(:,:,count_star_temp);
        A_v(:,:,count_star_temp-1)=A_v(:,:,count_star_temp);
        D_u(:,:,count_star_temp-1)=D_u(:,:,count_star_temp);
        D_v(:,:,count_star_temp-1)=D_v(:,:,count_star_temp);
        count_star_temp=3;
        
        
else
    count_star_temp=count_star_temp+1;
        
end if
 
end do

 close(9)
 close(10)
end program 

