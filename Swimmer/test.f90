program test
use mod_partial_volume_fractions
implicit none 
real,allocatable,dimension(:)::r,z,r_general,z_general;
real, allocatable, dimension(:,:)::phi_a,phi_v; 
integer ilow,ihigh,jlow,jhigh,n,imax,jmax,m,i,j,count_star,n_time_steps,k;
real dr,dz,L_r,L_z,nu,dt,F_x,F_y,a1,a2;
real r_c,z_c,h,r_x,r_y!circular IB
real act_area,area_tot,vol_tot,vol_act; 
!real x_c,y_c,wi,he,area_tot,vol_tot,vol_act;      !square IB, ENABLE if you want a square and comment out circular 
!real, parameter :: pi = 3.1415927

!--------------------------PART -1 ----------------------------------------------
imax=1000;!no of INTERIOR grid points in the Y direction
jmax=1000; !no of INTERIOR grid points in the X direction
r_x=1.0; ! Horizontal radius of the IB
r_y=1.0 ! Vertical radius of the IB  (Useful if we change the IB to an ellipse -(case of deforming cylinder)
h=0.001; ! step sizing for the number of points DEFINING the IB 
n=4*atan(1.0)/(h); !number of points DEFINING the IB
L_r=30; !length of the domain in the X direction
L_z=36; !length of the domain in the Y direction
r_c=0; !x coordinate of centroid of the IB
z_c=10; !y coordinate of centroid of the IB
dr=L_r/(jmax);
dz=L_z/(imax);
!wi=2; !WIDTH of the square IB 
!he=2; !HEIGHT of the square IB

open(unit=1,file='r_grid.txt');
open(unit=2,file='z_grid.txt');
open(unit=3,file='r_general.txt');
open(unit=4,file='z_general.txt');
count_star=1;!not used here 
allocate (phi_a(1:imax+2,1:jmax+2),phi_v(1:imax+2,1:jmax+2),r(1:jmax+2),z(1:imax+2),r_general(1:n),z_general(1:n));
!phi = volume fractions defined at the cell centres 
!phi (could also be defined at the velocity nodes)
! The module used to compute phi in the latter case remains same.
r(:)=0.0;r(:)=0.0;phi_a(:,:)=0.0;
!pguess(:,:)=0.0;pcorr(:,:)=0.0;
r(1)=-dr/2; z(1)=-dz/2;
act_area=2*atan(1.0)*r_x*r_y;
!act_area=wi*he;
vol_act=(4.0/3.0)*4*atan(1.0)*r_x*r_x*r_x;
!------------------------------x and y (cell centers values)----------------------------------- 
do m=2,jmax+2
r(m)=r(m-1)+dr; !x- array
end do
write(1,*) r
do m=2,imax+2
z(m)=z(m-1)+dz; !y-array
end do
write(2,*) z
!-------------------------------x_general and y_general for a cirlce-------------------------------
r_general(1)=r_c+r_x*cos(-2*atan(1.0));
z_general(1)=z_c+r_y*sin(-2*atan(1.0));
do i=2,n
r_general(i)=r_c+r_x*cos(i*h-2*atan(1.0));
z_general(i)=z_c+r_y*sin(i*h-2*atan(1.0));
end do

!-----------UNCOMMENT IF you want a SQUARE IB----------------------
!x_general=[x_c-(wi/2), x_c+(wi/2), x_c+(wi/2), x_c-(wi/2)];
!y_general=[y_c-(he/2), y_c-(he/2), y_c+(he/2), y_c+(he/2)];
!n=4;
write(3,*)r_general
write(4,*)z_general



!note, vertices of the square should be defined in anti-clockwise dirn

!-------------------------Finding Bounds for the Bounding Box Circle /Ellipse----------------------------------------
ilow=maxloc(z,dim=1,mask=(z<z_c-r_y));
ihigh=minloc(z,dim=1,mask=(z>z_c+r_y));
jlow=1;
jhigh=minloc(r,dim=1,mask=(r>r_c+r_x));
!-------------------------Finding Bounds Square---------------------------
!---------UNCOMMENT if you want a square IB
!ilow=maxloc(y,dim=1,mask=(y<y_c-(he/2)));
!ihigh=minloc(y,dim=1,mask=(y>y_c+(he/2)));
!jlow=maxloc(x,dim=1,mask=(x<x_c-(wi/2)));
!jlow=1;
!jhigh=minloc(x,dim=1,mask=(x>x_c+(wi/2)));

!--------------------------Setting Phi------------------------------------------
print*,'dr=',dr,'dz=',dz;
do k=1,kbmax
call deformation(ibmax,jbmax,kbmax,M,N,a,Tcycle,r_general,z_general,t_small,k) 
!print*,'calling routine for computing  partial volume fractions';
call set_phi(r,z,ilow,ihigh,jlow,jhigh,dr,dz,phi_a,phi_v,r_general,z_general,n);
!print*,'end of computing the volume fractions';
area_tot=0.0;
vol_tot=0.0;
do j=1,jmax+2
 do i=1,imax+2
   area_tot=area_tot+phi_a(i,j)*dr*dz;
   vol_tot=vol_tot+phi_v(i,j)*2*4*atan(1.0)*(abs((r(1)-r(j))))*dr*dz
end do
end do
print*,'area_tot=',area_tot; !check for area
print*,'vol_tot=',vol_tot;
print*,'actual area=',act_area;
print*,'actual vol=',vol_act;
print*,'error (%)=', abs(act_area-area_tot)*100.0/(act_area);
print*,'error (%) in vol=',abs(vol_tot-vol_act)*100.0/(vol_act); 
end program test
