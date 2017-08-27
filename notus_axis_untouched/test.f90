program test
use mod_partial_volume_fractions
implicit none 
real,allocatable,dimension(:)::z,r,z_general,r_general;
real, allocatable, dimension(:,:)::phi; 
integer ilow,ihigh,jlow,jhigh,n,imax,jmax,m,i,j,count_star,n_time_steps,k;
real dz,dr,L_z,L_r,nu,dt,F_x,F_y,a1,a2;
real c_z,c_r,h,R_z,R_r,vol_tot,area_tot;
real exact_surface,error,percent_error; !circular IB
!real x_c,y_c,wi,he,area_tot;      !square IB, ENABLE if you want a square and comment out circular 
real, parameter :: pi = 3.1415927

!--------------------------PART -1 ----------------------------------------------
imax=1152;!no of INTERIOR grid points in the z direction
jmax=960; !no of INTERIOR grid points in the r direction
R_z=1; ! Horizontal radius of the IB
R_r=1;! Vertical radius of the IB  (Useful if we change the IB to an ellipse -(case of deforming cylinder)
h=0.001; ! step sizing for the number of points DEFINING the IB 
n=pi/(h); !number of points DEFINING the IB ! CHANGE TO 2*PI for the torus, and PLUG BACK pi for the SPHERE
 L_z=30; !length of the domain in the z direction (z is horizontal )
 L_r=36; !length of the domain in the r direction (r is the radius)  
 c_z=L_z/2; !z coordinate of centroid of the IB
 c_r=0; !r coordinate of centroid of the IB
dz=L_z/(jmax);
dr=L_r/(imax);
!exact_surface=(4.0/3.0)*pi*(R_r**3); !sphere
exact_surface=(4.0/3.0)*pi*(R_r*R_z*R_r) !Ellipsoid
!exact_surface=2.0*abs(c_r)*r_R**2*pi**2 !torus
!wi=10; !WIDTH of the square IB 
!he=1; !HEIGHT of the square IB

count_star=1;!not used here 
allocate (phi(1:imax+2,1:jmax+2),z(1:jmax+2),r(1:imax+2),z_general(1:n),r_general(1:n));
!phi = volume fractions defined at the cell centres 
!phi (could also be defined at the velocity nodes)
! The module used to compute phi in the latter case remains same.
z(:)=0.0;r(:)=0.0;phi(:,:)=0.0;
!pguess(:,:)=0.0;pcorr(:,:)=0.0;
z(1)=-dz/2; r(1)=-dr/2;
!------------------------------z and r (cell centers values)----------------------------------- 
do m=2,jmax+2
z(m)=z(m-1)+dz; !z- array
end do
do m=2,imax+2
r(m)=r(m-1)+dr; !r-array
end do
!-------------------------------x_general and y_general for a cirlce-------------------------------
do i=1,n
z_general(i)=c_z+R_z*cos(i*h);
r_general(i)=c_r+R_r*sin(i*h);
end do
!z_general(n+1)=z_general(n);
!r_general(n+1)=r_general(1);
!print*,'z_general',z_general;
!print*,'r_generaĺ',r_general
!-----------UNCOMMENT IF you want a SQUARE IB----------------------
!x_general=[x_c-(wi/2), x_c+(wi/2), x_c+(wi/2), x_c-(wi/2)];
!y_general=[y_c-(he/2), y_c-(he/2), y_c+(he/2), y_c+(he/2)];
!n=4;
!note, vertices of the square should be defined in anti-clockwise dirn

!-------------------------Finding Bounds for the Bounding Box Circle /Ellipse----------------------------------------
ilow=maxloc(r,dim=1,mask=(r<c_r-R_r));ihigh=minloc(r,dim=1,mask=(r>c_r+R_r));
jlow=maxloc(z,dim=1,mask=(z<c_z-R_z));jhigh=minloc(z,dim=1,mask=(z>c_z+R_z));
!-------------------------Finding Bounds Square---------------------------
!---------UNCOMMENT if you want a square IB
!ilow=maxloc(y,dim=1,mask=(y<y_c-(he/2)));ihigh=minloc(y,dim=1,mask=(y>y_c+(he/2)));
!jlow=maxloc(x,dim=1,mask=(x<x_c-(wi/2)));jhigh=minloc(x,dim=1,mask=(x>x_c+(wi/2)));

!--------------------------Setting Phi------------------------------------------
!print*,'ihigh',ihigh;
print*,'calling routine for computing  partial volume fractions';
call set_phi(z,r,ilow,ihigh,jlow,jhigh,dz,dr,phi,z_general,r_general,n,c_z,c_r);
print*,'end of computing the volume fractions';
vol_tot=0.0;
area_tot=0.0;
do j=1,jmax+2
 do i=1,imax+2
   !area_tot=area_tot+phi(i,j)*dr*dz;
   vol_tot=vol_tot+phi(i,j)*2*pi*r(i)*dr*dz;
end do
end do
!print*,'area_tot=',area_tot;
print*,'vol_tot=',vol_tot; !check for area 
!print*,'exact surface',2.*abs(c_r)*r_R**2*pi**2 !torus 
!print*,'error=',2.*abs(c_r)*r_R**2*pi**2 -(vol_tot)
print*,'exact surface',exact_surface
print*,'error=',exact_surface-vol_tot!sphere
print*,'percentage error',(exact_surface-vol_tot)*100.0/(exact_surface)
!print*,'exact surface',(2.0/3.0)*pi*R_r**3;
!print*,'error=',((2.0/3.0)*pi*R_r**3)-(vol_tot)!sphere
end program test
