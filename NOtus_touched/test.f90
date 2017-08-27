program test
use mod_partial_volume_fractions
use mod_deformation
implicit none 
real,allocatable,dimension(:)::z,r,z_general,r_general;
real, allocatable, dimension(:,:)::phi; 
integer ilow,ihigh,jlow,jhigh,n,imax,jmax,m,i,j,count_star,n_time_steps,k_current;
integer ibmax,jbmax,kbmax,Mmax,Nmax
parameter(ibmax=200,jbmax=200,kbmax=200,Mmax=200,Nmax=200)
real(kind=8) a,Tcycle,t_small(1:kbmax)
real start,finish,start_1,finish_1
real dz,dr,L_z,L_r,nu,dt,F_x,F_y,a1,a2;
real c_z,c_r,h,R_z,R_r,vol_tot,area_tot;
real exact_surface;
real ,allocatable,dimension(:):: z_cyl,r_cyl,z_rect,r_rect
real time_d,time_cv,time_sp,time_tot
real percent_error(1:kbmax);

!--------------------------PART -1 ----------------------------------------------
imax=576;!no of INTERIOR grid points in the z direction
jmax=488; !no of INTERIOR grid points in the r direction
R_z=1.0; ! Vertical radius of the IB
R_r=1.0;! Horizontal radius of the IB  (Useful if we change the IB to an ellipse -(case of deforming cylinder)
h=0.001; ! step sizing for the number of points DEFINING the IB 
n=4*atan(1.0)/(h); !number of points DEFINING the IB ! CHANGE TO 2*PI for the torus, and PLUG BACK pi for the SPHERE
 L_z=36+5.5; !length of the domain in the z direction (z is vertical )
 L_r=30; !length of the domain in the r direction (r is the horizontal)  
 c_z=0.0; !z coordinate of centroid of the IB
 c_r=0.0; !r coordinate of centroid of the IB
dz=L_z/(imax);
dr=L_r/(jmax);
!exact_surface=(4.0/3.0)*pi*(R_r**3); !sphere
exact_surface=(4.0/3.0)*4*atan(1.0)*(R_r*R_z*R_r) !Ellipsoid
!exact_surface=2.0*abs(c_r)*r_R**2*pi**2 !torus
!wi=10; !WIDTH of the square IB 
!he=1; !HEIGHT of the square IB
a=1.0; 
Tcycle=0.5; 
count_star=1;!not used here 
allocate (phi(1:imax+2,1:jmax+2),z(1:imax+2),r(1:jmax+2));
!phi = volume fractions defined at the cell centres 
!phi (could also be defined at the velocity nodes)
! The module used to compute phi in the latter case remains same.
z(:)=0.0;r(:)=0.0;phi(:,:)=0.0;
!pguess(:,:)=0.0;pcorr(:,:)=0.0;
z(1)=-5.5; r(1)=-dr/2.0;
allocate(z_rect(1:5),r_rect(1:5))
!print*,r(1);
open(unit=11,file='r_cyl.txt');
open(unit=12,file='z_cyl.txt');
open(unit=13,file='z_c_bb.txt'); 
open(unit=14,file='r_c_bb.txt');
open(unit=15,file='z_general.txt');
open(unit=16,file='r_general.txt');
open(unit=17,file='percent_error.txt');
open(unit=18,file='r_grid.txt');
open(unit=19,file='z_grid.txt');
open(unit=20,file='coordinates.txt');
!------------------------------z and r (cell centers values)----------------------------------- 
do m=2,imax+2
z(m)=z(m-1)+dz; !z- array
end do
write(19,*) z

do m=2,jmax+2
r(m)=r(m-1)+dr; !r-array
end do
write(18,*) r


!print*,dr;
!print*,r(1)
!-------------------------------x_general and y_general for a cirlce-------------------------------
allocate(z_cyl(1:Mmax),r_cyl(1:Nmax))
!do i=1,n
!z_general(i)=c_z+R_z*cos(i*h);
!r_general(i)=c_r+R_r*sin(i*h);
!end do
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
jlow=1;
jhigh=minloc(r,dim=1,mask=(r>c_r+1.2*R_r));
ilow=maxloc(z,dim=1,mask=(z<c_z-1.3*R_z));
ihigh=minloc(z,dim=1,mask=(z>c_z+1.3*R_z));
!print*,r(jlow),r(jhigh),r(jhigh),r(jlow),r(jlow)
write(20,*) jlow,jhigh,ilow,ihigh, dr,dz
write(13,*) z(ilow),z(ilow),z(ihigh),z(ihigh),z(ilow)
write(14,*) r(jlow),r(jhigh),r(jhigh),r(jlow),r(jlow)
!write(15,*) z_general
!write(16,*) r_general


!-------------------------Finding Bounds Square---------------------------
!---------UNCOMMENT if you want a square IB
!ilow=maxloc(y,dim=1,mask=(y<y_c-(he/2)));ihigh=minloc(y,dim=1,mask=(y>y_c+(he/2)));
!jlow=maxloc(x,dim=1,mask=(x<x_c-(wi/2)));jhigh=minloc(x,dim=1,mask=(x>x_c+(wi/2)));

!--------------------------Setting Phi------------------------------------------
!print*,'ihigh',ihigh;

do k_current=1,1
call cpu_time(start_1)
t_small(k_current) = ((k_current-1.0)/(kbmax-1.0))*Tcycle
call cpu_time(start)
call deformation(ibmax,jbmax,kbmax,Mmax,Nmax,a,Tcycle,r_cyl,z_cyl,&
t_small(k_current),k_current)
call cpu_time(finish)
time_d=finish-start;
 write(11,*) r_cyl
 write(12,*) z_cyl
call cpu_time(start)

z_rect=[ 0.0, 0.0, 1.0 , 1.0 ,0.0];
r_rect=[ 1.0, 1.0, 1.0, 1.0 ,1.0];
call set_phi(z,r,1,imax+2,1,jmax+2,dz,dr,phi,z_rect,r_rect,5);
call cpu_time(finish)
time_sp=finish-start;
vol_tot=0.0;
area_tot=0.0;
do j=1,jmax+2
 do i=1,imax+2
   area_tot=area_tot+phi(i,j)*dr*dz;
   vol_tot=vol_tot+phi(i,j)*2*4*atan(1.0)*(abs(r(1)-(r(j)+0.5*dr)))*dr*dz;
end do
end do
print*,'area_tot=',area_tot;
print*,'time in(%) of cycle time',(k_current*1.0/kbmax*1.0)*Tcycle*200.0;
print*,'vol_tot=',vol_tot,'exact surface',exact_surface; !check for area 
print*,'error=',exact_surface-vol_tot,'percentage error',abs(exact_surface-vol_tot)*100.0/(exact_surface);
percent_error(k_current)=abs(exact_surface-vol_tot)*100.0/(exact_surface);
!print*,'exact surface',(2.0/3.0)*pi*R_r**3;
!print*,'error=',((2.0/3.0)*pi*R_r**3)-(vol_tot)!sphere
end do
write(17,*) percent_error
end program test
