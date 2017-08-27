!---def_test= deformation_test_program -----------------!
program test
use mod_partial_volume_fractions
use mod_deformation
implicit none 
real,allocatable,dimension(:)::z_p,r_p
real, allocatable, dimension(:,:)::phi; 
integer ilow,ihigh,jlow,jhigh,n,imax,jmax,m,i,j,count_star,n_time_steps,k;
integer ibmax,jbmax,kbmax,Mmax,Nmax
parameter(ibmax=80,jbmax=80,kbmax=80,Mmax=80,Nmax=80)
real(kind=8) a,Tcycle,t_small(1:kbmax)
real start,finish,start_1,finish_1
real dz,dr,L_z,L_r,nu,dt,F_x,F_y,a1,a2;
real c_z,c_r,h,R_z,R_r,vol_tot,area_tot;
real exact_surface,error,percent_error; !circular IB
!character(len=15)::str,str2;
!character(len=4)::str1;
!character(len=5)::str_3
real ,allocatable,dimension(:):: z_cyl,radius_cyl
real time_d,time_cv,time_sp,time_tot

open(unit=11,file='radius_cyl.txt'); !radius
open(unit=12,file='z_cyl.txt'); !z 
open(unit=13,file='z_c_bb.txt'); 
open(unit=14,file='r_c_bb.txt');
open(unit=15,file='z_grid.txt');
open(unit=16,file='r_grid.txt');
!--------------------------PART -1 ----------------------------------------------
!.................z---vertical
!.................r---horizontal

imax=100;!no of INTERIOR grid points in the z direction
jmax=100; !no of INTERIOR grid points in the r direction
R_z=1; ! Horizontal radius of the IB
R_r=1! Vertical radius of the IB  (Useful if we change the IB to an ellipse -(case of deforming cylinder)
L_r=10; !length of the domain in the z direction (z is the vertical )
L_z=4.5+5.5; !length of the domain in the r direction (r is the radius)  
 c_z=0; !z coordinate of centroid of the IB
a=1.0; 
Tcycle=0.5; 
dz=L_z/(imax);
dr=L_r/(jmax);
 c_r=L_r/2; !r coordinate of centroid of the IB
exact_surface=(4.0/3.0)*4*atan(1.0)*(a**3); !sphere
count_star=1;!not used here 
allocate (phi(1:imax+2,1:jmax+2),z_p(1:imax+2),r_p(1:jmax+2))
!z_p(:)=0.0;r_p(:)=0.0;phi(:,:)=0.0;
z_p(1)=-5.5; r_p(1)=0;
!------------------------------z and r (cell centers values)----------------------------------- 
do m=2,imax+2
z_p(m)=z_p(m-1)+dz; !z- array
end do
write(15,*) z_p

do m=2,jmax+2
r_p(m)=r_p(m-1)+dr; !r-array
end do
write(16,*) r_p
!-------------------------Bounding Box----------------------------------------
ilow=maxloc(z_p,dim=1,mask=(z_p<c_z-(1.5*R_z)));
ihigh=minloc(z_p,dim=1,mask=(z_p>c_z+(1.5*R_z)));
jlow=1;
jhigh=minloc(r_p,dim=1,mask=(r_p>c_r+(1.1*R_r)));
!----------------------------------------------
write(13,*) z_p(ilow),z_p(ilow),z_p(ihigh),z_p(ihigh),z_p(ilow)
write(14,*) r_p(jlow),r_p(jhigh),r_p(jhigh),r_p(jlow),r_p(jlow)
allocate(z_cyl(1:ibmax*jbmax),radius_cyl(1:ibmax*jbmax))
!-----------------------------------Main time loop----------------------------
do k=1,1
call cpu_time(start_1)
t_small(k) = ((k-1.0)/(kbmax-1.0))*Tcycle
call cpu_time(start)
call deformation(ibmax,jbmax,kbmax,Mmax,Nmax,a,Tcycle,radius_cyl,z_cyl,&
t_small(k))
call cpu_time(finish)
time_d=finish-start;
 write(11,*) radius_cyl
 write(12,*) z_cyl
call cpu_time(start)
call set_phi(r_p,z_p,ilow,ihigh,jlow,jhigh,dr,dz,phi,radius_cyl,z_cyl,ibmax*jbmax);
call cpu_time(finish)
time_sp=finish-start;
!print*,'end of computing the volume fractions';
vol_tot=0.0;
area_tot=0.0;
call cpu_time(start)
do j=1,jmax+2
do i=1,imax+2
   area_tot=area_tot+phi(i,j)*dr*dz;
   vol_tot=vol_tot+phi(i,j)*2*4*atan(1.0)*(abs(r_p(1)-(r_p(j)+dr/2)))*dr*dz
end do
end do
call cpu_time(finish)
time_cv=finish-start;
call cpu_time(finish_1)
time_tot=finish_1-start_1;
print*,'area_tot=',area_tot;
print*,'vol_tot=',vol_tot
print*,'error(%)=',(vol_tot-exact_surface)*100.0/(exact_surface)
print*,'time(%) in deformation routine=',time_d*100.0/time_tot
print*,'time(%) spent in setting phi=',time_sp*100.0/time_tot
!'time spent in computing volume=',time_cv*100.0/time_tot !check for area 
end do
!------------------------- Main Loop Ends----------------------------------------
end program test
