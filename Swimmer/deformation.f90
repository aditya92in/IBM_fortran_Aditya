module mod_deformation
contains
subroutine deformation(imax,jmax,kmax,M,N,a,Tcycle,radius_cyl,z_cyl,t_small,k_current) 
integer(kind=4) M,N,i,j,k,n_small,exit_status,m_small,k_current;
parameter(n_small=3) 
integer imax,jmax,kmax
!parameter(imax=80,jmax=80,kmax=80) 
real(kind=8) phi(1:N,1:M)
real(kind=8) a,Tcycle,t_small
real  radius_cyl(1:M),z_cyl(1:M)
real(kind=8) rho1(1:N,1:M),rho2(1:N,1:M),rho3(1:N,1:M)
real(kind=8) rho4(1:N,1:M),rho5(1:N,1:M),rho6(1:N,1:M),rho7(1:N,1:M)
real(kind=8) rho8(1:N,1:M),rho9(1:N,1:M),rho10(1:N,1:M),rho11(1:N,1:M)
real(kind=8) rho12(1:N,1:M),rho13(1:N,1:M),rho14(1:N,1:M),rho15(1:N,1:M)
real(kind=8) rho16(1:N,1:M),rho17(1:N,1:M),rho18(1:N,1:M),rho19(1:N,1:M)
real(kind=8) rho20(1:N,1:M),rho21(1:N,1:M),rho22(1:N,1:M)
real(kind=8) R1(1:N,1:M),R2(1:N,1:M),R3(1:N,1:M)
real(kind=8) R4(1:N,1:M),T(1:N,1:M)
character(len=15)::str,str2;
character(len=4)::str1;
character(len=5)::str_3
real(kind=8) omega,eta,mu1c,mu2c,mu3c,kappa2c,kappa3c,pi;
parameter(eta=500.0,mu1c=0.05,mu2c=0,kappa2c=0,mu3c=0)
real(kind=8) mu1s,mu2s,kappa2s,mu3s,kappa3s
parameter(mu1s=0.0,mu2s=0.0,mu3s=0.0,kappa3s=0.0)
real(kind=8) leg1(1:N,1:M),leg2(1:N,1:M),d_phi;
real(kind=8) output(0:n_small);  
str1='pic';
phi(1,1:M)=4*atan(1.0)*0.5; 
phi(N,1:M)=-4*atan(1.0)*0.5;
d_phi=(phi(N,1)-phi(1,1))*(1.0/(N-1));
omega=2*4*atan(1.0)/Tcycle;
kappa3c=-(27.0/59.0)*mu1c
kappa2s=(5.0/3.0)*(sqrt(230.0/413.0))*mu1c
pi=4*atan(1.0);
do i=2,N
 phi(i,1:M)=phi(i-1,1:M)+d_phi;
end do
output(:)=0.0;
do j=1,M
 do i=1,N
call legendre_associated(n_small,0,cos(4*atan(1.0)*0.5-phi(i,j)),output);
leg1(i,j)=output(size(output)-1);
call legendre_associated(n_small,1,cos(4*atan(1.0)*0.5-phi(i,j)),output);
leg2(i,j)=output(size(output)-1);
end do
end do 
U=7*mu1c*kappa2s-7*mu1s*kappa2c+6*kappa2c*kappa3s-6*kappa2s*kappa3c
D=8*pi*eta*(omega**2)*(a**3)*(1.0/420.0)*(630*(mu1c**2+mu1s**2)+ &
567*(kappa2c**2+kappa2s**2)+1180*(kappa3c**2+kappa3s**2))

do j=1,M
 do i=1,N
rho1(i,j)=a*sin(pi/2.0-phi(i,j));
rho2(i,j)=a*cos(pi/2.0-phi(i,j));
rho3(i,j)=mu1c*2*cos(pi/2.0-phi(i,j))*sin(pi/2-phi(i,j))+mu1c*sin(pi/2-phi(i,j))*cos(pi/2-phi(i,j));
rho4(i,j)=mu1c*2*cos(pi/2.0-phi(i,j))*cos(pi/2-phi(i,j))-mu1c*sin(pi/2-phi(i,j))*sin(pi/2-phi(i,j));
rho5(i,j)=mu2c*0.75*(1+3*cos(2*(pi/2-phi(i,j)))*sin(pi/2-phi(i,j)))+ &
mu2c*0.75*2*sin(2*(pi/2-phi(i,j)))*cos(pi/2-phi(i,j));
rho6(i,j)=mu2c*0.75*(1+3*cos(2*(pi/2-phi(i,j)))*cos(pi/2-phi(i,j))) &
-mu2c*0.75*2*sin(2*(pi/2-phi(i,j)))*sin(pi/2-phi(i,j));
rho7(i,j)=kappa2c*0.75*(1+3*cos(2*(pi/2-phi(i,j)))*sin(pi/2-phi(i,j)));
rho8(i,j)=kappa2c*0.75*(1+3*cos(2*(pi/2-phi(i,j)))*cos(pi/2-phi(i,j)));
rho9(i,j)=mu3c*4*leg1(i,j)*sin(pi/2-phi(i,j))+mu3c*leg2(i,j)*cos(pi/2-phi(i,j));
rho10(i,j)=mu3c*4*leg1(i,j)*cos(pi/2-phi(i,j))-mu3c*4*leg2(i,j)*sin(pi/2-phi(i,j));
rho11(i,j)=kappa3c*4*leg1(i,j)*sin(pi/2-phi(i,j))+kappa3c*(1/3)*leg2(i,j)*cos(pi/2-phi(i,j));
rho12(i,j)=kappa3c*4*leg1(i,j)*cos(pi/2-phi(i,j))-kappa3c*(1/3)*leg2(i,j)*sin(pi/2-phi(i,j));
rho13(i,j)= mu1s*2*cos(pi/2-phi(i,j))*sin(pi/2-phi(i,j))+ &
mu1s*sin(pi/2-phi(i,j))*cos(pi/2-phi(i,j));
rho14(i,j)=mu1s*2*cos(pi/2-phi(i,j))*cos(pi/2-phi(i,j)) &
-mu1s*sin(pi/2-phi(i,j))*sin(pi/2-phi(i,j));
rho15(i,j)=mu2s*0.75*(1+3*cos(2*(pi/2-phi(i,j)))*sin(pi/2-phi(i,j))) &
+mu2s*0.75*2*sin(2*(pi/2-phi(i,j)))*cos(pi/2-phi(i,j));
rho16(i,j)=mu2s*0.75*(1+3*cos(2*(pi/2-phi(i,j)))*cos(pi/2-phi(i,j))) &
-mu2s*0.75*2*sin(2*(pi/2-phi(i,j)))*sin(pi/2-phi(i,j));
rho17(i,j)=kappa2s*0.75*(1+3*cos(2*(pi/2-phi(i,j)))*sin(pi/2-phi(i,j)));
rho18(i,j)=kappa2s*0.75*(1+3*cos(2*(pi/2-phi(i,j)))*cos(pi/2-phi(i,j)));
rho19(i,j)=mu3s*4*leg1(i,j)*sin(pi/2-phi(i,j))+mu3s*leg2(i,j)*cos(pi/2-phi(i,j));
rho20(i,j)=mu3s*4*leg1(i,j)*cos(pi/2-phi(i,j))-mu3s*4*leg2(i,j)*sin(pi/2-phi(i,j));
rho21(i,j)=kappa3s*4*leg1(i,j)*sin(pi/2-phi(i,j))+kappa3s*(1/3)*leg2(i,j)*cos(pi/2-phi(i,j));
rho22(i,j)=kappa3s*4*leg1(i,j)*cos(pi/2-phi(i,j))-kappa3s*(1/3)*leg2(i,j)*sin(pi/2-phi(i,j));
end do
end do 

do k=1,k_current
do j=1,jmax
do i=1,imax
R1(i,j)=-(rho3(i,j)+rho5(i,j)+rho7(i,j)+rho9(i,j)+rho11(i,j));
R2(i,j)=+(rho13(i,j)+rho15(i,j)+rho17(i,j)+rho19(i,j)+rho21(i,j));
R3(i,j)=-(rho4(i,j)+rho6(i,j)+rho8(i,j)+rho10(i,j)+rho12(i,j));
R4(i,j)=+(rho14(i,j)+rho16(i,j)+rho18(i,j)+rho20(i,j)+rho22(i,j));
T(i,j)=  ((rho1(i,j)+R1(i,j)*sin(t_small*2*pi/Tcycle)+R2(i,j) &
*cos(t_small*2*pi/Tcycle))**2+ &
(rho2(i,j)+R3(i,j)*sin(t_small*2*pi/Tcycle)+ &
R4(i,j)*cos(t_small*2*pi/Tcycle))**2)**(1.0/2.0);      
end do
end do
!radius_cyl=reshape(T*sin((4*atan(1.0)*0.5)-phi),(/imax*jmax/));
!print*,'point a';
radius_cyl=T(:,1)*sin((4*atan(1.0)*0.5)-phi(:,1));
!print*,'size radius_cyl',size(radius_cyl)
z_cyl=T(:,1)*cos((4*atan(1.0)*0.5)-phi(:,1));
!print*,'size z_cyl',size(z_cyl)
!print*,'point b';
!z_cyl=reshape(T*cos(4*atan(1.0)*0.5-phi),(/imax*jmax/));
end do
end subroutine
end module mod_deformation
