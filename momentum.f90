subroutine momentum(u_old,v_old,A_u,D_u,A_v,D_v,dx,dy,imax,jmax,count_star,count_star_temp)
implicit none
integer imax,jmax,count_star,count_star_temp;
real A_u(1:imax+2,1:jmax+1,count_star),D_u(1:imax+2,1:jmax+1,count_star),u_old(1:imax+2,1:jmax+1)
real A_v(1:imax+1,1:jmax+2,count_star),D_v(1:imax+1,1:jmax+2,count_star),v_old(1:imax+1,1:jmax+2)
!real u_old(1:imax+2,1:jmax+1),v_old(1:imax+1,1:jmax+2)
real dt,dx,dy
real ,dimension (:,:),allocatable :: I_uxx,I_uyy,V_uxx,V_uyy,U_e,U_s,U_n,U_w,du_e,du_s,du_w,du_n
!re V_uxx(1:imax+2,1:jmax+1),V_uyy(1:imax+2,1:jmax+1)

!real I_uxx(1:imax+2,1:jmax+1),I_uyy(1:imax+2,1:jmax+1)
!real V_uxx(1:imax+2,1:jmax+1),V_uyy(1:imax+2,1:jmax+1)
!real U_e(1:imax+2,1:jmax+1),U_s(1:imax+2,1:jmax+1),U_w(1:imax+2,1:jmax+1)
!real U_n(1:imax+2,1:jmax+1)
!real du_e(1:imax+2,1:jmax+1),du_s(1:imax+2,1:jmax+1)
!real du_w(1:imax+2,1:jmax+1), du_n(1:imax+2,1:jmax+1)

!real I_vyy(1:imax+1,1:jmax+2),I_vxx(1:imax+1,1:jmax+2)
!real V_vyy(1:imax+1,1:jmax+2),V_vxx(1:imax+1,1:jmax+2)
!real V_e(1:imax+1,1:jmax+2),V_s(1:imax+1,1:jmax+2),V_n(1:imax+1,1:jmax+2)
!real V_w(1:imax+1,1:jmax+2)
!real dv_e(1:imax+1,1:jmax+2),dv_s(1:imax+1,1:jmax+2)
!real dv_w(1:imax+1,1:jmax+2), dv_n(1:imax+1,1:jmax+2)
real ,dimension (:,:),allocatable :: I_vyy,I_vxx,V_vyy,V_vxx,V_e,V_s,V_n,V_w,dv_e,dv_s,dv_w,dv_n
!real V_vyy(1:imax+1,1:jmax+2),V_vxx(1:imax+1,1:jmax+2)
!real pguess(1:imax+2,1:jmax+2)


!print*,'before allocate';
print*,'count_star_temp',count_star_temp;
print*,'count_star',count_star;
allocate(I_uxx(1:imax+2,1:jmax+1),I_uyy(1:imax+2,1:jmax+1),V_uxx(1:imax+2,1:jmax+1), &
V_uyy(1:imax+2,1:jmax+1),U_e(1:imax+2,1:jmax+1),U_s(1:imax+2,1:jmax+1),U_n(1:imax+2,1:jmax+1), &
U_w(1:imax+2,1:jmax+1))



allocate(du_e(1:imax+2,1:jmax+1),du_s(1:imax+2,1:jmax+1),du_w(1:imax+2,1:jmax+1), &
du_n(1:imax+2,1:jmax+1))

allocate(I_vyy(1:imax+1,1:jmax+2),I_vxx(1:imax+1,1:jmax+2),V_vyy(1:imax+1,1:jmax+2),&
V_vxx(1:imax+1,1:jmax+2),V_e(1:imax+1,1:jmax+2),V_s(1:imax+1,1:jmax+2), &
V_n(1:imax+1,1:jmax+2),V_w(1:imax+1,1:jmax+2))

allocate(dv_e(1:imax+1,1:jmax+2),dv_s(1:imax+1,1:jmax+2),dv_w(1:imax+1,1:jmax+2), &
dv_n(1:imax+1,1:jmax+2))
!print* ,'after allocate';


!print*,'I_uxx',allocated(I_uxx);
!print*,'I_uyy',allocated(I_uyy);
!print*,'V_uxx',allocated(V_uxx);
!print*,'V_uyy',allocated(V_uyy);
!print*,'U_e',allocated(U_e);
!print*,'U_s',allocated(U_s);
!print*,'U_n',allocated(U_n);
!print*,'U_w',allocated(U_w);

!print*,'du_e',allocated(du_e);
!print*,'du_s',allocated(du_s);
!print*,'du_w',allocated(du_w);
!print*,'du_n',allocated(du_n);
!print*,'I_vyy',allocated(I_vyy);
!print*,'I_vxx',allocated(I_vxx);
!print*,'V_vyy',allocated(V_vyy);
!print*,'V_vxx',allocated(V_vxx);
!print*,'V_e',allocated(V_e);
!print*,'V_s',allocated(V_s);
!print*,'V_n',allocated(V_n);
!print*,'V_w',allocated(V_w);
!print*,'dv_e',allocated(dv_e);
!print*,'dv_s',allocated(dv_s);
!print*,'dv_w',allocated(dv_w);
!print*,'dv_n',allocated(dv_n);







!I_uxx(:,:)=0.0; U_e(:,:)=0.0;
!I_uyy(:,:)=0.0; U_w(:,:)=0.0;
!V_uxx(:,:)=0.0; U_s(:,:)=0.0;
!V_uyy(:,:)=0.0; U_n(:,:)=0.0;
!I_vyy(:,:)=0.0; V_e(:,:)=0.0;
!I_vxx(:,:)=0.0; V_s(:,:)=0.0;
!V_vyy(:,:)=0.0; V_n(:,:)=0.0;
!V_vxx(:,:)=0.0; V_w(:,:)=0.0;
!du_e(:,:)=0.0; dv_e(:,:)=0.0;
!du_s(:,:)=0.0; dv_s(:,:)=0.0;
!du_w(:,:)=0.0; dv_n(:,:)=0.0;
!du_n(:,:)=0.0; dv_w(:,:)=0.0;
!!!!!!!!!!!!----------------------U---------------------------begins
du_e(1:imax+2,1:jmax)=(u_old(1:imax+2,2:jmax+1)-u_old(1:imax+2,1:jmax)) *(1/(dx));
du_w(1:imax+2,2:jmax)=du_e(1:imax+2,1:jmax);
!(uold(1:end,2:end)-uold(1:end,1:end-1))/(dx);
du_s(2:imax+2,1:jmax+1)=(u_old(2:imax+2,1:jmax+1)-u_old(1:imax+1,1:jmax+1))*(1/(dy));
du_n(1:imax+1,1:jmax+1)=du_s(2:imax+2,1:jmax+1);
U_e(1:imax+2,1:jmax+1-1)=(u_old(1:imax+2,2:jmax+1)+u_old(1:imax+2,1:jmax+1-1))*0.5;
U_w(1:imax+2,2:jmax+1)=U_e(1:imax+2,1:jmax+1-1);
!%(uold(1:end,2:end)+uold(1:end,1:end-1))*0.5;
U_s(2:imax+2,1:jmax+1)=(u_old(2:imax+2,1:jmax+1)+u_old(1:imax+1,1:jmax+1))*0.5;
U_n(1:imax+1,1:jmax+1)=U_s(2:imax+2,1:jmax+1);
!U_s(2:imax+2,1:jmax+1)=(u_old(2:imax+2,1:jmax+1)+u_old(1:imax+1,1:jmax+1))*0.5;
!U_n(1:imax+1,1:jmax+1)=U_s(2:imax+2,1:jmax+1);

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%***************
I_uxx(2:imax+1,2:jmax)=((U_e(2:imax+1,2:jmax)**2-U_w(2:imax+1,2:jmax)**2))*dy;
!axial component of x momentum
I_uyy(2:imax+1,2:jmax)=((-U_s(2:imax+1,2:jmax)*V_w(1:imax,3:jmax+1))+ &
U_n(2:imax+1,2:jmax)*V_w(2:imax+1,3:jmax+1))*dx;

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_uxx(2:imax+1,2:jmax)=(du_e(2:imax+1,2:jmax)-du_w(2:imax+1,2:jmax))*dy;
!viscous fluxes in axial direction
V_uyy(2:imax+1,2:jmax)=(du_n(2:imax+1,2:jmax)-du_s(2:imax+1,2:jmax))*dx;
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! U------------------ends

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!-V-begins-----------------------------
!V_e(1:imax+1,1:jmax+2)
!print*,'V_begins';
V_s(2:imax+1,1:jmax+2)=(v_old(2:imax+1,1:jmax+2)+v_old(1:imax+1-1,1:jmax+2))*0.5; !%U_w
V_n(1:imax+1-1,1:jmax+2)=V_s(2:imax+1,1:jmax+2);
!%(vold(1:end-1,1:end)+vold(2:end,1:end))*0.5;
V_e(1:imax+1,1:jmax+2-1)=(v_old(1:imax+1,2:jmax+2)+v_old(1:imax+1,1:jmax+2-1))*0.5;
V_w(1:imax+1,2:jmax+2)=V_e(1:imax+1,1:jmax+2-1);
!print*,'V-mom ends';
!dv_s(2:end,1:end)=(vold(2:end,1:end)-vold(1:end-1,1:end))./(dy); %U_w
!dv_n(1:end-1,1:end)=dv_s(2:end,1:end);
!%(vold(1:end-1,1:end)+vold(2:end,1:end))*0.5;
!dv_e(1:end,1:end-1)=(vold(1:end,2:end)-vold(1:end,1:end-1))./(dx);
!dv_w(1:end,2:end)=dv_e(1:end,1:end-1);
!print*,'dv_begins';
dv_s(2:imax+1,1:jmax+2)=(v_old(2:imax+1,1:jmax+2)-v_old(1:imax+1-1,1:jmax+2))*(1/(dy));! %U_w
dv_n(1:imax+1-1,1:jmax+2)=dv_s(2:imax+1,1:jmax+2);
!%(vold(1:end-1,1:end)+vold(2:end,1:end))*0.5;
dv_e(1:imax+1,1:jmax+2-1)=(v_old(1:imax+1,2:jmax+2)-v_old(1:imax+1,1:jmax+2-1))*(1/(dx));
dv_w(1:imax+1,2:jmax+2)=dv_e(1:imax+1,1:jmax+2-1);

!print*,'dv_ends';

I_vyy(2:imax+1-1,2:jmax+2-1)=((V_n(2:imax+1-1,2:jmax+2-1)**2-V_s(2:imax+1-1,2:jmax+2-1)**2))*dx;
!%axial component of v-momentum
I_vxx(2:imax+1-1,2:jmax+2-1)=(-V_w(2:imax+1-1,2:jmax+2-1)*U_s(3:imax+2-1,1:jmax+1-1) &
+V_e(2:imax+1-1,2:jmax+2-1)*U_s(3:imax+1-1,2:jmax+2))*dy;
!%transverse component of v-momentum
V_vyy(2:imax+1-1,2:jmax+2-1)=(dv_n(2:imax+1-1,2:jmax+2-1)-dv_s(2:imax+1-1,2:jmax+2-1))*dx;
!%viscous fluxes in axial direction
V_vxx(2:imax+1,2:jmax+2-1)=(dv_e(2:imax+1-1,2:jmax+2-1)-dv_w(2:imax+1-1,2:jmax+2-1))*dy;

!!!!!!!!!!!!!!!!!!!!!!!!!!-----------V-ends

!!!!!!!!!!!!!!!!!!!!-----------put it all together-------------!

print*,'put u together';
A_u(2:imax+1,2:jmax,count_star_temp)=I_uxx(2:imax+1,2:jmax)+I_uyy(2:imax+1,2:jmax);
!print*,'A_u(50,50,count_star_temp',A_u(50,50,count_star_temp);
D_u(2:imax+1,2:jmax,count_star_temp)=V_uxx(2:imax+1,2:jmax)+V_uyy(2:imax+1,2:jmax);
print*,'put v together';
A_v(2:imax,2:jmax+1,count_star_temp)=I_vyy(2:imax,2:jmax+1)+I_vxx(2:imax,2:jmax+1);
D_v(2:imax,2:jmax+1,count_star_temp)=V_vyy(2:imax,2:jmax+1)+V_vxx(2:imax,2:jmax+1);

print*,'end of v together';
!%transverse component of y momentum

!%viscous fluxes in transverse direction

!%%%%%%%%%%%%%%%%%%%%
!%**************-V star%

!visoucs fluxes in transverse direction
!v_star(2:imax,2:jmax+1)=(vold(2:imax,2:jmax+1)-dt/(rho*dx*dy)*(I_vyy(2:imax,2:jmax+1) &
!+I_vxx(2:imax,2:jmax+1))+dt/(dx*dy)*(V_vyy(2:imax,2:jmax+1)+V_vxx(2:imax,2:jmax+1)));
!end

!print*,'I_uxx',size(I_uxx,dim=1),size(I_uxx,dim=2)
!print*,'I_uyy',size(I_uyy,dim=1),size(I_uyy,dim=2)
!print*,'V_uxx',size(V_uxx,dim=1),size(V_uxx,dim=2)
!print*,'V_uyy',size(V_uyy,dim=1),size(V_uyy,dim=2)
!print*,'U_e',size(U_e,dim=1),size(U_e,dim=2)
!print*,'U_s',size(U_s,dim=1),size(U_s,dim=2)

!print*,'deallocate'

deallocate(I_uxx,I_uyy,V_uxx, &
V_uyy,U_e,U_s,U_n, &
U_w)

deallocate(du_e,du_s,du_w, &
du_n)

deallocate(I_vyy,I_vxx,V_vyy,&
V_vxx,V_e,V_s, &
V_n,V_w)

deallocate(dv_e,dv_s,dv_w, &
dv_n)

!print*,'end of deallocate'



!print*,'I_uxx',allocated(I_uxx);
!print*,'I_uyy',allocated(I_uyy);
!print*,'V_uxx',allocated(V_uxx);
!print*,'V_uyy',allocated(V_uyy);
!print*,'U_e',allocated(U_e);
!print*,'U_s',allocated(U_s);
!print*,'U_n',allocated(U_n);
!print*,'U_w',allocated(U_w);

!print*,'du_e',allocated(du_e);
!print*,'du_s',allocated(du_s);
!print*,'du_w',allocated(du_w);
!print*,'du_n',allocated(du_n);
!print*,'I_vyy',allocated(I_vyy);
!print*,'I_vxx',allocated(I_vxx);
!print*,'V_vyy',allocated(V_vyy);
!print*,'V_vxx',allocated(V_vxx);
!print*,'V_e',allocated(V_e);
!print*,'V_s',allocated(V_s);
!print*,'V_n',allocated(V_n);
!print*,'V_w',allocated(V_w);
!print*,'dv_e',allocated(dv_e);
!print*,'dv_s',allocated(dv_s);
!print*,'dv_w',allocated(dv_w);
!print*,'dv_n',allocated(dv_n);
end subroutine


!psuedo code

!imax1=number of columns in v
!jmax1=number of rows in v
!imax=number of columns in u
!jmax=number of columns in u


!note 
!               x   
!                   
!           x . x . x     stencil for uold

!               x
  
