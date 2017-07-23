subroutine solve_poisson(pcorr,imax,jmax,dx,dy,res_p,rhs)
integer imax,jmax,m_p,n_p
real pcorr(1:imax+2,1:jmax+2),rhs(1:imax+2,1:jmax+2)
real ,dimension(:,:),allocatable::piter(:,:)
real dx,dy,res_p,count_star,dx2,dy2;
allocate (piter(1:imax+2,1:jmax+2))
dx2=1/dx**2;
dy2=1/dy**2;

m_p=size(pcorr,dim=1);
n_p=size(pcorr,dim=2);

  res_p=1;
 piter=pcorr;
!% %pnew=piter;
 count_star=1;
! tic
!omega=1.5;
 do while(res_p>1e-6.and.count_star<50) 
!boundary conditions
     piter(:,1)=piter(:,2); !%dp/dx=0; west boundary
     piter(:,n_p)=-piter(:,n_p-1); !%p=0 east boundary
     piter(1,:)=piter(2,:);    !%dp/dy=0 !; south boundary
     piter(m_p,:)=piter(m_p-1,:);  !dp/dy=0 ;north boundary
    
  do i=2,imax+1
      do j=2,jmax+1
!          %piter(i,j)=piter(i,j)*(1-omega)+omega*0.5*((piter(i,j-1)+pnew(i,j+1)).*power(dy2,-1)+(pnew(i+1,j)+piter(i-1,j)).*power(dx2,-1)-rhs(i,j)*power(dx2,-1)*power(dy2,-1))./(power(dx2,-1)+power(dy2,-1));
          piter(i,j)=0.5*(((piter(i,j-1)+pcorr(i,j+1))*(1/dy2)+ &
(pcorr(i+1,j)+piter(i-1,j))*(1/(dx2))-rhs(i,j)*(1/dx2)*(1/dy2)))/((1/dx2)+(1/dy2));
      end do
   end do
 
!piter=adi_fourier(pnew',rhs',dx,dy);
 res_p=norm2(reshape(pcorr(2:m_p-1,2:n_p-1),(/1,imax*(jmax)/)) &
-reshape(piter(2:m_p-1,2:n_p-1),(/1,imax*(jmax)/)));
!res_y=norm2(reshape(v_new(2:m_v-1,2:n_v-1),(/1,(imax-1)*jmax/)) &
!-reshape(v_old(2:m_u-1,2:n_u-1),(/1,(imax-1)*jmax /)));
 !err_pos=norm(pnew(2:end-1,2:end-1)-piter(2:end-1,2:end-1));

! %plot(count,err,'.','MarkerSize',10);
! %hold all
 count_star=count_star+1;
!pause(0.1)
 pcorr=piter;
end  do
    pcorr(:,1)=pcorr(:,2); !dp/dx=0;
    pcorr(:,n_p)=-pcorr(:,n_p-1); !P=0
    pcorr(1,:)=pcorr(2,:);!dp/dy=0; south boundary
    pcorr(m_p,:)=pcorr(m_p-1,:);!dp/dy=0 north boundary
    !pcorr=pnew;
    !disp(['time spent in poissson=',num2str(toc)]);
    deallocate (piter)
end subroutine
