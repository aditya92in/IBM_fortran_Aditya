clear all;
close all
radiuscyl=dlmread('r_polygon_1.txt');
zcyl=dlmread('z_polygon_1.txt');
rp2=dlmread('r_polygon_2.txt');
zp2=dlmread('z_polygon_2.txt');
%radiuscyl=dlmread('r_cyl.txt');
%zcyl=dlmread('z_cyl.txt');
z_bb=dlmread('z_c_bb.txt');
r_bb=dlmread('r_c_bb.txt');
r_grid=dlmread('r_grid.txt');
z_grid=dlmread('z_grid.txt');
co_d=dlmread('coordinates.txt');
%r_inside_general=dlmread('r_general_inside.txt');
%z_inside_general=dlmread('z_general_inside.txt');
%r_polygon2=dlmread('r_polygon2.txt');
%z_polygon2=dlmread('z_polygon2.txt');

%for i=1:6400
 %   M=csvread('deka.dat',i);
%end
%%
[m,n]=size(r_grid);
%for i=1:m
 %   plot(r_polygon2(i,:),z_polygon2(i,:))
  %  hold on
  
  
    plot(radiuscyl,zcyl)
   hold on
    plot(r_bb,z_bb,'k--','LineWidth',3)
   hold on
for i=1:n
  % line(r_grid(i)*ones(length(z_grid)),z_grid)
   %hold on
   %line(r_grid,z_grid(i)*ones(length(r_grid)))
   
   
  %  hold off
   % xlim([0 2]);
    %ylim([-1.5 1.5]);
    %pause(0.1)
end
ilow=co_d(3);
ihigh=co_d(4);
jhigh=co_d(2);
jlow=co_d(1);
dr=co_d(5);
dz=co_d(6);
for i=1:1154*962
    y_test=zp2(i,:);
    x_test=rp2(i,:);
    %[z_grid(i), z_grid(i), z_grid(i+1) ,z_grid(i+1),z_grid(i)];
    %  x_test= [r_grid(j), r_grid(j+1), r_grid(j+1), r_grid(j),r_grid(j)];
   % x_test=[z_grid(i)-dz/2, z_grid(i)-dz/2, z_grid(i+1)-dz/2, z_grid(i+1)-dz/2];
    
    plot(x_test,y_test,'b--')
    hold on
    
     %[x_test, y_test]=poly2cw(x_test,y_test);
      %   centroid1=xycentroid(x_test,y_test);
       %  [x_patch,y_patch]=polybool('intersection',x_test,y_test,radiuscyl,zcyl);
        %  phi(i,j)=polyarea(x_patch,y_patch)/(dr*dz);
         % centroid2=xycentroid(x_patch,y_patch);
    
         % if(phi(i,j)==0)
         % patch(x_test,y_test,'green');
        %phi(i,j)=(polyarea(x_patch,y_patch)/(dr*dz)) ...
%*(2*pi*(abs(r_grid(1)-centroid2(1))))/ ...
%(2*pi*(abs(r_grid(1)-centroid1(1))));
        
   %     else
         %  patch(x_patch,y_patch,'red');
          %  phi(i,j)=(polyarea(x_patch,y_patch)/(dr*dz)) ...
%*(2*pi*(abs(r_grid(1)-centroid2(1))))/ ...
%(2*pi*(abs(r_grid(1)-centroid1(1))));
            
        end 
    
   % end
%end
%vol_tot=0;
for i=ilow:ihigh
 %   for j=jlow:jhigh
  %      vol_tot=vol_tot+phi(i,j)*2*pi*abs(r_grid(1)-(r_grid(j)+dr*0.5))*dr*dz;
  %  end
end


