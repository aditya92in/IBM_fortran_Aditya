clear all;
close all
radiuscyl=dlmread('r_general.txt');
zcyl=dlmread('z_general.txt');
r_grid=dlmread('r_grid.txt');
z_grid=dlmread('z_grid.txt');
%z_bb=dlmread('z_c_bb.txt');
%r_bb=dlmread('r_c_bb.txt');
%r_inside_general=dlmread('r_general_inside.txt');
%z_inside_general=dlmread('z_general_inside.txt');
r_polygon2=dlmread('r_polygon_2.txt');
z_polygon2=dlmread('z_polygon_2.txt');

%for i=1:6400
 %   M=csvread('deka.dat',i);
%end
%%
[m,n]=size(r_grid);

  %  plot(r_polygon2(i,:),z_polygon2(i,:))
   % hold on
   
    plot(radiuscyl,zcyl,'r')
    for i=1:n
    hold on
    line(r_grid(i)*ones(size(z_grid)),z_grid)
    end
    [m,n]=size(z_grid);
    for i=1:n
    hold on
    line(r_grid,z_grid(i)*ones(size(r_grid)))
    end
    hold on
    [m,n]=size(r_polygon2);
    for i=1:m
        plot(r_polygon2(i,:),z_polygon2(i,:),'g--')
    end
    
    
    
    %plot(r_bb,z_bb)
    %hold on
    %hold off
   % xlim([0 2]);
    %ylim([-1.5 1.5]);
    %pause(0.1)
%end