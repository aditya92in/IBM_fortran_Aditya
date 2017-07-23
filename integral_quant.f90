subroutine integral_quant(F_x,F_y,d,U_w,rho,C_d,C_l)

real F_x,F_y,d,U_w,rho,C_d,C_l

!function [C_d,C_l]=integral_quant(F_x,F_y,d,U_w,rho)

 C_d=F_x/(rho*U_w*U_w*0.5*d); 
 C_l=F_y/(rho*U_w*U_w*0.5*d);
end subroutine 
