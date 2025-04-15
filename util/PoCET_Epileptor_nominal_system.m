function dXdt = PoCET_Epileptor_nominal_system(t,X,PAR)

 x_1 = X(1);
 x_2 = X(2);
 x_3 = X(3);
 x_4 = X(4);
 x_5 = X(5);
 u = X(6);
 
if x_1 < 0
    f_1 = x_1^3-3*x_1^2;
else
    f_1 = (x_4-0.6*(x_3-4)^2-PAR.m)*x_1;
end

if x_4 < -0.25
    f_2 = 0;
else
    f_2 = 6*(x_4+0.25);
end
 
 ddt_x_1 = x_2 - f_1 - x_3 + PAR.I_1;
 ddt_x_2 = PAR.r_2 - 5*x_1^2 - x_2;
 ddt_x_3 = 1/PAR.tau_0*(4*(x_1 - PAR.r_1) - x_3);
 ddt_x_4 = -x_5 + x_4 -x_4^3 + PAR.I_2 + 2*u - 0.3*(x_3 - 3.5);
 ddt_x_5 = 1/PAR.tau_2*(-x_5 + f_2);
 ddt_u = -PAR.gamma*(u-0.1*x_1);

 dXdt = [ddt_x_1; ddt_x_2; ddt_x_3; ddt_x_4; ddt_x_5; ddt_u];
end