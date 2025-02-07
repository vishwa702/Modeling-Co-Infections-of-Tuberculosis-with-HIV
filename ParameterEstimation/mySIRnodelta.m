function [Y] = mySIRnodelta(par, initial_values, time)
% 6 Compartment Model with 0-delta for I compartment
%parameter list
% [beta_n epsilon_n gamma_n beta_p epsilon_p gamma_p delta S_n_0 I_n_0 R_n_0 S_p_0 I_p_0 R_p_0]
%widths = [1e-9 1e-5;1e-9 1e-1; 1e-7 1e1; 1e-8 1e-6;1e-11 1e-9; 1e-2 1; 1e-10 1e-8; 1e-2 1];

beta_n = 10^par(1);
epsilon_n = 10^par(2);
gamma_n = 10^par(3);
beta_p = 10^par(4);
epsilon_p = 10^par(5);
gamma_p = 10^par(6);
delta = 10^par (7);

beta_n = 1e-20;
epsilon_n = 1e-9;
% gamma_n = 0.0396;
% lambda_n = 1e-2;

 beta_p = 1.8e-8;
 epsilon_p = 1.5e-10;
 gamma_p = 0.8e-1;
% lambda_p = 0.8e-2;
delta = 0;

S_n_0 = initial_values(1);  I_n_0 = initial_values(2);  R_n_0 = initial_values(3);
S_p_0 = initial_values(4); I_p_0 = initial_values(5); R_p_0 = initial_values(6);


% x1 = S-   x2 = I-     x3 = R-     
% x4 = S+   x5 = I+     x6 = R+

H = @(t, x)[ -beta_n*x(1)*x(2) + epsilon_n*x(3) - delta*x(1);         +beta_n*x(1)*x(2) - gamma_n*x(2) - delta*x(2);       +gamma_n*x(2) - epsilon_n*x(3) - delta*x(3); 
             -beta_p*x(4)*x(5) + epsilon_p*x(6) + delta*x(1);         +beta_p*x(4)*x(5) - gamma_p*x(5) + delta*x(2);       +gamma_p*x(5) - epsilon_p*x(6) + delta*x(3);];

[~, out] = ode45(H, 1:1:time, [S_n_0, I_n_0, R_n_0, S_p_0, I_p_0, R_p_0]);

% plot(1:1:time,out);
% legend('1','2','3','4','5','6');

% out(:,1) = S_n_0;

Y(1,:) = (out(:,1))';
Y(2,:) = (out(:,2))';
Y(3,:) = (out(:,3))';
Y(4,:) = (out(:,4))';
Y(5,:) = (out(:,5))';
Y(6,:) = (out(:,6))';


% Y = out;
end
