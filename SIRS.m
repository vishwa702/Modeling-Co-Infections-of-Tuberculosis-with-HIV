function [Y] = SIRS(par, pars,time)
% 6 Compartment Model with 0-delta for I compartment
%parameter list
% [beta_n epsilon_n gamma_n beta_p epsilon_p gamma_p delta S_n_0 I_n_0 R_n_0 S_p_0 I_p_0 R_p_0]

beta_n = par(1);
epsilon_n = par(2);
gamma_n = par(3);
beta_p = par(4);
epsilon_p = par(5);
gamma_p = par(6);
delta = par (7);

S_n_0 = pars(1);  I_n_0 = pars(2);  R_n_0 = pars(3);
S_p_0 = pars(4); I_p_0 = pars(5); R_p_0 = pars(6);


% x1 = S-   x2 = I-     x3 = R-     
% x4 = S+   x5 = I+     x6 = R+

H = @(t, x)[ -beta_n*x(1)*x(2) + epsilon_n*x(3) - delta*x(1);         +beta_n*x(1)*x(2) - gamma_n*x(2) - delta*x(2);       +gamma_n*x(2) - epsilon_n*x(3) - delta*x(3); 
             -beta_p*x(4)*x(5) + epsilon_p*x(6) + delta*x(1);         +beta_p*x(4)*x(5) - gamma_p*x(5) + delta*x(2);       +gamma_p*x(5) - epsilon_p*x(6) + delta*x(3);];

[~, Y] = ode45(H, 0:1:time, [S_n_0, I_n_0, R_n_0, S_p_0, I_p_0, R_p_0]);

end
%newest version