%Parameter Estimation for SLIRS Model
function [final_size] = SLIR_sens(par, initial_values, time)
%parameter list
% [beta_n epsilon_n gamma_n beta_p epsilon_p gamma_p delta S_n_0 I_n_0 R_n_0 S_p_0 I_p_0 R_p_0]

beta_n = 10^par(1);
epsilon_n = 10^par(2);
gamma_n = 10^par(3);
beta_p = 10^par(4);
epsilon_p = 10^par(5);
gamma_p = 10^par(6);
delta = 10^par (7);
lambda_n= 10^par (9);
lambda_p= 10^par (10);

S_n_0 = initial_values(1); L_n_0 = 10^initial_values(7); I_n_0 = initial_values(2);  R_n_0 = initial_values(3);
S_p_0 = initial_values(4); L_p_0 = 10^initial_values(8); I_p_0 = initial_values(5);  R_p_0 = initial_values(6); 


% % % % % % Fixed Values (Uncomment below lines to fix any value)
% beta_n = 9.18E-09;
% epsilon_n = 0.292;
% gamma_n = 0.2691;
% beta_p = 4e-5;
% epsilon_p = 0.172;
% gamma_p = 0.0818;
% delta = 8.85E-10;
% lambda_n = 0.0664;
% lambda_p = 0.002;
% 
% L_n_0 = 4.66E+04;
% L_p_0 = 4.58E+03;  

% beta_n = 0.0000000027665;
% epsilon_n = 0.0271;
% gamma_n = 0.025;
% beta_p = 0.000048186;
% epsilon_p = 0.6038;
% gamma_p = 0.432;
% delta = 0.000000000594;
% lambda_n = 0.4403;
% lambda_p = 0.5277;

S_n_0 = S_n_0 - L_n_0;
S_p_0 = S_p_0 - L_p_0;

% x1 = S-   x2 = L-   x3 = I-     x4 = R-     
% x5 = S+   x6 = L+   x7 = I+     x8 = R+

% ( x(5)+x(6)+x(7)+x(8) ) = ( x(5)+x(6)+x(7)+x(8) )

H = @(t, x)[ -beta_n*x(1)*x(3) + epsilon_n*x(4) - delta*x(1)*( x(5)+x(6)+x(7)+x(8) );         +beta_n*x(1)*x(3) - lambda_n*x(2) - delta*x(2)*( x(5)+x(6)+x(7)+x(8) );       +lambda_n*x(2) - gamma_n*x(3);       +gamma_n*x(3) - epsilon_n*x(4) - delta*x(4)*( x(5)+x(6)+x(7)+x(8) ); 
             -beta_p*x(5)*x(7) + epsilon_p*x(8) + delta*x(1)*( x(5)+x(6)+x(7)+x(8) );         +beta_p*x(5)*x(7) - lambda_p*x(6) + delta*x(2)*( x(5)+x(6)+x(7)+x(8) );       +lambda_p*x(6) - gamma_p*x(7);       +gamma_p*x(7) - epsilon_p*x(8) + delta*x(4)*( x(5)+x(6)+x(7)+x(8) );   ];

[~, out] = ode45(H, 1:1:time, [S_n_0, L_n_0, I_n_0, R_n_0, S_p_0, L_p_0, I_p_0, R_p_0]);


Y(1,:) = (out(:,1))'; % S-
Y(2,:) = (out(:,3))'; % I-
Y(3,:) = (out(:,4))'; % R-
Y(4,:) = (out(:,5))'; % S+
Y(5,:) = (out(:,7))'; % I+
Y(6,:) = (out(:,8))'; % R+
Y(7,:) = (out(:,2))'; % L-
Y(8,:) = (out(:,6))'; % L+


final_size = Y(1, time) + Y(4, time);

end