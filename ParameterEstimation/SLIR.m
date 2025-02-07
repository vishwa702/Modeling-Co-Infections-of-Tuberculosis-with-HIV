function [Y] = SLIR(par, initial_values, time)
% 6 Compartment Model with 0-delta for I compartment
%parameter list
% [beta_n epsilon_n gamma_n beta_p epsilon_p gamma_p delta S_n_0 I_n_0 R_n_0 S_p_0 I_p_0 R_p_0]
%widths = [1e-9 1e-5;1e-9 1e-1; 1e-7 1e1; 1e-8 1e-6;1e-11 1e-9; 1e-2 1; 1e-10 1e-8; 1e-2 1];


% beta_n = par(1);
% epsilon_n = par(2);
% gamma_n = par(3);
% beta_p = par(4);
% epsilon_p = par(5);
% gamma_p = par(6);
% delta = par (7);
% lambda_n= par (9);
% lambda_p= par (10);


beta_n = 10^par(1);
epsilon_n = 10^par(2);
gamma_n = 10^par(3);
beta_p = 10^par(4);
epsilon_p = 10^par(5);
gamma_p = 10^par(6);
delta = 10^par (7);
lambda_n= 10^par (9);
lambda_p= 10^par (10);


beta_n = 3.7454e-09; %beta_n = 4.5e-9;              % DP to R and I, IP to S
epsilon_n = 3.77e-2;          % DP to I, IP to R
gamma_n =  5.66e-2;%0.0396;     % DP to R, IP to I .  
lambda_n =  4.403e-1;

beta_p = 4.0721e-05;             % DP to , IP to
epsilon_p = 8.603e-1;           % DP to S, IP to R
gamma_p = 6.437e-1;             % DP to , IP to I
lambda_p = 5.277e-1;
delta = 5.5e-10;

% 
% L_n_0 = 0.5e3;
% L_p_0 = 8.5e1;  

% beta_n = 2.5e-9;
% epsilon_n = 8.5e-8;
% gamma_n = 0.25e-1;
%gamma_n = 0.0396;
 % beta_p = 4e-5;
 % epsilon_p = 7e-1;
 % gamma_p = 6e-6;
%delta = 6e-6;
%delta = 10^-2.2227;
% beta_n = 2.5e-9;
% epsilon_n = 8.5e-8;
% gamma_n = 0.25e-1;
% %gamma_n = 0.0396;
% % % lambda_n = 1e-2;
% % 
%  beta_p = 4e-5;
%  epsilon_p = 7e-1;
%  gamma_p = 6e-6;
% lambda_p = 0.8e-2;
%delta = 0;

S_n_0 = initial_values(1); L_n_0 = 10^initial_values(7); I_n_0 = initial_values(2);  R_n_0 = initial_values(3);
S_p_0 = initial_values(4); L_p_0 = 10^initial_values(8); I_p_0 = initial_values(5); R_p_0 = initial_values(6); 

L_n_0 = 384.3134;
L_p_0 = 94.4834;  

S_n_0 = S_n_0 - L_n_0;
S_p_0 = S_p_0 - L_p_0;

% x1 = S-   x2 = L-   x3 = I-     x4 = R-     
% x5 = S+   x6 = L+   x7 = I+     x8 = R+

% ( x(5)+x(6)+x(7)+x(8) ) = ( x(5)+x(6)+x(7)+x(8) )

H = @(t, x)[ -beta_n*x(1)*x(3) + epsilon_n*x(4) - delta*x(1)*( x(5)+x(6)+x(7)+x(8) );         +beta_n*x(1)*x(3) - lambda_n*x(2) - delta*x(2)*( x(5)+x(6)+x(7)+x(8) );       +lambda_n*x(2) - gamma_n*x(3);       +gamma_n*x(3) - epsilon_n*x(4) - delta*x(4)*( x(5)+x(6)+x(7)+x(8) ); 
             -beta_p*x(5)*x(7) + epsilon_p*x(8) + delta*x(1)*( x(5)+x(6)+x(7)+x(8) );         +beta_p*x(5)*x(7) - lambda_p*x(6) + delta*x(2)*( x(5)+x(6)+x(7)+x(8) );       +lambda_p*x(6) - gamma_p*x(7);       +gamma_p*x(7) - epsilon_p*x(8) + delta*x(4)*( x(5)+x(6)+x(7)+x(8) );   ];

[~, out] = ode23(H, 1:1:time, [S_n_0, L_n_0, I_n_0, R_n_0, S_p_0, L_p_0, I_p_0, R_p_0]);


Y(1,:) = (out(:,1))'; % S-
Y(2,:) = (out(:,3))'; % I-
Y(3,:) = (out(:,4))'; % R-
Y(4,:) = (out(:,5))'; % S+
Y(5,:) = (out(:,7))'; % I+
Y(6,:) = (out(:,8))'; % R+
Y(7,:) = (out(:,2))'; % L-
Y(8,:) = (out(:,6))'; % L+

% Y = out;
end