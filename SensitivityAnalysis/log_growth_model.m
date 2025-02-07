function [out] = log_growth_model(x0, r, N)
%The function "SIV_brain_model" returns the x(t) amd y(t) solution of the 
% %SIV brain model.

t = 0:1:15;


f = @(t,x) [x(1, :) * (r - r/N*x(1, :))];

[~, xa] = ode45(f, t, x0);


%Returning only the final value of the time-series data
out = (xa(16,1));

end



