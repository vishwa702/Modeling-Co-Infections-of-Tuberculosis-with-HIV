function [loglike,m_out] = loglike_SLIR_model(theta,data)
% Parameters In Model:
%   beta_n, epsilon_n, gamma_n, lambda_n, 
%   beta_p, epsilon_p, gamma_p, lambda_p,
% %   delta;


%beta_n = theta(1);
%epsilon_n = theta(2);
%gamma_n = theta(3);
%beta_p = theta(4);
%epsilon_p = theta(5);
%gamma_p = theta(6);
%delta = theta (7);
%p = 10^theta (8);


initial_values = [19005317 1090 1037 10318 110 91 theta(11) theta(12)];
time = 22;
m_out = SLIR(theta, initial_values,time);
%fprintf('FirstRow: %d %d %d %d', m_out(1, 1), m_out(1, 5), m_out(1, 10), m_out(1, 15        ))

check_size = size(m_out);

% fprintf('\nCheck_size Rows = %d', check_size(1))
% fprintf('\nCheck_size Cols = %d', check_size(2))

%This if statement ensures that the resulting numerical solution vector is 
%the correct length. If the numerical solution vector has the incorrect 
%length, the parameter vector that produced this solution vector is 
%considered a bad guess and this parameter vector is assigned the log 
%likelihood value of -Inf.
if(check_size(2) == length(1:1:22))
    % fprintf('\nP = %d', p);



    t1 = 6.75e-6;
    t2 = 1.86e-5;
    t3=1.80e-5;
    t4=1.52e-6;
    t5=1.43e-3;
    t6=2.38e-3;

    % %Gaussian 
    loglike_1 = (22/2 * log(t1)) - 1/2*t1 * sum( (double(data(1,:))-(m_out(1,:) + m_out(7,:))).^2 );
    loglike_2 = (22/2 * log(t2)) - 1/2*t2 * sum( (double(data(2,:))-m_out(2,:)).^2 );
    loglike_3 = (22/2 * log(t3)) - 1/2*t3 * sum( (double(data(3,:))-m_out(3,:)).^2 );
    loglike_4 =(22/2 * log(t4)) - 1/2*t4 * sum( (double(data(4,:))-(m_out(4,:) + m_out(8,:))).^2 );
    loglike_5 = (22/2 * log(t5)) - 1/2*t5 * sum( (double(data(5,:))-m_out(5,:)).^2 );
    loglike_6 = (22/2 * log(t6)) - 1/2*t6 * sum( (double(data(6,:))-m_out(6,:)).^2 );


    loglike = loglike_1 + loglike_2 + loglike_3 + loglike_4 + loglike_5 + loglike_6;
    %fprintf('\nLogLikes: %d | %d | %d | %d | %d | %d', loglike_1, loglike_2, loglike_3, loglike_4, loglike_5, loglike_6)

else

    loglike = -Inf;

end

end