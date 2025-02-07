%Intervention Analysis

load("../Ali/data_australia.mat");
time = 22; t = 1:time;


beta_n = 0.00000000337454;
epsilon_n = 0.0377;
gamma_n = 0.0566;
beta_p = 0.000040721;
epsilon_p = 0.8603;
gamma_p = 0.6437;
delta = 0.00000000055335;
lambda_n = 0.4403;
lambda_p = 0.5277;
L_n_0 = 384;
L_p_0 = 94;


% Increase in awareness of HIV precautions
close all;
alphas = [0 0.25 0.5 0.75];

outs = zeros(8, time, length(alphas));

for i = 1:length(alphas)

    
    alpha = alphas(i);

    
    par = [beta_n*(1-alpha), epsilon_n gamma_n beta_p*(1-alpha) epsilon_p gamma_p delta 0 lambda_n lambda_p];
    inital_values = [data(1,1) data(2,1) data(3,1) data(4,1) data(5,1) data(6,1) L_n_0 L_p_0];
    
    model_out = SLIR_free(par, inital_values, time);

    outs(:, :, i) = model_out;
    
end


subplot(2, 4, 1)
hold on;
for i = 1:length(alphas)
    plot(t, outs(1, :, i) + outs(7, :, i));
end
plot(t, data(1,:), 'k.')
title("S-")
xlabel('Years'), ylabel('Population')
hold off;

subplot(2, 4, 2)
hold on;
for i = 1:length(alphas)
    plot(t, outs(7, :, i));
end
% plot(t, data(2,:), 'k.')
xlabel('Years'), ylabel('Population')
title("L-")
hold off;



subplot(2, 4, 3)
hold on;
for i = 1:length(alphas)
    plot(t, outs(2, :, i));
end
plot(t, data(2,:), 'k.')
title("I-")
xlabel('Years'), ylabel('Population')
hold off;


subplot(2, 4, 4)
hold on;
for i = 1:length(alphas)
    plot(t, outs(3, :, i));
end
plot(t, data(3,:), 'k.')
title("R-")
xlabel('Years'), ylabel('Population')
hold off;

subplot(2, 4, 5)
hold on;
for i = 1:length(alphas)
    plot(t, outs(4, :, i) + outs(8, :, i));
end
plot(t, data(4,:), 'k.')
title("S+")
xlabel('Years'), ylabel('Population')
hold off;

subplot(2, 4, 6)
hold on;
for i = 1:length(alphas)
    plot(t, outs(8, :, i));
end
% plot(t, data(3,:), 'k.')
title("L+")
xlabel('Years'), ylabel('Population')
hold off;

subplot(2, 4, 7)
hold on;
for i = 1:length(alphas)
    plot(t, outs(5, :, i));
end
plot(t, data(5,:), 'k.')
title("I+")
xlabel('Years'), ylabel('Population')
hold off;

subplot(2, 4, 8)
hold on;
for i = 1:length(alphas)
    plot(t, outs(6, :, i));
end
plot(t, data(6,:), 'k.')
title("R+")
xlabel('Years'), ylabel('Population')
hold off;

legend('0', "0.25", "0.5", "0.75", "Location", [0.85, 0.85, 0.15, 0.15])
% leg = legend('show');
% set(leg, 'Title', 'Redution in Delta');
sgtitle('Intevention Simulation for Reduction in Beta')


% hold on;

% for i = 1:length(a)
%     plot(t, outs(3) );
% end
% hold off;

% savefig(".fig");


%%


subplot(1, 2, 1)
hold on;
for i = 1:length(alphas)
    % plot(t, outs(1, :, i) + outs(7, :, i));
    plot(t, outs(1,:, i) + outs(2,:, i) + outs(3,:, i))
end
plot(t, data(1,:) + data(2,:) + data(3,:), 'k.')
xlabel('Years')
ylabel('Population')
title('Total HIV-negative Population')
hold off;


subplot(1, 2, 2)
hold on;
for i = 1:length(alphas)
    % plot(t, outs(1, :, i) + outs(7, :, i));
    plot(t, outs(4,:, i) + outs(5,:, i) + outs(6,:, i))
end
plot(t, data(4,:) + data(5,:) + data(6,:), 'k.')
xlabel('Years')
ylabel('Population')
title('Total HIV-positive Population')
hold off;



% 
% % Total HIV+
% subplot(1, 2, 2)
% hold on;
% plot(t, data(4,:) + data(5,:) + data(6,:), 'k.')
% plot(t, model_out(4,:) + model_out(5,:) + model_out(6,:), 'red')
% xlabel('Years')
% ylabel('Population')
% title('Total HIV-positive Population')
% hold off;

