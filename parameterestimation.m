%years
t = 1:1:22; time = 22;


%Load the given data
load ('Matrices/AusData.mat')
%data = table2array(data);
data = double(data);
fprintf("done")
%% Parameters to Fit

% Parameters In Model:
%   beta_n, epsilon_n, gamma_n,
%   beta_p, epsilon_p, gamma_p,
% %   delta;
% %   p;
%%lambda_n, lambda_p
%%L_n L_p
%10^-9 10^-5 og beta_n values
widths = [log10(1e-10) log10(1e-8);log10(1e-4) log10(1e0); log10(1e-4) log10(1e1); log10(1e-7) log10(1e-3);log10(1e-3) log10(1e1); log10(1e-3) log10(1e1); log10(1e-12) log10(1e-8); log10(1e-3) log10(1e1) ; log10(1e-3) log10(1e1); log10(1e-3) log10(1e1); log(1e0) log(1e5); log(1e0) log(1e5)];



%We suppose that the parameter beta lies between 1e-9 and 1e-4. Since this 
%range spans several orders of magnitude in terms of base 10, it is 
%beneficial to use the loguniform distribution as the prior distribution 
%for the parameter beta. We do this by fitting log base 10 of beta with a 
%uniform distribution prior and then exponentiating the resulting samples 
%of log base 10 of beta by 10.
%
%It is useful to use the loguniform distribution as the prior distribution 
%for parameters that take positive values and have an uncertainty that 
%spans several orders of magnitude. This technique can help speed up
%parameter estimation.

%Using the loguniform distribution as the prior distribution for a 
%parameter is discussed in the following book chapter:
%Brewer BJ (2018) Bayesian inference and computation: a beginner's guide. 
%In: Ramos AA & Arregui I (Eds.). Bayesian Astrophysics (pp. 221-232). 
%Cambridge University Press.

%Prior for x0, y0, log10_beta, and p are the following uniform 
%distributions
prior1 = @(log10_beta_n) unifpdf(log10_beta_n, widths(1,1), widths(1,2));
%prior1 = @(beta_n) unifpdf(beta_n, widths(1,1), widths(1,2));
prior2 = @(log10_epsilon_n) unifpdf(log10_epsilon_n, widths(2,1), widths(2,2));
prior3 = @(log10_gamma_n) unifpdf(log10_gamma_n, widths(3,1), widths(3,2));
%prior4 = @(lambda_n) unifpdf(lambda_n, widths(4,1), widths(4,2));
prior4 = @(log10_beta_p) unifpdf(log10_beta_p, widths(4,1), widths(4,2));
prior5 = @(log10_epsilon_p) unifpdf(log10_epsilon_p, widths(5,1), widths(5,2));
prior6 = @(log10_gamma_p) unifpdf(log10_gamma_p, widths(6,1), widths(6,2));
%prior8 = @(lambda_p) unifpdf(lambda_p, widths(8,1), widths(8,2));
prior7 = @(log10_delta) unifpdf(log10_delta, widths(7,1), widths(7,2));
prior8 = @(log10_p) unifpdf(log10_p, widths(8,1), widths(8,2));

%lambda_n lambda_p L_n L_p
prior9 = @(log10_lambda_n) unifpdf(log10_lambda_n, widths(9,1), widths(9,2));
prior10 = @(log10_lambda_p) unifpdf(log10_lambda_p, widths(10,1), widths(10,2));
prior11 = @(log10_L_n) unifpdf(log10_L_n, widths(11,1), widths(11,2));
prior12 = @(log10_L_p) unifpdf(log10_L_p, widths(12,1), widths(12,2));

%'par_names' is a string array holding the names of the parameters
par_names = ["beta_n", "epsilon_n", "gamma_n", "beta_p", "epsilon_p", "gamma_p", "delta" ,"p","lambda_n", "lambda_p", "L_n" ,"L_p"];
fprintf("\nParameters to Fit done")
%% Parameter Fit Setup
% clc;

%Given the vector of parameters theta = [x0, y0, log10_beta, p], logl 
%returns the log likelihood distribution
logl = @(theta) (loglike_SLIR_model(theta, data));

%'logp' is the natural logorithm of the prior distribution for the 
%parameters. The prior distribution is a product of uniform
%distributions, the resulting joint distribution does integrate to one.
%p taken out: + log(prior8(theta(8)))
%logp = @(theta) (log(prior1(theta(1))) + log(prior2(theta(2))) + log(prior3(theta(3))) + log(prior4(theta(4))) + log(prior5(theta(5))) + log(prior6(theta(6))) + log(prior7(theta(7))) + log(prior9(theta(9))) + log(prior10(theta(10))) + log(prior11(theta(11))) + log(prior12(theta(12))));
%logp = @(theta) (log(prior3(theta(3)))+log(prior6(theta(6))));
logp = @(theta) (log(prior2(theta(2))) + log(prior4(theta(4))) + log(prior9(theta(9))) + log(prior10(theta(10))) + log(prior11(theta(11))) + log(prior12(theta(12))));

theta0 = zeros(12, 24); % Good Rule of Thumb: chains = 2 times parameters that you will fit
% 
%Generate the vector of initial conditions

k = 0;
j = 1;

while(j < (24+1))

while(k == 0 || isnan(k) || isinf(k))
    %widths = [1e-9 1e-5;1e-9 1e-1; 1e-7 1e1; 1e-8 1e-6;1e-11 1e-9; 1e-2 1; 1e-10 1e-8; 1e-2 1];
    % initial(4) <= 0.001
    initial = [widths(1,1)+(widths(1,2)-widths(1,1))*rand widths(2,1)+(widths(2,2)-widths(2,1))*rand widths(3,1)+(widths(3,2)-widths(3,1))*rand widths(4,1)+(widths(4,2)-widths(4,1))*rand widths(5,1)+(widths(5,2)-widths(5,1))*rand widths(6,1)+(widths(6,2)-widths(6,1))*rand widths(7,1)+(widths(7,2)-widths(7,1))*rand widths(8,1)+(widths(8,2)-widths(8,1))*rand widths(9,1)+(widths(9,2)-widths(9,1))*rand widths(10,1)+(widths(10,2)-widths(10,1))*rand widths(11,1)+(widths(11,2)-widths(11,1))*rand widths(12,1)+(widths(12,2)-widths(12,1))*rand];       
    %initial = [0.1 0.0327 9.5447 0.1 0.1 0.0545 0 0.5524];

    [ll,m_out] = logl(initial);
    lp = logp(initial);
    k = ll + lp;
    % k = logl(initial)
    fprintf("\nLogP: %d", lp );
    fprintf("\nLogL: %d", ll );
    fprintf("\nK: %d", k);
    fprintf("\nJ: %d", j);

end

theta0(1,j) = initial(1);
theta0(2,j) = initial(2);
theta0(3,j) = initial(3);
theta0(4,j) = initial(4);
theta0(5,j) = initial(5);
theta0(6,j) = initial(6);
theta0(7,j) = initial(7);
theta0(8,j) = initial(8);
theta0(9,j) = initial(9);
theta0(10,j) = initial(10);
theta0(11,j) = initial(11);
theta0(12,j) = initial(12);

j = j + 1;
k = 0;

display(j);

end
disp(theta0);
save theta0.mat theta0

%'logfuns' holds one function for the log of the priors and another 
%function for the log of the likelihood function
logfuns = {@(theta)logp(theta) @(theta)logl(theta)};
fprintf("\nParameter Fit Setup done")
%% Affine Fitting Algorithm

%'gwmcmc' runs the MCMC sampler with affine invariance
[models, logP] = gwmcmc(theta0, logfuns, 2000000);

%'models' holds the proposed guesses for the parameter values
save parameter_estimation_models.mat models

%'logP' holds the corresponding value of the log prior and the value of the
%log likelihood for each guess of the parameter values
save parameter_estimation_logP.mat logP
fprintf("\nAffine Fitting Algorithm done")
%% Assess results

%Let 'l' hold the values of the log prior and log likelihood
l = logP;

%Set each of the eight chains to have a burn-in of 2500 samples
l(:,:,1:2000)=[];

%Extract the values of the log prior
lprior = l(1,:,:);

%Collapse the samples of all chains into a single larger sample
lprior = lprior(:,:);

%Extract the values of the log likelihood
llik = l(2,:,:);

%Collapse the samples of all chains into a single larger sample
llik = llik(:,:);

%Add the log prior and log likelihood together to get the log unnormalized 
%posterior samples
l_un_post_samples = lprior + llik;
%l_un_post_samples = l_un_post_samples/max(l_un_post_samples);

% l_un_post_samples = l_un_post_samples/1000; %SCALING
%Let 'm' hold the values of the proposed guesses for the parameter values
m = models;

%Set each of the eight chains to have a burn-in of 2500 samples
m(:,:,1:2000)=[];

%Collapse the samples of all chains into a single larger sample
m = m(:,:);

theta_samples = m;

% Parameters In Model:
%   beta_n, epsilon_n, gamma_n,
%   beta_p, epsilon_p, gamma_p,
% %   delta;
% %   p;

%'k1' contains all the proposed values for parameter 'beta_n'
k1 = m(1,:);

%'k2' contains all the proposed values for parameter 'epsilon_n'
k2 = m(2,:);

%'k3' contains all the proposed values for parameter 'gamma_n'
k3 = m(3,:);

%'k4' contains all the proposed values for parameter 'beta_p'
k4 = m(4,:);

%'k4' contains all the proposed values for parameter 'epsilon_p'
k5 = m(5,:);

%'k4' contains all the proposed values for parameter 'gamma_p'
k6 = m(6,:);

%'k4' contains all the proposed values for parameter 'delta'
k7 = m(7,:);

%'k4' contains all the proposed values for parameter 'p'
k8 = m(8,:);
%'k4' contains all the proposed values for parameter 'lambda_n'
k9 = m(9,:);

%'k4' contains all the proposed values for parameter 'lambda_p'
k10 = m(10,:);

%'k4' contains all the proposed values for parameter 'L_n'
k11 = m(11,:);

%'k4' contains all the proposed values for parameter 'L_p'
k12 = m(12,:);

%parameter values at the maximum log unnormalized posterior value
[M, I] = max(l_un_post_samples);
fprintf("\nAssess results done")
%% Parameter Posterior Plots

%Beta_n

close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 1)
savefig("post_surf_beta_n.fig");

[l_HPD1, u_HPD1] = HPD(k1, 0.05);
k1_atmax=theta_samples(1,I);
10^k1_atmax;
k1_median=median(k1);

% 95% CI: (2.2917e+06, 2.3019e+06)
% Parameter at Maximum log unnormalized posterior value: 2.2967e+06
% Parameter at median: 2.2968e+06

% Epsilon_n
close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 2)
savefig("post_surf_epsilon_n.fig");

[l_HPD2, u_HPD2] = HPD(k2, 0.05);
k2_atmax=theta_samples(2,I);
10^k2_atmax;
k2_median=median(k2);

% 95% CI: (4.6069,19.1019)
% Parameter at Maximum log unnormalized posterior value: 8.8646
% Parameter at median: 10.7186

%gamma_n
close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 3)
savefig("post_surf_gamma_n.fig");

[l_HPD3, u_HPD3] = HPD(k3, 0.05);
k3_atmax=theta_samples(3,I);
10^k3_atmax;
k3_median=median(k3);
% Optimal log value = -2.60939 using elbow method while keeping only delta and
% gamma_n free. Actual value = 10^-2.60939 = 0.0025
% With fixed loglikelihood for other compartments, we get a unimodal curve
% posterior for gamma_n. It peaks at -1.4022. i.e. 0.0396


% 95% CI: (4.2836e-06,7.2034e-06)
% Parameter at Maximum log unnormalized posterior value: 6.0419e-06
% Parameter at median: 5.6579e-06



% beta_p
close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 4)
savefig("post_surf_beta_p.fig");

[l_HPD4, u_HPD4] = HPD(k4, 0.05);
k4_atmax=theta_samples(4,I);
10^k4_atmax;
k4_median=median(k4);



% epsilon_p
close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 5)
savefig("post_surf_epsilon_p.fig");

[l_HPD5, u_HPD5] = HPD(k5, 0.05);
k5_atmax=theta_samples(5,I);
10^k5_atmax;
k5_median=median(k5);



% gamma_p
close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 6)
savefig("post_surf_gamma_p.fig");

[l_HPD6, u_HPD6] = HPD(k6, 0.05);
k6_atmax=theta_samples(6,I);
10^k6_atmax;
k6_median=median(k6);


% delta
close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 7)
savefig("post_surf_delta.fig");

[l_HPD7, u_HPD7] = HPD(k7, 0.05);
k7_atmax=theta_samples(7,I);
k7_median=median(k7);
10^k7_atmax;
% delta atmax is outside > u_HPD. Using median log value = -9.0649 i.e. 8.6119e-10 

% lambda_n
close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 9)
savefig("post_surf_lambda_n.fig");

[l_HPD9, u_HPD9] = HPD(k9, 0.05);
k9_atmax=theta_samples(9,I);
10^k9_atmax;
k9_median=median(k9);

% lambda_p
close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 10)
savefig("post_surf_lambda_p.fig");

[l_HPD10, u_HPD10] = HPD(k10, 0.05);
k10_atmax=theta_samples(10,I);
10^k10_atmax;
k10_median=median(k10);

% L_n
close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 11)
savefig("post_surf_L_n.fig");

[l_HPD11, u_HPD11] = HPD(k11, 0.05);
k11_atmax=theta_samples(11,I);
10^k11_atmax;
k11_median=median(k11);

% L_p
close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 12)
savefig("post_surf_L_p.fig");

[l_HPD12, u_HPD12] = HPD(k12, 0.05);
k12_atmax=theta_samples(12,I);
k12_median=median(k12);



close all;
% post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 8)
% savefig("post_surf_p.fig");
% 
% [l_HPD, u_HPD] = HPD(k8, 0.05);
% k8_atmax=theta_samples(8,I);
% k8_median=median(k8);

% 95% CI: (0.0070,0.0156)
% Parameter at Maximum log unnormalized posterior value: 0.0113
% Parameter at median: 0.0110
fprintf("\n Parameter Posterior Plotsdone")
%% Graphs (Parameters at max log posterior value)
close all;
beta_n = 1e-11;
epsilon_n = 1e-8;
gamma_n = 0.0396;
% lambda_n = NA;

beta_p = 1.8e-8;
epsilon_p = 1.5e-10;
gamma_p = 0.8e-1;
% lambda_p = NA;
delta = 0e-12;
par = [beta_n, epsilon_n gamma_n beta_p epsilon_p gamma_p delta];
inital_values = [19005317 1090 1037 10318 110 91];

model_out = mySIR(par, inital_values, time);

% subplot(2,1,1)
close all;
hold on
set(0,'defaultLineLineWidth',1.5);   
set(0,'defaultLineMarkerSize',9);
set(gca, 'FontSize', 10, 'LineWidth', 1);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

plot(t, data(3,:), 'ko')
plot(t, model_out(3,:))

xlabel('Years'), ylabel('S-')

subplot(2,1,2)

set(0,'defaultLineLineWidth',1.5);   
set(0,'defaultLineMarkerSize',9);
set(gca, 'FontSize', 10, 'LineWidth', 1);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

plot(t, data(6,:), 'ko')
plot(t, model_out(6,:))

xlabel('Years'), ylabel('R+')

hold off

savefig("max_post_solution_SIR.fig");
close all
fprintf("\ndone")
%% Posterior Predictive distribution Loop

%The following 'for loop' estimates the posterior predictive distribution
%for the model solution x

t1 = 6.75e-6;
t2 = 1.86e-5;
t3=1.80e-5;
t4=1.52e-6;
t5=1.43e-3;
t6=2.38e-3;

variance = [1/t1 1/t2 1/t3 1/t4 1/t5 1/t6];
Sigma = diag(variance);

x_predict = zeros(length(data),length(l_un_post_samples));
y_predict = zeros(length(data),length(l_un_post_samples));

discrepancy = zeros(1,length(l_un_post_samples));

discrepancy_pred = zeros(1,length(l_un_post_samples));

ind_D_pred_exceeds = zeros(1,length(l_un_post_samples));

h = waitbar(0,'Initialize...');
for i = 1:length(l_un_post_samples)
    %fprintf("\ni = %d",i);
    % L_n_0 = 384.3134;
    % L_p_0 = 94.4834;  
    
    theta = theta_samples(:,i);
    initial_values = [19005317 1090 1037 10318 110 91 theta(11) theta(12)];
    
    out = SLIR(theta, initial_values,time);
    out = out(1:6,:);
    % varience = ..   We need the variance here to use in discrepancy
    
    discrepancy(1,i) = sum(((data(1,:) - out(1,:)).^2)./(out(1,:)/(variance(1)))) + sum(((data(2,:) - out(2,:)).^2)./(out(2,:)/(variance(2))));
    %+ sum(((data(3,:) - out(3,:)).^2)./(out(3,:)/(variance(3)))) + sum(((data(4,:) - out(4,:)).^2)./(out(4,:)/(variance(4)))) + sum(((data(5,:) - out(5,:)).^2)./(out(5,:)/(variance(5)))) + sum(((data(6,:) - out(6,:)).^2)./(out(6,:)/(variance(6))));
    
    %Sigma=... we need the Sigma matrix here to use in  mvnrnd

    x_predict(:,i) = mvnrnd(out(1,:)', Sigma(1), length(data(1,:)));
    y_predict(:,i) = mvnrnd(out(2,:)', Sigma(2), length(data(2,:)));

    
    discrepancy_pred(1,i) = sum((((x_predict(:,i))' - out(1,:)).^2)./(out(1,:)/(variance(1)))) + sum((((y_predict(:,i))' - out(2,:)).^2)./(out(2,:)/(variance(2))));
    
    if(discrepancy_pred(1,i) > discrepancy(1,i))
        
        ind_D_pred_exceeds(1,i) = 1;
        
    else
        
        ind_D_pred_exceeds(1,i) = 0;
        
    end
    
   waitbar(i/length(l_un_post_samples),h,sprintf('%d%%',(i/length(l_un_post_samples))*100))
   
end
close(h)

save x_predict_SIV.mat x_predict
save y_predict_SIV.mat y_predict
save discrepancy_SIV.mat discrepancy
save discrepancy_pred_SIV.mat discrepancy_pred
save ind_D_pred_exceeds_SIV.mat ind_D_pred_exceeds

%%
%Bayesian p-value (using the chi squared realized discrepancy function)
%Given theta, the discrepancy function chosen has an approximate chi 
%squared distribution with n degrees of freedom.

%(Note: if the data and predicted data had came from a normal distribution 
%with constant variance, then the discrepancy function chosen would have a 
%chi squared distribution with n degrees of freedom. In this example, the 
%data and predicted data came from a negative binomial distribution; so, 
%the discrepancy function chosen has an approximate chi squared 
%distribution with n degrees of freedom.)

p_value = sum(ind_D_pred_exceeds)/length(ind_D_pred_exceeds);

%This Bayesian p-value measures how extreme is the realized value of the
%chi squared discrepancy among all possible values that could have been
%realized under this SIV brain model with the same value of theta that 
%generates the current data.

%The Bayesian p-value (using the chi squared realized discrepancy function)
%is 0.6955, which indicates that there is no evidence against the null 
%hypothesis that the model predictions fit the data.

%Note: If this Bayesian p-value was close to 0 or 1, then there would be 
%evidence against the null hypothesis that the model predictions fit the 
%data. Extreme tail-area probabilities, less than 0.01 or more than 0.99, 
%indicate a major failure of the model predictions to fit the data.

%%

%Predictive mean for the Susceptible brain macrophage population
mean_x_predict = (1/(length(l_un_post_samples)))*sum(x_predict');

%Generate 95% prediction interval
lower_predict_x = zeros(1, length(data));

upper_predict_x = zeros(1, length(data));

h = waitbar(0,'Initialize...');
for i = 1:length(data)
    
    lower_predict_x(1, i) = prctile(x_predict(i,:),2.5);
    upper_predict_x(1, i) = prctile(x_predict(i,:),97.5);
    
    waitbar(i/length(data),h,sprintf('%d%%',(i/length(data))*100))
end
close(h)

%Predictive mean for the SIV-infected brain macrophage population
mean_y_predict = (1/(length(l_un_post_samples)))*sum(y_predict');

%Generate 95% prediction interval
lower_predict_y = zeros(1, length(data));

upper_predict_y = zeros(1, length(data));

h = waitbar(0,'Initialize...');
for i = 1:length(data)
    
    lower_predict_y(1, i) = prctile(y_predict(i,:),2.5);
    upper_predict_y(1, i) = prctile(y_predict(i,:),97.5);
    
    waitbar(i/length(data),h,sprintf('%d%%',(i/length(data))*100))
end;
close(h)

%%

model_true = SIV_brain_model(2.3e6, 8, 5.5e-6);

%Plot the model fitting and prediction
hold on
set(0,'defaultLineLineWidth',1.5);   
set(0,'defaultLineMarkerSize',9);
set(gca, 'FontSize', 10, 'LineWidth', 1);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

plot(t, model_true(1,:), 'b')
plot(t, data(1,:), 'ro')
plot(t, model_out(1,:), 'r')
plot(t, mean_x_predict, 'k')
plot(t, lower_predict_x, '--k')
plot(t, upper_predict_x, '--k')

legend('True Model', 'Data', 'Best Fit Model', 'Mean', '95% Prediction Interval')

xlabel('Years'), ylabel('Susceptible brain macrophages / gram')

hold off

savefig("post_predictive_solution_SIV_Susceptible.fig");
close all

hold on
set(0,'defaultLineLineWidth',1.5);   
set(0,'defaultLineMarkerSize',9);
set(gca, 'FontSize', 10, 'LineWidth', 1);
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

plot(t, model_true(2,:), 'b')
plot(t, data(2,:), 'ro')
plot(t, model_out(2,:), 'r')
plot(t, mean_y_predict, 'k')
plot(t, lower_predict_y, '--k')
plot(t, upper_predict_y, '--k')

legend('True Model', 'Data', 'Best Fit Model', 'Mean', '95% Prediction Interval')

xlabel('Years'), ylabel('SIV-infected brain macrophages / gram')

hold off

savefig("post_predictive_solution_SIV_infected.fig");
close a