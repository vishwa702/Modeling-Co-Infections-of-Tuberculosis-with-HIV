%years
t = 1:1:22; time = 22;


%Load the given data
load ('Matrices/AusData.mat')
%data = table2array(data);
fprintf("done")
%% Parameters to Fit

% Parameters In Model:
%   beta_n, epsilon_n, gamma_n,
%   beta_p, epsilon_p, gamma_p,
% %   delta;
% %   p;
%10^-9 10^-5 og beta_n values
widths = [log10(1e-9) log10(1e-8);log10(1e-9) log10(1e-1); log10(1e-4) log10(1e-1); log10(1e-8) log10(1e-6);log10(1e-11) log10(1e-9); log10(1e-5) log10(1); log10(1e-7) log10(1e-5); log10(1e-2) log10(1)];

% % % % % % % % % % % % % Trial with some fixed values
% beta_n = 1.2e-7;
% epsilon_n = 1e-10;
% gamma_n = 1e-1;
% lambda_n = 1e-2;
% 
% beta_p = 1.8e-7;
% epsilon_p = 1.5e-10;
% gamma_p = 0.8e-1;
% lambda_p = 0.8e-2;
% 
% delta = 1.5e-9;

% widths = [1e-6 2e-5; 1e-10 2e-10; 1e-1 2e1;    1e-7 2e-7; 1e-10 1e-10; 1e-1 2e-1; 1e-2 1e-1; 1e-10 1e-7;];

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

%'par_names' is a string array holding the names of the parameters
par_names = ["beta_n", "epsilon_n", "gamma_n", "beta_p", "epsilon_p", "gamma_p", "delta" ,"p"];
fprintf("\nParameters to Fit done")
%% Parameter Fit Setup
% clc;

%Given the vector of parameters theta = [x0, y0, log10_beta, p], logl 
%returns the log likelihood distribution
logl = @(theta) (loglike_SIR_model(theta, data));

%'logp' is the natural logorithm of the prior distribution for the 
%parameters. The prior distribution is a product of uniform
%distributions, the resulting joint distribution does integrate to one.
logp = @(theta) (log(prior1(theta(1))) + log(prior2(theta(2))) + log(prior3(theta(3))) + log(prior4(theta(4))) + log(prior5(theta(5))) + log(prior6(theta(6))) + log(prior7(theta(7))));

theta0 = zeros(7, 15); % Good Rule of Thumb: chains = 2 times parameters that you will fit
% 
%Generate the vector of initial conditions

k = 0;
j = 1;

while(j < (16+1))

while(k == 0 || isnan(k) || isinf(k))
    %widths = [1e-9 1e-5;1e-9 1e-1; 1e-7 1e1; 1e-8 1e-6;1e-11 1e-9; 1e-2 1; 1e-10 1e-8; 1e-2 1];
    % initial(4) <= 0.001
    initial = [widths(1,1)+(widths(1,2)-widths(1,1))*rand widths(2,1)+(widths(2,2)-widths(2,1))*rand widths(3,1)+(widths(3,2)-widths(3,1))*rand widths(4,1)+(widths(4,2)-widths(4,1))*rand widths(5,1)+(widths(5,2)-widths(5,1))*rand widths(6,1)+(widths(6,2)-widths(6,1))*rand widths(7,1)+(widths(7,2)-widths(7,1))*rand widths(8,1)+(widths(8,2)-widths(8,1))*rand];       
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
[models, logP] = gwmcmc(theta0, logfuns, 500000);

%'models' holds the proposed guesses for the parameter values
save parameter_estimation_models2.mat models

%'logP' holds the corresponding value of the log prior and the value of the
%log likelihood for each guess of the parameter values
save parameter_estimation_logP2.mat logP
fprintf("\nAffine Fitting Algorithm done")
%% Assess results

%Let 'l' hold the values of the log prior and log likelihood
l = logP;

%Set each of the eight chains to have a burn-in of 2500 samples
l(:,:,1:2500)=[];

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
m(:,:,1:2500)=[];

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
% k8 = m(8,:);


%parameter values at the maximum log unnormalized posterior value
[M, I] = max(l_un_post_samples);
fprintf("\nAssess results done")
%% Parameter Posterior Plots

%new delta

close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 7)
savefig("post_surf_delta.fig");

[l_HPD, u_HPD] = HPD(k7, 0.05);
k7_atmax=theta_samples(7,I);
k7_median=median(k7);
% delta atmax is outside > u_HPD. Using median log value = -9.0649 i.e. 8.6119e-10 


%Beta_n

close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 1)
savefig("post_surf_beta_n.fig");
    %widths = [1e-9 1e-5;1e-9 1e-1; 1e-7 1e1; 1e-8 1e-6;1e-11 1e-9; 1e-2 1; 1e-10 1e-8; 1e-2 1];

[l_HPD, u_HPD] = HPD(k1, 0.05);
k1_atmax=theta_samples(1,I);
k1_median=median(k1);

% 95% CI: (2.2917e+06, 2.3019e+06)
% Parameter at Maximum log unnormalized posterior value: 2.2967e+06
% Parameter at median: 2.2968e+06

% Epsilon_n
close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 2)
savefig("post_surf_epsilon_n.fig");

[l_HPD, u_HPD] = HPD(k2, 0.05);
k2_atmax=theta_samples(2,I);
k2_median=median(k2);

% 95% CI: (4.6069,19.1019)
% Parameter at Maximum log unnormalized posterior value: 8.8646
% Parameter at median: 10.7186

%gamma_n
close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 3)
savefig("post_surf_gamma_n.fig");

[l_HPD, u_HPD] = HPD(k3, 0.05);
k3_atmax=theta_samples(3,I);
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

[l_HPD, u_HPD] = HPD(k4, 0.05);
k4_atmax=theta_samples(4,I);
k4_median=median(k4);



% epsilon_p
close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 5)
savefig("post_surf_epsilon_p.fig");

[l_HPD, u_HPD] = HPD(k5, 0.05);
k5_atmax=theta_samples(5,I);
k5_median=median(k5);



% gamma_p
close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 6)
savefig("post_surf_gamma_p.fig");

[l_HPD, u_HPD] = HPD(k6, 0.05);
k6_atmax=theta_samples(6,I);
k6_median=median(k6);


% delta
close all;
post_silhouette_plot(theta_samples,l_un_post_samples, par_names, 7)
savefig("post_surf_delta.fig");

[l_HPD, u_HPD] = HPD(k7, 0.05);
k7_atmax=theta_samples(7,I);
k7_median=median(k7);
% delta atmax is outside > u_HPD. Using median log value = -9.0649 i.e. 8.6119e-10 

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
% beta_n = 2.5e-9;
% epsilon_n = 8.5e-8;
% gamma_n = 0.25e-1;
% %gamma_n = 0.0396;
% % % lambda_n = 1e-2;
% % 
%  beta_p = 4e-5;
%  epsilon_p = 7e-1;
%  gamma_p = 6e-6;
% % lambda_p = NA;
% delta = 6e-6;
par = [beta_n, epsilon_n gamma_n beta_p epsilon_p gamma_p delta];
inital_values = [19005317 1090 1037 10318 110 91];


model_out = mySIR(par, inital_values, time);

% subplot(2,1,1)
close all;
subplot(1,3,1)
hold on
plot(t, data(1,:), 'ko')
plot(t, model_out(1,:))
hold off
subplot(1,3,2)
hold on
plot(t, data(2,:), 'ko')
plot(t, model_out(2,:))
hold off
subplot(1,3,3)
hold on
plot(t, data(3,:), 'ko')
plot(t, model_out(3,:))
hold off

%widths = [log10(1e-9) log10(1e-8);log10(1e-9) log10(1e-1); log10(1e-4) log10(1e-1); log10(1e-8) log10(1e-6);log10(1e-11) log10(1e-9); log10(1e-5) log10(1); log10(1e-10) log10(1e-8); log10(1e-2) log10(1)];

% Parameters In Model:
%   beta_n, epsilon_n, gamma_n,
%   beta_p, epsilon_p, gamma_p,
% %   delta;
% %   p;

% beta_p = -6;
% gamma_p = -0.5;
% epsilon_p = -0.5;
% par = [beta_n, epsilon_n gamma_n beta_p epsilon_p gamma_p delta];
% model_out = mySIR(par, inital_values, time);

close all
subplot(1,3,1)
hold on
plot(t, data(4,:), 'ko')
plot(t, model_out(4,:))
hold off
subplot(1,3,2)
hold on
plot(t, data(5,:), 'ko')
plot(t, model_out(5,:))
hold off
subplot(1,3,3)
hold on
plot(t, data(6,:), 'ko')
plot(t, model_out(6,:))
hold off


% set(0,'defaultLineLineWidth',1.5);   
% set(0,'defaultLineMarkerSize',9);
% set(gca, 'FontSize', 10, 'LineWidth', 1);
% set(0, 'DefaultAxesFontName', 'Arial');
% set(0, 'DefaultTextFontName', 'Arial');
% 
% plot(t, data(3,:), 'ko')
% plot(t, model_out(3,:))

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

x_predict = zeros(length(data),length(l_un_post_samples));
y_predict = zeros(length(data),length(l_un_post_samples));

discrepancy = zeros(1,length(l_un_post_samples));

discrepancy_pred = zeros(1,length(l_un_post_samples));

ind_D_pred_exceeds = zeros(1,length(l_un_post_samples));

h = waitbar(0,'Initialize...');
for i = 1:length(l_un_post_samples)
    
    theta = theta_samples(:,i);
    
    out = SIV_brain_model(theta(1), theta(2), 10^theta(3));
    
    discrepancy(1,i) = sum(((data(1,:) - out(1,:)).^2)./(out(1,:)/(theta(4)))) + sum(((data(2,:) - out(2,:)).^2)./(out(2,:)/(theta(4))));
    
    x_predict(:,i) = multnbinrnd(out(1,:)', theta(4), length(data(1,:)));
    y_predict(:,i) = multnbinrnd(out(2,:)', theta(4), length(data(2,:)));
    
    discrepancy_pred(1,i) = sum((((x_predict(:,i))' - out(1,:)).^2)./(out(1,:)/(theta(4)))) + sum((((y_predict(:,i))' - out(2,:)).^2)./(out(2,:)/(theta(4))));
    
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