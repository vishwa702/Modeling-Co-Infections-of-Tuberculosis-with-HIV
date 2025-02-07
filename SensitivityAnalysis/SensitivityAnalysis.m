%years
t = 1:1:22; 
time = 22;


%Load the given data
load ('data_australia.mat')
%data = table2array(data);
data = double(data);
fprintf("Load Data Done")

%%


% par =  [beta_n epsilon_n gamma_n beta_p epsilon_p gamma_p delta S_n_0 I_n_0 R_n_0 S_p_0 I_p_0 R_p_0]
                  % beta_n                  epsilon_n             gamma_n                     beta_p                  epsilon_p              gamma_p                    delta                    p                         lambda_n               lambda_p              L_n_0              L_p_0 

%PRCC_var = ["beta_n","epsilon_n","gamma_n","beta_p","epsilon_p","gamma_p","delta","lambda_n","lambda_p","L_n_0","L_p_0"];
PRCC_var = ["beta_n","lambda_n","gamma_n","epsilon_n","beta_p","lambda_p","gamma_p","epsilon_p","delta","L_n_0","L_p_0"];

%generating the parameter values
widths = [log10(1e-10) log10(1e-6); log10(1e-4) log10(1e0); log10(1e-4) log10(1e0); log10(1e-8) log10(1e-1); log10(1e-3) log10(1e3); log10(1e-3) log10(1e3); log10(1e-10) log10(1e-6); log10(1e-3) log10(1e1) ; log10(1e-3) log10(1e1); log10(1e-8) log10(1e1); log(1e1) log(1e6); log(1e1) log(1e8)];
[sens_param,~]=lhsdesign_modified(10000,[-10 -4 -4 -8 -3 -3 -10 -3 -3 -8 1 1],[-8 -2 -2 -6 3 3 -6 1 1 1 3 3]);
%[sens_param,~]=lhsdesign_modified(10000,[-10 -4 -4 -8 -3 -3 -10 -3 -3 -8 1 1],[-6 0 0 -1 3 3 -6 1 1 1 6 8]);


beta_n_col = sens_param(:,1);
epsilon_n_col = sens_param(:,2);
gamma_n_col = sens_param(:,3);
beta_p_col = sens_param(:,4);
epsilon_p_col = sens_param(:,5);
gamma_p_col = sens_param(:,6);
delta_col = sens_param(:,7);
lambda_n_col = sens_param(:,8);
lambda_p_col = sens_param(:,9);
L_n_0_col = sens_param(:,10);
L_p_0_col = sens_param(:,11);

LHSmatrix=[beta_n_col,epsilon_n_col,gamma_n_col,beta_p_col,epsilon_p_col,gamma_p_col,delta_col,lambda_n_col,lambda_p_col,L_n_0_col,L_p_0_col];

final_size=zeros(1,10000);
time=22;

L_n_0 = 1770;
L_p_0 = 1120;



%generate output for each parameter value 
for i=1:10000
    par = [beta_n_col(i), epsilon_n_col(i), gamma_n_col(i),beta_p_col(i),epsilon_p_col(i),gamma_p_col(i),delta_col(i),0,lambda_n_col(i),lambda_p_col(i)];
    initial_values = [data(1,1) data(2,1) data(3,1) data(4,1) data(5,1) data(6,1) L_n_0_col(i) L_p_0_col(i)];
    final_size(1,i)=SLIR_sens(par, initial_values,time);
    fprintf("i = %d\n",i);
end

[prcc,sign,sign_label]=PRCC(LHSmatrix,final_size,1,PRCC_var,0.05);
%PRCC_var = ["beta_p","lambda_n","gamma_n","epsilon_p","gamma_p","beta_n","epsilon_n","delta","lambda_p","L_n_0","L_p_0"];
PRCC_var = ["beta_n","lambda_n","gamma_n","epsilon_n","beta_p","lambda_p","gamma_p","epsilon_p","delta","L_n_0","L_p_0"];

X=categorical(PRCC_var);
X=reordercats(X,PRCC_var);

%newprcc = prcc;
X = categorical(X);
newprcc(1) = prcc(1); %beta_n
newprcc(8) = prcc(2); %epsilon_n
newprcc(3) = prcc(3); %gamma_n
newprcc(2) = prcc(4); %beta_p
newprcc(4) = prcc(5); %epsilon_p
newprcc(9) = prcc(6); %gamma_p
newprcc(6) = prcc(7); %delta
newprcc(5) = prcc(8); %lambda_n
newprcc(7) = prcc(9); %lambda_p
newprcc(10) = prcc(10); %L_n_0
newprcc(11) = prcc(11); %L_p_0
close all
bar(X,newprcc);
%val = ["1" "2" "3"];
%h=bar(val,newprcc);
%h(2).FaceColor = 'red';
%h.CData(2,:) = [.5 0 .5];

ylabel('PRCC final size')
xlabel('Parameters')
title('Global Sensitivity Analysis of Logistic Growth Model')

newprcc

