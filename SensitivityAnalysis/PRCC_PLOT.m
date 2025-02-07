%% Plot the residual of the partial regression of X (input: LHS matrix) and Y (output)
%% at column s (one time point saved). RCC Coefficients are calculated on these
%% var: labels of the parameters varied in the X (as legend)
%% The Title of the plot is the Pearson correlation coefficient of the
%% transformed data, that is  the PRCC calculated on the original data.
%% The p-value is also showed in the title
%% By Simeone Marino, June 5 2007 %%
%% Modified by Marissa Renardy and Paul Wolberg, June 8, 2022.
%
% Inputs:
%    X: The LHS Matrix, N x k, where N is the number of runs and k is the
%       number of varied parameters.
%
%   Y: The model outputs, T x N, where T is the number of time points and
%      N is the number of runs.
%
%  s: A single time point ordinal. If T is the number of time points then s is
%     a single value in the range [1, T], i.e. 1 <= s <= T. For example if T is
%     10 then s is in the range [1, 10], i.e. 1 <= s <= 10.
%
% PRCC_var: A cell array of string names of the k varied parameters. This is
%           from the settings file, and is in the Matlab workspace and result
%           .mat file (Model_LHS.mat) after running Model_LHS.m.
%
% y_var: A cell array of string names of the model outputs. This is from the
%        settings file,  and is in the Matlab workspace and result
%        .mat file (Model_LHS.mat) after running Model_LHS.m, with name
%        y_var_label.
%
% For example:
% PRCC_PLOT(LHSmatrix, V_lhs, 1, PRCC_var, y_var_label)
% PRCC_PLOT(LHSmatrix, V_lhs, 2, PRCC_var, y_var_label)

function PRCC_PLOT(X, Y, s, PRCC_var, y_var)

Y=Y(s,:);
[a k]=size(X); % Define the size of LHS matrix
Xranked=rankingN(X);
Yranked=ranking1(Y);
for i=1:k  % Loop for the whole submatrices, Zi
    c1=['LHStemp=Xranked;LHStemp(:,',num2str(i),')=[];Z',num2str(i),'=[ones(a,1) LHStemp];LHStemp=[];'];
    eval(c1);
end
for i=1:k
    c2=['[b',num2str(i),',bint',num2str(i),',r',num2str(i),']= regress(Yranked,Z',num2str(i),');'];
    c3=['[b',num2str(i),',bint',num2str(i),',rx',num2str(i),']= regress(Xranked(:,',num2str(i),'),Z',num2str(i),');'];
    eval(c2);
    eval(c3);
end
for i=1:k
    c4=['r',num2str(i)];
    c5=['rx',num2str(i)];
    [r p]=corr(eval(c4),eval(c5));
    a=['[PRCC , p-value] = ' '[' num2str(r) ' , '  num2str(p) '].'];% ' Time point=' num2str(s-1)];
    figure,plot((eval(c4)),(eval(c5)),'.'),title(a),...
            legend(PRCC_var{i}),xlabel(PRCC_var{i}),ylabel(y_var);%eval(c6);

end
