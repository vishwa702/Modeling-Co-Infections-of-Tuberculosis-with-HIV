% beta=par(1);
% gamma=par(2);
% S0=par(3);
% I0=par(4);
% R0=par(5);

%beta range: 4e-8 4e-7
%gamma range: 0.01 0.1
%I0 range: 10 100

PRCC_var = ["beta","gamma","I0"];

[sens_param,~]=lhsdesign_modified(10000,[4e-8 0.01 10],[4e-7 0.1 100]);

bcol=sens_param(:,1);
gcol=sens_param(:,2);
Icol=sens_param(:,3);

LHSmatrix=[bcol,gcol,Icol];

final_size=zeros(1,10000);
time=150;

for i=1:10000

    par=[bcol(i) gcol(i) 4000000-Icol(i) Icol(i) 0];
    final_size(1,i)=SIR_se(par,time);
end

[prcc,sign,sign_label]=PRCC(LHSmatrix,final_size,1,PRCC_var,0.05);
X=categorical(PRCC_var);
X=reordercats(X,PRCC_var);
bar(X,prcc)
ylabel('PRCC final size')
ylim([-1 1])