
load('ParamsFINAL.mat')


fib = [];
GA = [];
GP = []; 
RAT = []; 

for i = 1 
    clfall
X = AllParams2(1,:);     


kfi = X(1)
kfg = X(2)
kpi = X(3)
kpg = X(4)
Sagg = X(5)
Salpha = X(6)
Sbeta = X(6)
Sgamma = X(7)
Rbiv = X(8)
Sbind = X(9)

% USE THIS FOR BASE PARAMETERS
% kfi = 10^(-21)*(6.022*10^23*10^(-6));
% kfg = 2*10^(-17)*(6.022*10^23*10^(-6));
% kpi = 4*10^(-18)*(6.022*10^23*10^(-6));
% kpg =  10^(-16)*(6.022*10^23*10^(-6));
% TESTSCALE = 0;%X(5)
% alpha = 1%X(6)
% beta = 1
% gamma = 1
% scl = 0
% kpscalegpi = 1



 


output = fminconplotscriptsteric(kfi,kfg,kpi,kpg,Salpha,Sbeta,Sgamma,Sagg,Rbiv,Sbind)



end
% 

