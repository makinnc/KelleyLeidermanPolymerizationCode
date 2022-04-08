function [output] = fminconplotscriptsteric(kfi,kfg,kpi,kpg,alpha,beta,gamma,TESTSCALE,scl,kpscalegpi)

%Smaller param sweep
% clfall
% clear all
close all
% filename='ParamSweep2_2_1.mat'; 
testind=1;
% 
% alpha=[1 0 1 .1 1]; 
% beta=[1 0 1 .1 1]; 
% gamma=[1 1 0 1 0.1];

% alpha=[1 0 1 1 1]; 
% beta=[1 0 1 1 1]; 
% % gamma=[1 1 0 100 0.1];
% alpha=[1 0 1 1 1]; 
% beta=[1 0 1 1 1]; 
% gamma=[1 1 0 100 0.1];
% alpha = 1;
% beta = 1; 
% % filename='2ParamSweep2_ALL_1_alt_gamma.mat'
% gamma =[.01 .1 .5 1 5 10 100];
% 

kacat=84*1;
kbcat=49*1;
kma=7.2/1;
kmb=7.5/1; 

% kpi=100; %100;
% kpg=150; %150; 

% 
% kpi= 1*(4*10^(-18)*(6.022*10^23*10^(-6))); %100;
% %Protofibril growth
% kpg=1*(10^(-16)*(6.022*10^23*10^(-6)));
% 
% % kfi= 5*10^(-21)*(6.022*10^23*10^(-6)) ; 
% % kfg=5*12;
% 
% 
% kfi= 1*10^(-21)*(6.022*10^23*10^(-6)); 
% kfg= 1*12;
% 
% % ThromScale = 1;
ThromScale = 1;

% params = [.01 .1 .5 1 5 10 100];

% TESTSCALE = [0 0.25 0.5 0.75 1]; 
% TESTSCALE = TESTSCALE(testind); 

% kfi = kfi * params; 
% kfg = kfg * params;
% kpi = kpi * params; 
% kpg = kpg * params;

% 
% kfi = X(1)
% kfg = X(2)
% kpi = X(3)
% kpg = X(4)
% TESTSCALE = X(5)
% alpha = X(6)
% beta = X(6)
% gamma = X(7)
% scl = X(8)
% kpscalegpi = X(9)



% ProtNumAll = zeros(length(gamma));

tscales= {'one','two','three','four','five','six','seven'}; 

% [AA,GP] = DomData; 

% load('inds.mat')
% load('ParamSweep2_2_1.mat')

% protn=[];

% figure(1)
% plot(AA)
% hold on
% size(Cind)
% Cind = Cind([1 8 15 22 29],:);

%   h=3; 
%   g=7; 
%  ceil(Ch/7);
%   i = h;%ceil(h/7); %Bind(h,1);
%   j = Cind(h,2);
%   k = Cind(h,3);
%   l = Cind(h,4);
%   p = Cind(h,5); 
% 
%   
%    i = 2
%   j = 6;
%   k = 7;
%   l = 3;
%   p = 4; 
%   
%   g = 5; 
% 
% i = ceil(h/7); %Bind(h,1);
%   j = Cind(2,h);
%   k = Cind(1,h);
%   l = Cind(1,h);
%   p = Cind(1,h); 



output = PolymModelRunBAse(kfi,kfg,kpi,kpg,kacat,kbcat,kma,kmb,alpha,beta,gamma,ThromScale,TESTSCALE,scl,kpscalegpi);


% pause

FpbTestDriverAriens(kfi,kfg,kpi,kpg,kacat,kbcat,kma,kmb,alpha,beta,gamma,ThromScale,TESTSCALE,scl,kpscalegpi)


paperplot11_21

% P1=FpbTestDriverAriensParamSearch(kfi(k),kfg(j),kpi(p),kpg(l),kacat,kbcat,kma,kmb,alpha,beta,gamma(g),ThromScale,TESTSCALE(i));
% P2=FpbTestDriverAriensParamSearch2(kfi(k),kfg(j),kpi(p),kpg(l),kacat,kbcat,kma,kmb,alpha,beta,gamma(g),ThromScale,TESTSCALE(i));

% 'Kfi'
% kfi(k)
% 'Kfg'
% kfg(j)
% 'Kpi'
% kpi(p)
% 'Kpg'
% kpg(p)
% 'AggScale'
% TESTSCALE(i)
% 'FpbEffic'
% gamma(g)


