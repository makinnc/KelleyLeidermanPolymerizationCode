clear all
% close all
 clfall; 

FS=32;
fignum=3;

%final time in seconds
tend=1/3*60*60;


load('ParamsFINAL.mat')
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


Thrombin=0.001;
% Thrombin=[.1,.05,.01,.001];
Fibrinogen=1  ;
fracgp = 0.3
% fracgp=[0,.3,1,2];
% 
% FibType= 1; % controls the type of fibrin 0 is AA, 2 is WT, 3 is gammaAgammaP, 4 is gammap/gammap




kacat=84*4;
kbcat=49*1;
kma=7.2/4;
kmb=7.5/1; 

kpi=100;
kpg=150; 
kfi= 5*10^(-21)*(6.022*10^23*10^(-6)) ; 
kfg=12*5;
ka=0.1; 

% params = [1e-2 1e-1 1 1e1 1e2];
params = [.2 .5 1 2 5];

kpi=params*kpi;
kpg=params*kpg; 
kfi= params*kfi; 
kfg=params*kfg;
ka=params*ka; 

%Throbin cleavage
ThromScale = 1;

%Which fibrinopeptides are cleaved fpa for just fpa cleavage, both for
%both, and none, for original fibrin activation (non-enzymatic)
fpABcleavage = 'both'; 


%relative fpb cleavage by E domain and bivalently bound thrombin alpha,
%beta are for E domain, gamma is for bivalently bound.
alpha=[1]; 
beta=[1]; 
gamma=[1];%linspace(1,50, 5); 

abg=1 ; % Which case to choose abg  is the index for alpha....
% 1 is default, 2 is all B, 3 is all E, 4 is B enhanced, 5 is B slowed 

% KPGSCALES = [1 0 1 0]; %1 for all olgimerization steps, 0 for OG
% KPISCALES = [1 0 0 1];
% 
% % both 0 is the orginal oligimerization, both 1 is the new oligomerization, can be toggled individually.  
% kind=1; 

%Do the KPI and KPG terms scale with size of oligomers 'yes'. 'no' are
%options
RatesScale = 'yes';




% fracgp

% GeneralWeiselDriverTest(tend, Thrombin, Fibrinogen, fracgp, kacat, kbcat, kma, kmb, kpi, kpg, kfi, kfg, ka, ThromScale, fpABcleavage, alpha, beta, gamma, KPGSCALES, KPISCALES, RatesScale, TESTSCALE,scl,kpscalegpi)
%select whether or not FpB and FpB are being considered
switch fpABcleavage
    case 'fpa'
% to turn off fpb cleavage set kbcat = 0; 

        kbcat = 0;
        ka=0; 
        
    case 'both'
        ka=0; 
        
    case 'none'   
%T0 turn off enzymatic cleavage of thrombin and make it kfa*fa set kbcat = kacat = 0
%set kfa too

    kbcat = 0; 
    kacat = 0; 
%     ka = .1;     
end

    
nICs=2; %number of initial cond supplied. 
% names_IC={};
Init_conds= {'T_IC',Thrombin;
             'fab_IC',Fibrinogen};


rates={'fracgp',fracgp;
    'kpi', kpi ;%scalekpi*4*10^(-18)*(6.022*10^23*10^(-6));
    'kfi', kfi ;%scalekpi*4*10^(-18)*(6.022*10^23*10^(-6));
    'kfg',kfg;
    'kpg', kpg;
    'alpha',Salpha;
    'beta',Sbeta;
    'gamma',Sgamma; 
    'ThromScale', ThromScale;
    'TESTSCALE', Sagg;
    'scl', Rbiv;
    'kpscalegpi', Sbind;};
[nrates,dummy]= size(rates);




sol=PolymModel(tend,nICs,nrates,Init_conds,rates, RatesScale);





 
 

function sol=PolymModel(tend,nICs,nrates,Init_conds,rates,scales)


for i=1:nICs
    
    names_IC{i}=Init_conds{i,1};
     params_IC.(names_IC{i})=Init_conds{i,2};
end

for i=1:nrates
    
   names_rates{i}=rates{i,1};
     params_rates.(names_rates{i})=rates{i,2};
end


[sol]= weiselmodelfcnnew(tend,names_IC,params_IC,names_rates,params_rates, scales);

end

% pf=[0, 30,60,90,120,150,180,210,240];
% diameter=[54, 75.5, 92, 105.5, 118, 128.7,139.1,148.7];
% numb_innac=[0, 4,17,35,54,75,96,120,141];
% frac_innac=numb_innac./pf;
% frac_innac(1)=0;