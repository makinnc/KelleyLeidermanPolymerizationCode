function FpbTestDriverAriens(kfi,kfg,kpi,kpg,kacat,kbcat,kma,kmb,alpha,beta,gamma,ThromScale,TESTSCALE,scl,kpscalegpi)
% clear all; 
% clfall; 
clf(figure(10))
% close all;


runname={'_or_fpball',...
        '_or_fpbEonly',...
        '_or_fpbBonly',...
        '_or_fpbBenhanced',...
        '_or_fpbBslowed'};

    
   FS=32;
fignum=3;

%final time in seconds
tend=60*60*24;%1/3*60*60;


Thrombin=U2nmole([0.005 0.01 0.05 0.1 1 2 10])*10^(-3);
Fibrinogen=fibrin_mgml2muM(1);
% Thrombin=[.1,.01,.001,.0001];%, .00005];
% Thrombin=[.1,.05,.01,.001];
% Fibrinogen=[0.1 1 3 8 15]; 
fracgp=[0,.3,1,2];

FibType= 1; % controls the type of fibrin 0 is AA, 2 is WT, 3 is gammaAgammaP, 4 is gammap/gammap




% kacat=84*2;
% kbcat=49*1;
% kma=7.2/2;
% kmb=7.5/1; 
% 
% kpi=100; %100;
% kpg=150; %150; 
% kfi= 5*10^(-21)*(6.022*10^23*10^(-6)) ; 
% kfg=5*12;
ka=10; 

% params = [1e-2 1e-1 1 1e1 1e2];
% params = [.2 .5 1 2 5];

% kpi=params*kpi;
% kpg=params*kpg; 
% kfi= params*kfi; 
% kfg=params*kfg;
% ka=params*ka; 

%Throbin cleavage
% ThromScale = 1;

%Which fibrinopeptides are cleaved fpa for just fpa cleavage, both for
%both, and none, for original fibrin activation (non-enzymatic)
fpABcleavage = 'both'; 


%relative fpb cleavage by E domain and bivalently bound thrombin alpha,
%beta are for E domain, gamma is for bivalently bound.
% alpha=[1 0 1 .1 1]; 
% beta=[1 0 1 .1 1]; 
% gamma=[1 1 0 100 0.1];%linspace(1,50, 5); 

abg=1 ; % Which case to choose abg  is the index for alpha....
% 1 is default, 2 is all B, 3 is all E, 4 is B enhanced, 5 is B slowed 

KPGSCALES = [1 0 1 0]; %1 for all olgimerization steps, 0 for OG
KPISCALES = [1 0 0 1];

% both 0 is the orginal oligimerization, both 1 is the new oligomerization, can be toggled individually.  
kind=1; 

%Do the KPI and KPG terms scale with size of oligomers 'yes'. 'no' are
%options
RatesScale = 'yes';

% pf=[0, 30,60,90,120,150,180,210,240];
% diameter=[54, 75.5, 92, 105.5, 118, 128.7,139.1,148.7];
% numb_innac=[0, 4,17,35,54,75,96,120,141];
% frac_innac=numb_innac./pf;
% frac_innac(1)=0;
% 
% FS=32;
% fignum=1;
% 
% kacat=84;%84*4;
% kbcat=49;%49*1;
% kma=7.2;%/4;
% kmb=7.5;%/4; 
% % scalekpi=100;
% 
% 
% kpi=100;
% kpg=150; 
% kfg=12;
% 
% 
% %final time in seconds
% tend=12*60*60; %12 hrs
% 
% 
% 
% Thrombin=U2nmole([0.005 0.01 0.05 0.1 1 2 10])*10^(-3);
% Fibrinogen=fibrin_mgml2muM(1); 
% fracgp=[0,.3,1,2];
% 
% 
% alpha=[1 0 1 1 1]; 
% beta=[1 0 1 1 1]; 
% gamma=[1 1 0 10 0.1];
% KPGSCALES = [1 0 1 0]; %1 for all olgimerization steps, 0 for OG
% KPISCALES = [1 0 0 1];
% 
% % 
% % 
% % for i=1:length(Thrombin)
% %      for j=1:length(fracgp)
% 
% names={'100 nM thrombin','10 nM thrombin','1 nM thormbin','0.1 nM thormbin'};


% ProteinDensity= ones(2,length(Thrombin)); 
% ProtNumber= ones(2,length(Thrombin));
% ProtDistance= ones(2,length(Thrombin));

for k=1;%:5
%     clfall
Thrombin=U2nmole([0.005 0.01 0.05 0.1 1 2 10])*10^(-3);
count=0; 
for j=[1 3] %:length(fracgp)
    count
    count=count+1; 
    for i=1:length(Thrombin)

    
% Init_conds= {'T_IC',Thrombin(i);
%              'fab_IC',Fibrinogen};
% [nICs,dummy]= size(Init_conds);
%          
%          
% rates={'fracgp',fracgp(j);
%         'kbcat',kbcat;
%     'kacat',kacat;
%     'kmb',kmb;
%     'kma',kma;
%     'kpi', kpi ;%scalekpi*4*10^(-18)*(6.022*10^23*10^(-6));
%     'kfg',kfg;
%     'kpg', kpg;
%     'alpha',alpha(k);
%     'beta',beta(k);
%     'gamma',gamma(k);
%     'kpiscale', KPISCALES(1)
%     'kpgscale', KPGSCALES(1)};
% [nrates,dummy]= size(rates);
% sol=PolymModel(tend,nICs,nrates,Init_conds,rates);
            Ti =i;
%             Fi = 1;
%             kpii= 3; 
%             kpgi=3;
%             kfii=3; 
%             kfgi=3;
%             kai=3;
            FibType=j; 

sol = GeneralWeiselDriverTest(tend, Thrombin(Ti), Fibrinogen, fracgp(FibType),...
        kacat, kbcat, kma, kmb, kpi, kpg, kfi, kfg, ka, ThromScale,...
        fpABcleavage, alpha, beta, gamma, KPGSCALES(kind), KPISCALES(kind), RatesScale,TESTSCALE,scl,kpscalegpi); 

% 
% ind=(j-1)*4;
% FpBTestPlot(sol,FS,)
% fignum=[1,2,3,4,5]; 
% fignum=[1+ind,2+ind,3+ind,4+ind]; 
% FpBTestPlot(sol,FS,fignum, names)
% FibrinSpeciesPlot(sol,FS,fignum,names)


R = [102 102 100 100 90 85 85; 90 90 95 94 90 91 90]; 


V = FiberVolume2(sol, R(count, i)) ; %nm^3
Vp =  (pi * (R(count, i)).^2)./sol.m(end) ; 


NumFib = sol.m(end);%*sol.l(end); 

% ProteinDensity(count, i) = NumFib*(134/6.022 * 10^(-20)) / (V *(10^(-21))); 
ProteinDensity(count, i) = NumFib*(340/6.022 * 10^(-20)) / (V *(10^(-21))); 

ProtNumber(count, i)= sol.m(end); %#protofibril
ProtDistance(count, i)= sqrt(Vp);%2*sqrt(Vp/pi); %nm


% sol.E + sol.BE
% pause
% sol.m(end)
% abs(sol.fab(1) - sol.cfr(end))
    end
%     resetcolorall
end


% Thrombin = categorical({'0.005 U/ml' '0.01 U/ml' '0.05 U/ml' '0.1 U/ml' '1 U/ml' '2 U/ml' '10 U/ml'}); 
Thrombin = categorical({'0.005' '0.01' '0.05' '0.1' '1' '2' '10'}); 
figure(10); %fignum(end)+1) 

subplot(2,2,1)
hold on
% bar(Thrombin, R')
bar(Thrombin, R')
ylabel('Fiber Radius (nm)')
xlabel('Thrombin (U/ml)' )
legend('\gamma_A/\gamma_A','\gamma_A/\gamma''')
set(gca,'FontSize',28)
title('Fiber Radius')


subplot(2,2,2)
hold on
% bar(Thrombin, ProtNumber' )
bar(Thrombin,ProtNumber' )
ylabel('Protofibril Number')
xlabel('Thrombin (U/ml)')
legend('\gamma_A/\gamma_A','\gamma_A/\gamma''')
set(gca,'FontSize',28)
title('Protofibril Number')


subplot(2,2,3)
hold on
bar(Thrombin, ProteinDensity')
ylabel('Protein Density (g/cm^3)')
xlabel('Thrombin (U/ml)')
legend('\gamma_A/\gamma_A','\gamma_A/\gamma''')
set(gca,'FontSize',28)
title('Protein Density')


subplot(2,2,4)
hold on
bar(Thrombin, ProtDistance')
ylabel('Protofibril Distance (nm)')
xlabel('Thrombin (U/ml)')
legend('\gamma_A/\gamma_A','\gamma_A/\gamma''')
set(gca,'FontSize',28)
title('Protofibril Distance')

% figure(1)
% set(gcf,'PaperPosition',[0 0 30 30])
% print(['AriensBarPlot', runname{k}],'-dpng','-r0')

end
end

function sol=PolymModel(tend,nICs,nrates,Init_conds,rates)


for i=1:nICs
    
    names_IC{i}=Init_conds{i,1};
     params_IC.(names_IC{i})=Init_conds{i,2};
end

for i=1:nrates
    
   names_rates{i}=rates{i,1};
     params_rates.(names_rates{i})=rates{i,2};
end


[sol]= weiselmodelfcnnew(tend,names_IC,params_IC,names_rates,params_rates);

end