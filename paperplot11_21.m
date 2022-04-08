function paperplot11_21
clear all 
% close all
% clfall

% 
% load('BestSol15minNEW','X','R')
load('ParamsFINAL')
X = AllParams2(1,:)
% X15=X
kfi = X(1)
kfg = X(2)
kpi = X(3)
kpg = X(4)
TESTSCALE = X(5)
alpha = X(6)
beta = X(6)
gamma = X(7)
scl = X(8)
kpscalegpi = X(9)


% load('BestSolOLD','X','R')
% kfi = X(1)
% kfg = X(2)
% kpi = X(3)
% kpg = X(4)
% TESTSCALE = X(5)
% alpha = 1
% beta = 1
% gamma = X(6)
% scl = X(7)
% kpscalegpi = X(8)


FS=32;


fignum=3;

%final time in seconds
tend=120*60;

Thrombin=[.1,.01,.001,.0001];%, .00005];
% Thrombin=[.1,.05,.01,.001];
Fibrinogen=[0.1 1 3 8 15]; 
% fracgp=[0,.3,1,2];
fracgp=[2,1,0.3,0];

FibType= 4; % controls the type of fibrin 0 is AA, 2 is WT, 3 is gammaAgammaP, 4 is gammap/gammap

kacat=84*1.5;
kbcat=49*1;
kma=7.2*1;
kmb=7.5/1; 


% kacat=84*2;
ka=10; 

% params = [1e-2 1e-1 1 1e1 1e2];
% params = [.2 .5 1 2 5];
% ThromScale = 1;
ThromScale = 1;

%Throbin cleavage
% ThromScale = 1;

%Which fibrinopeptides are cleaved fpa for just fpa cleavage, both for
%both, and none, for original fibrin activation (non-enzymatic)
fpABcleavage = 'both'; 


%relative fpb cleavage by E domain and bivalently bound thrombin alpha,





abg=1 ; % Which case to choose abg  is the index for alpha....
% 1 is default, 2 is all B, 3 is all E, 4 is B enhanced, 5 is B slowed 

KPGSCALES = [1];% 0 1 0]; %1 for all olgimerization steps, 0 for OG
KPISCALES = [1];% 0 0 1];

% both 0 is the orginal oligimerization, both 1 is the new oligomerization, can be toggled individually.  
kind=1; 

%Do the KPI and KPG terms scale with size of oligomers 'yes'. 'no' are
%options
RatesScale = 'yes';


figure(1)
[XX,b] =sort(Rs1);
Xs1= Xs1(b,:);
boxplot(Xs1(1:100,:))

set(gca,'XTickLabel',{'Kfi','Kfg','Kpi','Kpg','S_agg','S_alpha / S_beta','S_gamma','R_biv','S_bind'})
% set(gca,'XTickLabel',{'Kfi','Kfg','Kpi','K_{pg}','S_{agg}','S_{\alpha} / S_{\beta}','S_{\gamma}','R_{biv}','S_{bind}'},'Interpreter','latex')

set(gca, 'yscale','log')
set(gca, 'linew',.5)
set(gca,'FontSize',36)
% , 'LineWidth',3

a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
% idx=strcmpi(t,'box');  % Find Box objects
% boxes=a(idx);          % Get the children you need
boxes = a; 
set(boxes,'linewidth',3);

% xlabel('Parameter')
ylabel('Parameter Value (log-scale)')
hold on 


% return
for i =[1 2 4]% 1:4
    FibType = i; 
%      
% if i == 1 
%     kacat=84*1.1;
% kma=7.2/1;
% else if i == 2
%     kacat=84*1.2;
% kma=7.2/1;
%     else
%      kacat=84 * 1.3;
% kma=7.2/1;
%     end
% end
% 
if i == 1 
   clr = [0 0.4470 0.7410];
else if i == 2
    clr = [0.8500 0.3250 0.0980];
    else
    clr = [0.9290 0.6940 0.1250];
    end
end


    
    
    if i  == 4 
        alpha =1;
        beta = 1;
    end
    

sol = GeneralWeiselDriverTest(tend, 0.0001, 0.29412, fracgp(FibType),...
        kacat, kbcat, kma, kmb, kpi, kpg, kfi, kfg, ka, ThromScale,...
        fpABcleavage, alpha, beta, gamma, KPGSCALES(kind), KPISCALES(kind), RatesScale,TESTSCALE,scl,kpscalegpi); 
    
    time = sol.time/60;

  
figure(5)
subplot(2,2,1)
hold on
plot(time,  (sol.fab(1)- sol.fab)/sol.fab(1),'Color',clr,'LineStyle','--')
title('fpa cleavage %')
xlabel('time(min)')
ylabel('% Fpa cleaved')
% legend('2 (\gamma'' BS''s)/(E domain BS''s)','0.3 (\gamma'' BS''s)/(E domain BS''s)','0 (\gamma'' BS''s)/(E domain BS''s)')
xlim([0 time(end)])
ylim([0 1])
set(gca,'FontSize',36)
 
% subplot(1,2,2)
hold on
plot(time, 1-(sol.fab(1)- sol.f - sol.cf - sol.cfn - sol.cfr)/sol.fab(1),'Color',clr)
title('fpa/fpb cleavage %')
xlabel('time(min)')
ylabel('% FpA/B cleavage')
% legend('FPA - 2 (\gamma'' BS''s)/(E domain BS''s)','FPB - 2 (\gamma'' BS''s)/(E domain BS''s)'...
%     ,'FPA  - 1 (\gamma'' BS''s)/(E domain BS''s)','FPB - 1 ( \gamma'' BS''s)/(E domain BS''s)'...
%     ,'FPA - 0 (\gamma'' BS''s)/(E domain BS''s)','FPB - 0 (\gamma'' BS''s)/(E domain BS''s)')
legend('FpA - 2:1 \gamma'' per E-domain','FPB - 2:1 \gamma'' per E-domains'...
    ,'FpA  - 1:1 \gamma'' per E-domain','FPB - 1:1 \gamma'' per E-domain'...
    ,'FpA - 0:1 \gamma'' per E-domain','FPB - 0:1 \gamma'' per E-domain')
% legend('FPA - kacat, km','FPB - kacat, km'...
%     ,'FPA - kacat/2, km','FPB - kacat/2, km'...
%     ,'FPA  - kacat ,km x 2','FPB - kacat,km x 2'...
%     )

xlim([0 time(end)])
ylim([0 1])
set(gca,'FontSize',36)


subplot(2,2,2)
hold on
plot(time, sol.T./sol.T(1))
title('Fraction Free Thrombin')
xlabel('time(min)')
ylabel('Fraction Free Thrombin')
legend('2:1  \gamma'' per E-domain'...
    ,'1:1 \gamma'' per E-domain'...
    ,' 0:1 \gamma'' per E-domain')
% legend('FPA - kacat, km'...
%      ,'FPA - kacat/2, km'...
%     ,'FPA  - kacat,km x 2'...
%    )
xlim([0 time(end)])
ylim([0 1])
set(gca,'FontSize',36)

subplot(2,2,3)
hold on
plot(time, (sol.Efab + sol.Efb + sol.BE1 + sol.BG1)./sol.T(1))
title('Low affintity associated thrombin')
xlabel('time(min)')
ylabel('Fraction LA-bound Thrombin')
legend('2:1  \gamma'' per E-domain'...
    ,'1:1 \gamma'' per E-domain'...
    ,' 0:1 \gamma'' per E-domain')
% legend('FPA - kacat, km'...
%      ,'FPA - kacat/2, km'...
%     ,'FPA  - kacat,km x 2'...
%    )
xlim([0 time(end)])

ylim([0 1])
set(gca,'FontSize',36)

subplot(2,2,4)
hold on
plot(time, (sol.BG + sol.BE + sol.B)./sol.T(1))
title('High affinity associated thrombin')
xlabel('time(min)')
ylabel('Fraction HA-bound Thrombin')
legend('2:1  \gamma'' per E-domain'...
    ,'1:1 \gamma'' per E-domain'...
    ,' 0:1 \gamma'' per E-domain')
% legend('FPA - kacat, km'...
%      ,'FPA - kacat/2, km'...
%     ,'FPA  - kacat,km x 2'...
%    )
xlim([0 time(end)])
ylim([0 1])
set(gca,'FontSize',36)

% 
% figure(6)
% hold on
% plot(time, (sol.fab/sol.fab(1)))
% title('Total Fibrinogen')
% xlabel('time(min)')
% ylabel('\% Fpb cleaved')
% legend('2 (\gamma'' BS''s)/(E domain BS''s)'...
%     ,'1 (\gamma'' BS''s)/(E domain BS''s)'...
%     ,'wt (\gamma'' BS''s)/(E domain BS''s)'...
%     ,' 0 (\gamma'' BS''s)/(E domain BS''s)')
% xlim([0 time(end)])
% ylim([0 1])
% 
% 
% figure(7)
% hold on
% plot(time, (sol.fab/sol.fab(1)))
% title('Total Fibrinogen')
% xlabel('time(min)')
% ylabel('\% Fpb cleaved')
% legend('2 (\gamma'' BS''s)/(E domain BS''s)'...
%     ,'1 (\gamma'' BS''s)/(E domain BS''s)'...
%     ,'0.3 (\gamma'' BS''s)/(E domain BS''s)'...
%     ,' 0 (\gamma'' BS''s)/(E domain BS''s)')
% xlim([0 time(end)])
% ylim([0 1])





  
'fibrinogen'
sol.fab(end)/sol.fab(1) + sol.Efab(end)/sol.fab(1)
'Thrombin'
sol.T(end)/sol.T(1)
% pause
end

end

% kacat=pe(1);
% kap=pe(2);
% kam=pe(3);
% kbcat=pe(4);
% kbp=pe(5);
% kbm=pe(6);
    