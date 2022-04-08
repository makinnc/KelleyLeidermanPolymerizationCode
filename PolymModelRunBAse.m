 function [output] = PolymModelRunBAse(kfi,kfg,kpi,kpg,kacat,kbcat,kma,kmb,alpha,beta,gamma,ThromScale,TESTSCALE,scl,kpscalegpi)

% clear all; 
% clfall; 

FS=32;
fignum=3;

%final time in seconds
tend=15*60;

Thrombin=[.1,.01,.001,.0001];%, .00005];
% Thrombin=[.1,.05,.01,.001];
Fibrinogen=[0.1 1 3 8 15]; 
% fracgp=[0,.3,1,2];
fracgp=[2,1,0.3,0];

FibType= 1; % controls the type of fibrin 0 is AA, 2 is WT, 3 is gammaAgammaP, 4 is gammap/gammap


output = []; 

% kacat=84*2;
% kbcat=49*1;
% kma=7.2/2;
% kmb=7.5/1; 

% kpi=100; %100;
% kpg=150; %150; 
% kfi= 5*10^(-21)*(6.022*10^23*10^(-6)) ; 
% kfg=5*12;
ka=10; 

% params = [1e-2 1e-1 1 1e1 1e2];
params = [.2 .5 1 2 5];

kpi=params*kpi;
kpg=params*kpg; 
kfi= params*kfi; 
kfg=params*kfg;
ka=params*ka; 

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



% for j=1:6
%     j
% figure(4)
% subplot(2,3,j)
% hold on

% for  i=1:5
%     i
%     
%         switch j 
%         case 1
%             Ti = 3;
%             Fi = i;
%             kpii= 3; 
%             kpgi=3;
%             kfii=3; 
%             kfgi=3;
%             kai=3; 
%         case 2
%             Ti = i;
%             Fi = 3;
%             kpii= 3; 
%             kpgi=3;
%             kfii=3; 
%             kfgi=3;
%             kai=3; 
%         case 3
%             Ti = 3;
%             Fi = 3;
%             kpii= 3; 
%             kpgi=3;
%             kfii=3; 
%             kfgi=3;
%             kai=3; 
%         case 4
%             Ti = 3;
%             Fi = 3;
%             kpii= 3; 
%             kpgi=3;
%             kfii=i; 
%             kfgi=3;
%             kai=3; 
%         case 5
%             Ti = 3;
%             Fi = 3;
%             kpii= 3; 
%             kpgi=3;
%             kfii=3; 
%             kfgi=i;
% % %             kai=3; 
% % %          case 6
% % %             Ti = 3;
% % %             Fi = 3;
% % %             kpii= 3; 
% % %             kpgi=i;
% % %             kfii=3; 
% % %             kfgi=3;
% % %             kai=3;           
% % %         end
% 
% %
% Fibrinogen = fibrin_mgml2muM(.15)
% pause
for i= 1:length(fracgp)
            Ti =3;
            Fi = 3;
            kpii= 3; 
            kpgi=3;
            kfii=3; 
            kfgi=3;
            kai=3;
            
            FibType = i; 
tend = 30*60;
sol = GeneralWeiselDriverTest(tend, U2nmole(.1)*10^(-3), fibrin_mgml2muM(.5), fracgp(FibType),...
        kacat, kbcat, kma, kmb, kpi(kpii), kpg(kpgi), kfi(kfii), kfg(kfgi), ka(kai), ThromScale,...
        fpABcleavage, alpha(abg), beta(abg), gamma(abg), KPGSCALES(kind), KPISCALES(kind), RatesScale,TESTSCALE,scl,kpscalegpi);
% sol = GeneralWeiselDriverTest(tend, Thrombin(Ti), fibrin_mgml2muM(.5), fracgp(FibType),...
%         kacat, kbcat, kma, kmb, kpi(kpii), kpg(kpgi), kfi(kfii), kfg(kfgi), ka(kai), ThromScale,...
%         fpABcleavage, alpha(abg), beta(abg), gamma(abg), KPGSCALES(kind), KPISCALES(kind), RatesScale,TESTSCALE,scl,kpscalegpi); 
    
fignum=[3,4]; 
legnd = 'BS';
fourpanelplot2alt(sol, fignum, FS,legnd) 

% figure(5)
% subplot(2,2,1)
% hold on
% plot(sol.time,  (sol.fab(1)- sol.fab)/sol.fab(1))
% title('fpa cleavage \%')
% xlabel('time(min)')
% ylabel('\% Fpa cleaved')
% switch legnd 
%      case 'BS'
%           legend('2 (\gamma'' BS''s)/(E domain BS''s)','1(\gamma'' BS''s)/(E domain BS''s)','0.3 (\gamma'' BS''s)/(E domain BS''s)','0 (\gamma'' BS''s)/(E domain BS''s)')
% 
%      case 'Thrombin'
%           legend('100 nM thrombin','10 nM thrombin','1 nM thrombin','0.1 nM thrombin')
% end
%  
% 
% subplot(2,2,2)
% hold on
% plot(sol.time, 1-(sol.fab(1)- sol.f - sol.cf - sol.cfn - sol.cfr)/sol.fab(1))
% title('fpa cleavage \%')
% xlabel('time(min)')
% ylabel('\% Fpb cleaved')
% switch legnd 
%      case 'BS'
%           legend('2 (\gamma'' BS''s)/(E domain BS''s)','1(\gamma'' BS''s)/(E domain BS''s)','0.3 (\gamma'' BS''s)/(E domain BS''s)','0 (\gamma'' BS''s)/(E domain BS''s)')
% 
%      case 'Thrombin'
%           legend('100 nM thrombin','10 nM thrombin','1 nM thrombin','0.1 nM thrombin')
%  end
% fourpanelplot2alt2(sol, 20, FS,legnd) 

% fignum=11; 
% fpareleaseplot(sol, fignum, FS,legnd) 

% figure(5)
% subplot(2,2,1)
% hold on
% subplot(2,2,2)
% hold on
% subplot(2,2,3)
% hold on
% subplot(2,2,4)
% hold on
% fignum=3; 
% fourpanelplot(sol, fignum, FS)
% figure(5)
% subplot(2,2,1)
% hold on
% plot(sol.time, sol.B)
% title('B')
% 
% subplot(2,2,2)
% hold on
% plot(sol.time, sol.BE1,sol.time, sol.BG1)
% title('BE1+BG1')
% 
% subplot(2,2,3)
% hold on
% plot(sol.time, sol.BG)
% title('BG')
% 
% subplot(2,2,4)
% hold on
% plot(sol.time, sol.BE)
% title('BE')

% pause
end
% pause


for i= 1:length(Thrombin)
            Ti =i;
            Fi = 3;
            kpii= 3; 
            kpgi=3;
            kfii=3; 
            kfgi=3;
            kai=3;
            
            FibType = 3; 
tend = 30*60;

sol = GeneralWeiselDriverTest(tend, Thrombin(Ti), Fibrinogen(Fi), fracgp(FibType),...
        kacat, kbcat, kma, kmb, kpi(kpii), kpg(kpgi), kfi(kfii), kfg(kfgi), ka(kai), ThromScale,...
        fpABcleavage, alpha(abg), beta(abg), gamma(abg), KPGSCALES(kind), KPISCALES(kind), RatesScale,TESTSCALE,scl,kpscalegpi); 
    
fignum=[3 4]; 
legnd = 'Thrombin';
fourpanelplot2alt(sol, fignum, FS,legnd) 
% fourpanelplot2alt2(sol, 21, FS,legnd) 


% figure(5)
% subplot(1,2,1)
% hold on
% plot(sol.time,  (sol.fab(1)- sol.fab)/sol.fab(1))
% title('fpa cleavage %')
% xlabel('time(min)')
% ylabel('% Fpa cleaved')
% set(gca,'FontSize',28)
% switch legnd 
%      case 'BS'
%           legend('2 (\gamma'' BS''s)/(E domain BS''s)','1(\gamma'' BS''s)/(E domain BS''s)','0.3 (\gamma'' BS''s)/(E domain BS''s)','0 (\gamma'' BS''s)/(E domain BS''s)')
% 
%      case 'Thrombin'
%           legend('100 nM thrombin','10 nM thrombin','1 nM thrombin','0.1 nM thrombin')
% end
%  
% subplot(1,2,2)
% hold on
% plot(sol.time, 1-(sol.fab(1)- sol.f - sol.cf - sol.cfn - sol.cfr)/sol.fab(1))
% title('fpb cleavage %')
% xlabel('time(min)')
% ylabel('% Fpb cleaved')
% set(gca,'FontSize',28)
% switch legnd 
%      case 'BS'
%           legend('2 (\gamma'' BS''s)/(E domain BS''s)','1(\gamma'' BS''s)/(E domain BS''s)','0.3 (\gamma'' BS''s)/(E domain BS''s)','0 (\gamma'' BS''s)/(E domain BS''s)')
% 
%      case 'Thrombin'
%           legend('100 nM thrombin','10 nM thrombin','1 nM thrombin','0.1 nM thrombin')
%  end
% fignum=12; 
% fpareleaseplot(sol, fignum, FS,legnd) 

% figure(6)
% subplot(2,2,1)
% hold on
% plot(sol.time, sol.B)
% title('B')

% subplot(2,2,2)
% hold on
% plot(sol.time, sol.BE1,sol.time, sol.BG1)
% title('BE1+BG1')
% 
% subplot(2,2,3)
% hold on
% plot(sol.time, sol.BG)
% title('BG')
% 
% subplot(2,2,4)
% hold on
% plot(sol.time, sol.BE)
% title('BE')
% fignum=3; 
% fourpanelplot(sol, fignum, FS)



% pause
end

% for i= [1 2 4]
%     tend = 120*60; 
%             Ti =3;
%             Fi = 3;
%             kpii= 3; 
%             kpgi=3;
%             kfii=3; 
%             kfgi=3;
%             kai=3;
%             
%             FibType = i; 
% 
%             
% %             if i == 4; 
% %                 abg=1; 
% %                 alpha = 1; 
% %                 beta = 1; 
% %             end 
% 
% sol = GeneralWeiselDriverTest(tend, 0.0001, 0.29412, fracgp(FibType),...
%         kacat, kbcat, kma, kmb, kpi(kpii), kpg(kpgi), kfi(kfii), kfg(kfgi), ka(kai), ThromScale,...
%         fpABcleavage, alpha(abg), beta(abg), gamma(abg), KPGSCALES(kind), KPISCALES(kind), RatesScale,TESTSCALE,scl,kpscalegpi); 
%     
%     
%  if i == 4 
%      output = sol;
%  else
%      output = 0; 
%  end
% % fignum=[3 4]; 
% % legnd = 'Thrombin';
% % fourpanelplot2alt(sol, fignum, FS,legnd) 
% % % fourpanelplot2alt2(sol, 21, FS,legnd) 
% time = sol.time/60;
% 
% figure(20)
% % subplot(1,2,1)
% hold on
% plot(time,  (sol.fab(1)- sol.fab)/sol.fab(1))
% title('fpa cleavage %')
% xlabel('time(min)')
% ylabel('% Fpa cleaved')
% % legend('2 (\gamma'' BS''s)/(E domain BS''s)','0.3 (\gamma'' BS''s)/(E domain BS''s)','0 (\gamma'' BS''s)/(E domain BS''s)')
% xlim([0 time(end)])
% ylim([0 1])
%  
% % subplot(1,2,2)
% hold on
% plot(time, 1-(sol.fab(1)- sol.f - sol.cf - sol.cfn - sol.cfr)/sol.fab(1))
% title('Fibrinopeptide cleavage \%')
% xlabel('time(min)')
% ylabel('% Fpb cleaved')
% legend('FPA - 2 (\gamma'' BS''s)/(E domain BS''s)','FPB - 2 (\gamma'' BS''s)/(E domain BS''s)'...
%     ,'FPA  - 1 (\gamma'' BS''s)/(E domain BS''s)','FPB - 1 ( \gamma'' BS''s)/(E domain BS''s)'...
%     ,'FPA - 0 (\gamma'' BS''s)/(E domain BS''s)','FPB - 0 (\gamma'' BS''s)/(E domain BS''s)')
% xlim([0 time(end)])
% ylim([0 1])
%   
% 'fibrinogen'
% sol.fab(end)
% end
%  


% 
% for i= [2 4]
%     tend = 30*60; 
%             Ti =3;
%             Fi = 3;
%             kpii= 3; 
%             kpgi=3;
%             kfii=3; 
%             kfgi=3;
%             kai=3;
%             
%             FibType = i; 
% 
%             
% %             if i == 4; 
% %                 abg=1; 
% %                 alpha = 1; 
% %                 beta = 1; 
% %             end 
% 
% sol = GeneralWeiselDriverTest(tend, 0.01, 5, fracgp(FibType),...
%         kacat, kbcat, kma, kmb, kpi(kpii), kpg(kpgi), kfi(kfii), kfg(kfgi), ka(kai), ThromScale,...
%         fpABcleavage, alpha(abg), beta(abg), gamma(abg), KPGSCALES(kind), KPISCALES(kind), RatesScale,TESTSCALE,scl, kpscalegpi); 
%     
%     
%  if i == 4 
%      output = sol;
%  else
%      output = 0; 
%  end
% % fignum=[3 4]; 
% % legnd = 'Thrombin';
% % fourpanelplot2alt(sol, fignum, FS,legnd) 
% % % fourpanelplot2alt2(sol, 21, FS,legnd) 
% time = sol.time/60;
% 
% figure(8)
% % subplot(1,2,1)
% hold on
% plot(time,  (sol.fab(1)- sol.fab)/sol.fab(1))
% title('fpa cleavage %')
% xlabel('time(min)')
% ylabel('% Fpa cleaved')
% % legend('2 (\gamma'' BS''s)/(E domain BS''s)','0.3 (\gamma'' BS''s)/(E domain BS''s)','0 (\gamma'' BS''s)/(E domain BS''s)')
% xlim([0 time(end)])
% ylim([0 1])
%  
% % subplot(1,2,2)
% hold on
% plot(time, 1-(sol.fab(1)- sol.f - sol.cf - sol.cfn - sol.cfr)/sol.fab(1))
% title('Kim-FPA/B cleavage \%')
% xlabel('time(min)')
% ylabel('% Fpb cleaved')
% legend('FPA  - 1 (\gamma'' BS''s)/(E domain BS''s)','FPB - 1 ( \gamma'' BS''s)/(E domain BS''s)'...
%     ,'FPA - 0 (\gamma'' BS''s)/(E domain BS''s)','FPB - 0 (\gamma'' BS''s)/(E domain BS''s)')
% xlim([0 time(end)])
% ylim([0 1])
%   
% 'fibrinogen'
% sol.fab(end)
% end
 
% 
% 

Fibrinogen = [2 4 6 8 10 14 18] ;
% F_max = zeros(7,2);
% F_clottime = zeros(7,2);
count= 0; 
for i= 1:length(Fibrinogen)
    for j = [2 4]
    tend = 15*60; 
            Ti =3;
            Fi = i;
            kpii= 3; 
            kpgi=3;
            kfii=3; 
            kfgi=3;
            kai=3;
            
            FibType = j; 

            
%             if i == 4; 
%                 abg=1; 
%                 alpha = 1; 
%                 beta = 1; 
%             end 
% 
sol = GeneralWeiselDriverTest(tend, 0.001, Fibrinogen(i), fracgp(FibType),...
        kacat, kbcat, kma, kmb, kpi(kpii), kpg(kpgi), kfi(kfii), kfg(kfgi), ka(kai), ThromScale,...
        fpABcleavage, alpha(abg), beta(abg), gamma(abg), KPGSCALES(kind), KPISCALES(kind), RatesScale,TESTSCALE,scl, kpscalegpi); 
    
%  solbatr = GeneralWeiselDriverTestBatr(tend, .0294, Fibrinogen(i), fracgp(FibType),...
%         kacat, kbcat, kma, kmb, kpi(kpii), kpg(kpgi), kfi(kfii), kfg(kfgi), ka(kai), ThromScale,...
%         fpABcleavage, alpha(abg), beta(abg), gamma(abg), KPGSCALES(kind), KPISCALES(kind), RatesScale,TESTSCALE,scl, kpscalegpi);     
  
% solbatr = GeneralWeiselDriverTestBatr(tend, .3, Fibrinogen(i), fracgp(FibType),...
%         kacat, kbcat, kma, kmb, kpi(kpii), kpg(kpgi), kfi(kfii), kfg(kfgi), ka(kai), ThromScale,...
%         fpABcleavage, 0, 0,0, KPGSCALES(kind), KPISCALES(kind), RatesScale,TESTSCALE,scl, kpscalegpi);     

  if j == 2 
     count = 1;
 else
     count = 2; 
 end
% % fignum=[3 4]; 
% % legnd = 'Thrombin';
% % fourpanelplot2alt(sol, fignum, FS,legnd) 
% % fourpanelplot2alt2(sol, 21, FS,legnd) 
F_clottime(i,count) = 0;

m= sol.m; 
time = sol.time; 
F_max = max(m)
% pause
% sol

for k = 1:(length(m)-5)
 if  m(k) == 0.5*F_max
     F_clottime(i,count)= time(k); 
 else if m(k) < 0.5*F_max & m(k+1) > 0.5*F_max
         F_clottime(i,count)= (time(k)+time(k+1))/2;
%      else
%          F_clottime= NaN;
     end
 end
% F_clottime(i,j) = 
end


% m= solbatr.m; 
% time = solbatr.time; 
% F_max = max(m)
% % pause
% % sol
% 
% for k = 1:(length(m)-5)
%  if  m(k) == 0.5*F_max
%      F_clottimeb(i,count)= time(k); 
%  else if m(k) < 0.5*F_max & m(k+1) > 0.5*F_max
%          F_clottimeb(i,count)= (time(k)+time(k+1))/2;
% %      else
% %          F_clottime= NaN;
%      end
%  end
% % F_clottime(i,j) = 
% end

  
% 'fibrinogen'
% sol.fab(end)
% 
% % if count == 1
% figure(27)
% subplot(1,2,1)
% hold on
% plot(sol.time,sol.m,'*')
% title('gp')
% xlabel('time')
% ylabel('Prot NUm')
% % else
% % figure(27)
% % subplot(1,2,2)
% hold on
% plot(sol.time,sol.m)
% title('ga')
% xlabel('time')
% ylabel('Prot NUm')

% end
    end
 end
% time = sol.time/60;

[F,CT] = KimData;
% [CTbt, Fbt] = kimdatascratch;
F_clottime/60
% pause
figure(7)
% figure
% subplot(1,2,1)
hold on
plot(Fibrinogen',  F_clottime/60)
plot(Fibrinogen', CT,'*')
% title('Thrombin activated')
xlabel('Fibrinogen')
ylabel('ClotTime')

 legend('\gamma_A/\gamma'' - model','\gamma_A/\gamma_A - model','\gamma_A/\gamma_A - data','\gamma_A/\gamma'' - data')
xlim([Fibrinogen(1) Fibrinogen(end)])
set(gca,'FontSize',28)

% subplot(1,2,2)
% hold on
% plot(Fibrinogen',  F_clottimeb/60)
% plot(Fibrinogen', CTbt,'*')
% title('Reptilase activated')
% xlabel('Fibrinogen')
% ylabel('ClotTime')
% 
%  legend('\gamma_A/\gamma'' - model','\gamma_A/\gamma_A - model','\gamma_A/\gamma_A - data','\gamma_A/\gamma'' - data')
% xlim([Fibrinogen(1) Fibrinogen(end)])
% set(gca,'FontSize',28)
% 



% Get coefficients of a line fit through the data.
coefficients = polyfit(Fibrinogen',F_clottime(:,1)./F_clottime(:,2) , 1);
% coefficientsb = polyfit(Fibrinogen',F_clottimeb(:,1)./F_clottimeb(:,2) , 1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit =Fibrinogen;% linspace(min(Fibrinogen), max(Fibrinogen), 1000);

% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);
% yFitb = polyval(coefficientsb , xFit);

y=F_clottime(:,1)./F_clottime(:,2);
% yb=F_clottimeb(:,1)./F_clottimeb(:,2);

[RatioBat, FbgnBat, RatioThr, FbgnThr] = KimRatios;

Fib = [2; 4; 6; 8; 10; 14; 18];
FibData = [Fib Fib];
FibRatio = [Fib Fib Fib];

figure(38)
% subplot(1,2,1)
hold on
plot(Fib,  y,'b')
% plot(Fib,  yb,'r')
errorbar(Fib, RatioThr(:,1), RatioThr(:,1)-RatioThr(:,2), RatioThr(:,3)-RatioThr(:,1),'b*')
% errorbar(Fib, RatioBat(:,1),RatioBat(:,1)- RatioBat(:,2), RatioBat(:,3)-RatioBat(:,1),'r*')
% plot(Fibrinogen',  RatioThr(:,1),'b*')
% plot(Fibrinogen',  RatioBat(:,1),'r*')
% plot(xFit,  yFit)
% plot(xFit,  yFitb
% plot(Fibrinogen', CT,'*')
xlim([0 20])
ylim([0 3])
title('clot time ratio')
xlabel('Fibrinogen')
ylabel('Ratio')
% legend('Thrombin','Batroxobin','Thrombin Data','Batroxabin Data')
% legend('Thrombin','Thrombin Data','Batroxabin Data')

% 'Changing with fibrinogen concentration'
% p_th_fbgn = anova1(F_clottime(:), FibData(:), 'off')
% pmeaning(p_th_fbgn,'fbgn_thr')
% % p_bat_fbgn = anova1(F_clottimeb(:),FibData(:) , 'off')
% % pmeaning(p_bat_fbgn,'fbgn_bat')
% 
% 'Changing with fibrin type'
% p_th_gp = anova1(F_clottime, {'gp','ga'}, 'off')
% pmeaning(p_th_gp,'gp_thr')
% % p_bat_gp = anova1(F_clottimeb, {'gp','ga'} , 'off')
% pmeaning(p_bat_gp,'gp_bat')

% 'Changing with fibrinogen'
% p_th_fbgn_ratio = anova1(y, Fib, 'off')
% pmeaning(p_th_fbgn,'fbgn_thr')
% p_bat_fbgn_ratio = anova1(y,Fib , 'off')
% pmeaning(p_bat_fbgn,'fbgn_bat')
% p_th_bat_fbgn = anova1([F_clottime(:); F_clottime(:)], [FibData(:); FibData(:)], 'off')
% pmeaning(p_th_bat_fbgn,'fbgn')
Fibrinogen = fibrin_mgml2muM(.5) ;
Thrombin = U2nmole(.1)*10^(-3) ;
% Bat = U2mumole_batr(5)
% TESTSCALE
% pause
% for i= 1:length(Fibrinogen)
%     for j = [2 4]
%     tend = 60*60; 
%             Ti =3;
%             Fi = 1;
%             kpii= 3; 
%             kpgi=3;
%             kfii=3; 
%             kfgi=3;
%             kai=3;
%             
%             FibType = j; 
% 
%             
% %             if i == 4; 
% %                 abg=1; 
% %                 alpha = 1; 
% %                 beta = 1; 
% %             end 
% % 100BU = 1 mg batrox
% % 5 BU  = 
% %
% TESTSCALE = 1; 
% 
% sol = GeneralWeiselDriverTest(tend, Thrombin, Fibrinogen, fracgp(FibType),...
%         kacat, kbcat, kma, kmb, kpi(kpii), kpg(kpgi), kfi(kfii), kfg(kfgi), ka(kai), ThromScale,...
%         fpABcleavage, alpha(abg), beta(abg), gamma(abg), KPGSCALES(kind), KPISCALES(kind), RatesScale,TESTSCALE,scl,kpscalegpi); 
% solbatr = GeneralWeiselDriverTestBatr(tend, Bat, Fibrinogen, fracgp(FibType),...
%         kacat, kbcat, kma, kmb, kpi(kpii), kpg(kpgi), kfi(kfii), kfg(kfgi), ka(kai), ThromScale,...
%         fpABcleavage, 0, 0, 0, KPGSCALES(kind), KPISCALES(kind), RatesScale,TESTSCALE,scl, kpscalegpi);
    
%     sol = GeneralWeiselDriverTest(tend, Thrombin, Fibrinogen, fracgp(FibType),...
%         kacat, kbcat, kma, kmb, kpi(kpii), kpg(kpgi), kfi(kfii), kfg(kfgi), ka(kai), ThromScale,...
%         fpABcleavage, alpha(abg), beta(abg), gamma(abg), KPGSCALES(kind), KPISCALES(kind), RatesScale,TESTSCALE,scl, 1);%kpscalegpi); 
% solbatr = GeneralWeiselDriverTestBatr(tend, Bat, Fibrinogen, fracgp(FibType),...
%         kacat, kbcat, kma, kmb, kpi(kpii), kpg(kpgi), kfi(kfii), kfg(kfgi), ka(kai), ThromScale,...
%         fpABcleavage, 0, 0, 0, KPGSCALES(kind), KPISCALES(kind), RatesScale,TESTSCALE,scl,1);% kpscalegpi,1);
%     
%     
% sol2 = GeneralWeiselDriverTest(tend, Thrombin, Fibrinogen, fracgp(FibType),...
%         kacat, kbcat, kma, kmb, kpi(kpii), kpg(kpgi), kfi(kfii), kfg(kfgi), ka(kai), ThromScale,...
%         fpABcleavage, alpha(abg), beta(abg), gamma(abg), KPGSCALES(kind), KPISCALES(kind), RatesScale,TESTSCALE,scl, 1); 
% %     
% solbatr = GeneralWeiselDriverTestBatr(tend, .026, Fibrinogen, fracgp(FibType),...
%         kacat, kbcat, kma, kmb, kpi(kpii), kpg(kpgi), kfi(kfii), kfg(kfgi), ka(kai), ThromScale,...
%         fpABcleavage, alpha(abg), beta(abg), gamma(abg), KPGSCALES(kind), KPISCALES(kind), RatesScale,TESTSCALE,scl, kpscalegpi);
    
    

    
%  figure(57)
% %  subplot(1,3,1)
%  hold on
%  plot(sol.time,sol.m)
%  
% %  plot(sol2.time,sol2.m)
%  
% %  subplot(1,3,2)
% %  hold on
%  plot(solbatr.time,solbatr.m)
%  
% %   subplot(1,3,3)
% %  hold on
% %  plot(sol2.time,sol2.m)
%     
%     end
    
%     figure(57)
% %  subplot(1,3,1)
%  set(gca,'FontSize',28)
%  xlabel('time')
% ylabel('Prot Num')
% title('Thrombin')
% legend('Thrombin -\gamma_A/\gamma''','Batroxobin-\gamma_A/\gamma''','Thrombin-\gamma_A/\gamma_A','Batroxobin-\gamma_A/\gamma_A')
% legend('Thrombin -\gamma''','Batr-\gamma''','No steric -\gamma''','Thrombin-\gamma_A','Batr-\gamma_A','No steric -\gamma_A')
%  subplot(1,3,2)
%  set(gca,'FontSize',28)
%  ylabel('Prot Num')
%  xlabel('time')
%  title('Batroxobin')
%  legend('\gamma''','\gamma_A')
%  
%   subplot(1,3,3)
%  set(gca,'FontSize',28)
%  ylabel('Prot Num')
%  xlabel('time')
%  title('No Steric Inhbition')
% legend('\gamma''','\gamma_A')
% % end
  
 end
% ylim([0 1])
 
% subplot(1,2,2)
% hold on
% plot(time, 1-(sol.fab(1)- sol.f - sol.cf - sol.cfn - sol.cfr)/sol.fab(1))
% title('Kim-FPA/B cleavage \%')
% xlabel('time(min)')
% ylabel('% Fpb cleaved')
% legend('FPA  - 1 (\gamma'' BS''s)/(E domain BS''s)','FPB - 1 ( \gamma'' BS''s)/(E domain BS''s)'...
%     ,'FPA - 0 (\gamma'' BS''s)/(E domain BS''s)','FPB - 0 (\gamma'' BS''s)/(E domain BS''s)')
% xlim([0 time(end)])
% ylim([0 1])
