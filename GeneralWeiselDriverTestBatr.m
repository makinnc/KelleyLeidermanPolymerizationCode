



function sol = GeneralWeiselDriverTest(tend, Thrombin, Fibrinogen, fracgp, kacat, kbcat, kma, kmb, kpi, kpg, kfi, kfg, ka, ThromScale, fpABcleavage, alpha, beta, gamma, KPGSCALES, KPISCALES, RatesScale, TESTSCALE,scl,kpscalegpi)
FS=32;

% scl
%T0 turn off enzymatic cleavage of thrombin and make it kfa*fa set kbcat = kacat = 0
%set kfa too

% kbcat = 0; 
%  kacat = 0; 
%  ka = .1; 
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
        'ka', ka; 
        'kbcat',kbcat;
    'kacat',kacat;
    'kmb',kmb;
    'kma',kma;
    'kpi', kpi ;%scalekpi*4*10^(-18)*(6.022*10^23*10^(-6));
    'kfi', kfi ;%scalekpi*4*10^(-18)*(6.022*10^23*10^(-6));
    'kfg',kfg;
    'kpg', kpg;
    'alpha',alpha;
    'beta',beta;
    'gamma',gamma; 
    'kpiscale', KPISCALES;
    'kpgscale', KPGSCALES;
    'ThromScale', ThromScale;
    'TESTSCALE', TESTSCALE;
    'scl', scl;
    'kpscalegpi', kpscalegpi;};
[nrates,dummy]= size(rates);




sol=PolymModel(tend,nICs,nrates,Init_conds,rates, RatesScale);


%  testplot(sol)

T= sol.T + sol.Efab + sol.Efb + sol.BE1 + sol.BE + sol.BG + sol.BG1 + sol.B;
F= (sol.fab +sol.Efab + sol.Efb + sol.fb) + sol.f  ...
    + 2*sol.f2 + 3*sol.f3 + 4*sol.f4 + 5*sol.f5 + 6*sol.f6 + 7*sol.f7 ...
    + 8*sol.f8 + 9*sol.f9 + 10*sol.f10 ...
    + sol.cfb + sol.cfn + sol.cfr;




 
%  max(abs(sol.T(1)-T))
%  
%  max(abs(sol.fab(1)-F))
%  figure(2)

end
 
 
 function testplot(sol)
 FS=32;
 
 T= sol.T + sol.Efab + sol.Efb + sol.BE1 + sol.BE + sol.BG + sol.BG1 + sol.B;
F= (sol.fab +sol.Efab + sol.Efb + sol.fb) + sol.f  ...
    + 2*sol.f2 + 3*sol.f3 + 4*sol.f4 + 5*sol.f5 + 6*sol.f6 + 7*sol.f7 ...
    + 8*sol.f8 + 9*sol.f9 + 10*sol.f10 ...
    + sol.cfb + sol.cfn + sol.cfr;

F2 = sol.cb+sol.cfb+sol.cf+sol.cfn+sol.cfr+(sol.fab+sol.Efab+sol.Efb+sol.fb)+sol.f;

Ffree = sol.fab + sol.fb +sol.Efab + sol.Efb + sol.f; 

Folig=  + 2*sol.f2 + 3*sol.f3 + 4*sol.f4 + 5*sol.f5 + 6*sol.f6 + 7*sol.f7 ...
    + 8*sol.f8 + 9*sol.f9 + 10*sol.f10;



Folig2= sol.cb + sol.cf; 

t=sol.time; 
tend= t(end); 
% plot(t,2*sol.f2,t,3*sol.f3,t,4*sol.f4,t,5*sol.f5,t,6*sol.f6,t,7*sol.f7,t,8*sol.f8,t,9*sol.f9,t,10*sol.f10)
% % plot(t,sol.f2,t,sol.f3,t,sol.f4,t,sol.f5,t,sol.f6,t,sol.f7,t,sol.f8,t,sol.f9,t,sol.f10)
% hold on
% plot(t,Folig,t,Folig2)

figure(1)
xl=tend;
subplot(1,2,1)
hold on
plot(t,T,t,T)
% legend('all','cfb/cfn')
title('Thrombin')
xlim([0 xl])
ylim ([ 0 2*sol.T(1)]);

% subplot(1,3,2)
% hold on
% plot(t,Folig,t,Folig2)
% xlim([0 xl])
% title('Should be the same')
% legend('all','cfb/cfn')

subplot(1,2,2)
hold on
plot(t,F,t,F2)
% legend('all','cfb/cfn')
title('fibrin')
xlim([0 xl])
ylim ([ 0 1.5*sol.fab(1)]);

ind1= find(max(F) == F); 


% plot(t,2*sol.f2,t,3*sol.f3,t,4*sol.f4,t,5*sol.f5,t,6*sol.f6,t,7*sol.f7,t,8*sol.f8,t,9*sol.f9,t,10*sol.f10)% %
% plot(t,sol.f2,t,sol.f3,t,sol.f4,t,sol.f5,t,sol.f6,t,sol.f7,t,sol.f8,t,sol.f9,t,sol.f10)

figure(2)

subplot(3,2,1)
hold on
plot(t,sol.cb)
xlim([0 xl])
ylim([0 max(sol.cb)])
legend('cb')

subplot(3,2,2)
hold on
plot(t,sol.cf)
xlim([0 xl])
% ylim([0 max(sol.cf)])
legend('cf')


subplot(3,2,3)
hold on
plot(t,sol.cfb)
xlim([0 xl])
ylim([0 sol.fab(1)])
legend('cfb')

subplot(3,2,4)
hold on
plot(t,sol.cfn)
xlim([0 xl])
ylim([0 sol.fab(1)])
legend('cfn')

subplot(3,2,5)
hold on
plot(t,sol.cfr)
xlim([0 xl])
ylim([0 sol.fab(1)])
legend('cfr')

subplot(3,2,6)
hold on
plot(t,sol.fab)
xlim([0 xl])
ylim([0 sol.fab(1)])
legend('fibrinogen')

% figure(3)
% fignum=3; 
% fourpanelplot2(sol,fignum,FS)
 
 end

function sol=PolymModel(tend,nICs,nrates,Init_conds,rates,scales)


for i=1:nICs
    
    names_IC{i}=Init_conds{i,1};
     params_IC.(names_IC{i})=Init_conds{i,2};
end

for i=1:nrates
    
   names_rates{i}=rates{i,1};
     params_rates.(names_rates{i})=rates{i,2};
end


[sol]= weiselmodelfcnnewbatr(tend,names_IC,params_IC,names_rates,params_rates, scales);

end

% pf=[0, 30,60,90,120,150,180,210,240];
% diameter=[54, 75.5, 92, 105.5, 118, 128.7,139.1,148.7];
% numb_innac=[0, 4,17,35,54,75,96,120,141];
% frac_innac=numb_innac./pf;
% frac_innac(1)=0;