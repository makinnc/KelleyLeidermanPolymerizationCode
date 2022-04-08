function [sol]= weiselmodelfcnnew(tend,names_IC,params_IC,names_rates,params_rates,scales)

rates.TESTSCAlE=0; 
rates.ThromScale=0; %1 is thrombin is kicked into fluid by fpb cleavage, 0 if not
rates.gamma=1; 
rates.alpha=1;
rates.beta=1; 
rates.scl = 0.5; 
rates.kpscalepgi=0.5; 

%All rates that can be changed as inputs
rates.ka = [0];% 10 1 .1]*10^(-2); %s^-1 [1 1 1 1];%
rates.kpi = 4*10^(-18)*(6.022*10^23*10^(-6)); %original number in L/(molecules*s)- put in 1/muM s (L/mumoles*s)- 
rates.kpg = 10^(-16)*(6.022*10^23*10^(-6)); %original number in L/(molecules*s)- put in 1/muM s (L/mumoles*s) - protofibrils grow in length
rates.kfi = 10^(-21)*(6.022*10^23*10^(-6)); %original number in L/(molecules*s)- put in 1/muM s (L/mumoles*s) - protofibrils aggregate into fiber
rates.kfg = 2*10^(-17)*(6.022*10^23*10^(-6)); %original number in L/(molecules*s)- put in 1/muM s (L/mumoles*s) - protofibrils add to fiber

rates.fracgp=[.3];
% rates.pf=[0, 30,60,90,120,150,180,210,240];
% rates.diameter=[54, 75.5, 92, 105.5, 118, 128.7,139.1,148.7];
% rates.numb_innac=[0, 4,17,35,54,75,96,120,141];
% rates.frac_innac=numb_innac./pf;
% rates.frac_innac(1)=0;
rates.kpgscale=1;
rates.kpiscale=1; 

rates.kacat=84;
rates.kma=7.2;
rates.kbcat=49; %7.4;
rates.kmb=7.5;




rates.kone=1;
rates.koffe=2.8;
rates.kong=1;
rates.koffg=9;
rates.kone2=1500;
% rates.koffe2=rates.koffe;
rates.kong2=3000;
% rates.koffg2= rates.koffg;

ICs.T_IC=[.100];
ICs.fab_IC = 3;











ICs.fb_IC=0;
ICs.Efab_IC=0;
ICs.Efb_IC=0;




ICs.f2_IC = 0; 
ICs.f3_IC = 0; 
ICs.f4_IC = 0; 
ICs.f5_IC = 0; 
ICs.f6_IC = 0; 
ICs.f7_IC = 0; 
ICs.f8_IC = 0; 
ICs.f9_IC = 0; 
ICs.f10_IC = 0; 
ICs.fn_IC = 0; 
ICs.fr_IC = 0; 
ICs.ftotn_IC=0;
ICs.cfb_IC=0;
ICs.cfn_IC=0;
ICs.cfr_IC=0;

ICs.cb_IC=0; 
ICs.cf_IC=0; 


ICfrac=1;


ICs.f_IC = ICs.fab_IC*(1-ICfrac); 
ICs.fab_IC=ICs.fab_IC*(ICfrac);


%%%%%%%%%%%%%%%%%%

% names_rates{1}
% names_IC{1}
% 
% 
% rates.fracgp
% rates.(names_rates{1})=.4

for i=1:length(names_rates) 
%     params_rates.(names_rates{i})
    rates.(names_rates{i})=params_rates.(names_rates{i});
end

for i=1:length(names_IC)   
    ICs.(names_IC{i})=params_IC.(names_IC{i}); 
end

%  totalpolym= ICs.cfn_IC+ICs.cfr_IC+2*ICs.f2_IC+3*ICs.f3_IC+4*ICs.f4_IC+5*ICs.f5_IC ...
%        +6*ICs.f6_IC+7*ICs.f7_IC+8*ICs.f8_IC+9*ICs.f9_IC+10*ICs.f10_IC;
% totalnotpolym= ICs.f_IC+ICs.fab_IC+ICs.fb_IC+ICs.Efab_IC+ICs.Efb_IC; %fa_IC+f_IC;
    
   fracgp=rates.fracgp;     
% ICs.bsE_IC=2*totalnotpolym+1.7*totalpolym;
% bsE_IC=1.7*totalpolym;
% bsEp_IC=0.3*totalpolym ;
% bsG_IC=0.3*totalpolym;
% bsGp_IC=0.3*totalnotpolym;
% bsE_IC=(1-fracgp)*totalpolym;
% bsEp_IC=fracgp*totalpolym ;
% bsG_IC=fracgp*totalpolym;
% bsGp_IC=fracgp*totalnotpolym;

% pause

ICs.E1_IC = 0;
ICs.BE1_IC = 0; 
ICs.E_IC =0;%0;% 0.3*fib; 
ICs.BE_IC = 0; 

ICs.G1_IC = ICs.fab_IC*fracgp;%.3*fib; 
ICs.BG1_IC = 0; 


ICs.G_IC = 0;%.3*fib; 
ICs.BG_IC = 0; 

ICs.B_IC = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%All rates that depend on other rates. 
if rates.kacat == 0  
    
     
rates.kam=0;%1;
rates.kap=0;
rates.kbm=0;%1;
rates.kbp=0; 

    
else
 
     
rates.kam=rates.koffe;%1;
rates.kap=(rates.kam+rates.kacat)/rates.kma;
rates.kbm=rates.koffe;%1;
rates.kbp=(rates.kbm+rates.kbcat)/rates.kmb; 

 

end

kpiscale=rates.kpiscale; 
kpgscale=rates.kpgscale; 


rates.kpi21 =kpiscale*OligScale( scales, 3) * rates.kpi;
rates.kpi22 =kpiscale*OligScale( scales, 4) * rates.kpi;
rates.kpi23 =kpiscale*OligScale( scales, 5) * rates.kpi;
rates.kpi24 =kpiscale*OligScale( scales, 6) * rates.kpi;
rates.kpi25 =kpiscale*OligScale( scales, 7) * rates.kpi;
rates.kpi26 =kpiscale*OligScale( scales, 8) * rates.kpi;
rates.kpi27 =kpiscale*OligScale( scales, 9) * rates.kpi;
rates.kpi28 =kpiscale*OligScale( scales, 10) * rates.kpi;
rates.kpi29 =kpiscale*OligScale( scales, 11) * rates.kpi;
rates.kpi31 =kpiscale*OligScale( scales, 4) * rates.kpi;
rates.kpi33 =kpiscale*OligScale( scales, 6) * rates.kpi;
rates.kpi34 =kpiscale*OligScale( scales, 7) * rates.kpi; 
rates.kpi35 =kpiscale*OligScale( scales, 8) * rates.kpi;
rates.kpi36 =kpiscale*OligScale( scales, 9) * rates.kpi;
rates.kpi37 =kpiscale*OligScale( scales, 10) * rates.kpi;
rates.kpi38 =kpiscale*OligScale( scales, 11) * rates.kpi;
rates.kpi41 =kpiscale*OligScale( scales, 5) * rates.kpi;
rates.kpi44 =kpiscale*OligScale( scales, 8) * rates.kpi;
rates.kpi45 =kpiscale*OligScale( scales, 9) * rates.kpi;
rates.kpi46 =kpiscale*OligScale( scales, 10) * rates.kpi;
rates.kpi47 =kpiscale*OligScale( scales, 11) * rates.kpi;
rates.kpi51 =kpiscale*OligScale( scales, 6) * rates.kpi;
rates.kpi55 =kpiscale*OligScale( scales, 10) * rates.kpi;
rates.kpi56 =kpiscale*OligScale( scales, 11) * rates.kpi;
rates.kpi61 =kpiscale*OligScale( scales, 7) * rates.kpi;
rates.kpi71 =kpiscale*OligScale( scales, 8) * rates.kpi;
rates.kpi81 =kpiscale*OligScale( scales, 9) * rates.kpi;
rates.kpi91 =kpiscale*OligScale( scales, 10) * rates.kpi;
rates.kpi101 =kpiscale*OligScale( scales, 11) * rates.kpi;






rates.kpi210 = kpiscale*OligScale( scales, 12)*rates.kpi; 
rates.kpi39 = kpiscale*OligScale( scales, 11)*rates.kpi;  
rates.kpi310 = kpiscale*OligScale( scales, 13)*rates.kpi; 
rates.kpi48 = kpiscale*OligScale( scales, 12)*rates.kpi;  
rates.kpi49 = kpiscale*OligScale( scales, 13)*rates.kpi; 
rates.kpi410 = kpiscale*OligScale( scales, 14)*rates.kpi;  
rates.kpi57 = kpiscale*OligScale( scales, 12)*rates.kpi; 
rates.kpi58 = kpiscale*OligScale( scales, 13)*rates.kpi;  
rates.kpi59 = kpiscale*OligScale( scales, 14)*rates.kpi; 
rates.kpi510 = kpiscale*OligScale( scales, 15)*rates.kpi; 
rates.kpi66 = kpiscale*OligScale( scales, 12)*rates.kpi; 
rates.kpi67 = kpiscale*OligScale( scales, 13)*rates.kpi; 
rates.kpi68 = kpiscale*OligScale( scales, 14)*rates.kpi;  
rates.kpi69 = kpiscale*OligScale( scales, 15)*rates.kpi;  
rates.kpi610 = kpiscale*OligScale( scales, 16)*rates.kpi;  
rates.kpi77 = kpiscale*OligScale( scales, 14)*rates.kpi; 
rates.kpi78 = kpiscale*OligScale( scales, 15)*rates.kpi; 
rates.kpi79 = kpiscale*OligScale( scales, 16)*rates.kpi; 
rates.kpi710 = kpiscale*OligScale( scales, 17)*rates.kpi; 
rates.kpi88 = kpiscale*OligScale( scales, 16)*rates.kpi;  
rates.kpi89 = kpiscale*OligScale( scales, 17)*rates.kpi; 
rates.kpi810 = kpiscale*OligScale( scales, 18)*rates.kpi; 
rates.kpi99 = kpiscale*OligScale( scales, 18)*rates.kpi; 
rates.kpi910 = kpiscale*OligScale( scales, 19)*rates.kpi; 
rates.kpi1010 = kpiscale*OligScale( scales, 20)*rates.kpi; 

rates.kpg2 = kpgscale*ProtScale( scales, 2)*rates.kpg;
rates.kpg3 = kpgscale*ProtScale( scales, 3)*rates.kpg;
rates.kpg4 = kpgscale*ProtScale( scales, 4)*rates.kpg;
rates.kpg5 = kpgscale*ProtScale( scales, 5)*rates.kpg;
rates.kpg6 = kpgscale*ProtScale( scales, 6)*rates.kpg;
rates.kpg7 = kpgscale*ProtScale( scales, 7)*rates.kpg;
rates.kpg8 = kpgscale*ProtScale( scales, 8)*rates.kpg;
rates.kpg9 = kpgscale*ProtScale( scales, 9)*rates.kpg;
rates.kpg10 = kpgscale*ProtScale( scales, 10)*rates.kpg;

% rates.kpi210 = 0; 
% rates.kpi39 = 0;  
% rates.kpi310 = 0; 
% rates.kpi48 = 0;  
% rates.kpi49 = 0; 
% rates.kpi410 = 0;  
% rates.kpi57 = 0; 
% rates.kpi58 = 0;  
% rates.kpi59 = 0; 
% rates.kpi510 = 0; 
% rates.kpi66 = 0; 
% rates.kpi67 = 0; 
% rates.kpi68 = 0;  
% rates.kpi69 = 0;  
% rates.kpi610 = 0;  
% rates.kpi77 = 0; 
% rates.kpi78 = 0; 
% rates.kpi79 = 0; 
% rates.kpi710 = 0; 
% rates.kpi88 = 0;  
% rates.kpi89 = 0; 
% rates.kpi810 = 0; 
% rates.kpi99 = 0; 
% rates.kpi910 = 0; 
% rates.kpi1010 = 0; 
% 
% rates.kpg2 = 0;
% rates.kpg3 = 0;
% rates.kpg4 = 0;
% rates.kpg5 = 0;
% rates.kpg6 = 0;
% rates.kpg7 = 0;
% rates.kpg8 = 0;
% rates.kpg9 = 0;
% rates.kpg10 = 0;


% rates.kpi210, rates.kpi39,  rates.kpi310, rates.kpi48,  rates.kpi49, rates.kpi410, ...
% rates.kpi57, rates.kpi58,  rates.kpi59, rates.kpi510, rates.kpi66, rates.kpi67, ...
% rates.kpi68,  rates.kpi69,  rates.kpi610,  rates.kpi77, rates.kpi78, rates.kpi79, rates.kpi710, ...
% rates.kpi88,  rates.kpi89, rates.kpi810, rates.kpi99, rates.kpi910, rates.kpi1010... 
% 
% 
% rates.kpg2, rates.kpg3, rates.kpg4, rates.kpg5, rates.kpg6, rates.kpg7, rates.kpg8, rates.kpg9, rates.kpg10, 

rates.koffe2=rates.koffe;
rates.koffg2= rates.koffg;


if fracgp ==0 
    rates.koffg=0;
    rates.kong=0; 
    rates.kong2=0; 
    rates.kone2=0; 
    rates.koffg2=0;
    rates.koffe2=0;
    
end

%%%%%%%%%%%%%%%%%%%%%%%
pw = [ rates.ka, rates.kpi, rates.kpi21, rates.kpi22, rates.kpi23, rates.kpi24, rates.kpi25, rates.kpi26, rates.kpi27, rates.kpi28, ...
     rates.kpi29, rates.kpi31, rates.kpi33, rates.kpi34, rates.kpi35, rates.kpi36, rates.kpi37, rates.kpi38, rates.kpi41, rates.kpi44, ...
     rates.kpi45, rates.kpi46, rates.kpi47, rates.kpi51, rates.kpi55, rates.kpi56, rates.kpi61, rates.kpi71, rates.kpi81,rates. kpi91,...
     rates.kpi101, rates.kpg, rates.kfi, rates.kfg...
     rates.kpi210, rates.kpi39,  rates.kpi310, rates.kpi48,  rates.kpi49, rates.kpi410, ...
     rates.kpi57, rates.kpi58,  rates.kpi59, rates.kpi510, rates.kpi66, rates.kpi67, ...
     rates.kpi68,  rates.kpi69,  rates.kpi610,  rates.kpi77, rates.kpi78, rates.kpi79, rates.kpi710, ...
     rates.kpi88,  rates.kpi89, rates.kpi810, rates.kpi99, rates.kpi910, rates.kpi1010...
     rates.kpg2, rates.kpg3, rates.kpg4, rates.kpg5, rates.kpg6, rates.kpg7, rates.kpg8, rates.kpg9, rates.kpg10];

pe=[rates.kacat,rates.kap,rates.kam,rates.kbcat,rates.kbp,rates.kbm];
pb=[rates.kone,rates.koffe,rates.kong,rates.koffg,rates.kone2,rates.koffe2,rates.kong2,rates.koffg2, rates.TESTSCALE,rates.scl,rates.kpscalegpi];


ICs.G1_IC;
init_cond = [ICs.fab_IC,ICs.fb_IC,ICs.Efab_IC,ICs.Efb_IC, ICs.f_IC, ICs.f2_IC, ICs.f3_IC, ICs.f4_IC, ICs.f5_IC, ICs.f6_IC, ICs.f7_IC,...
    ICs.f8_IC, ICs.f9_IC, ICs.f10_IC, ICs.fn_IC, ICs.fr_IC,ICs.ftotn_IC,ICs.cfb_IC, ICs.cfn_IC, ICs.cfr_IC,...
    ICs.T_IC,ICs.E1_IC,ICs.BE1_IC,ICs.E_IC,ICs.BE_IC,ICs.G1_IC,ICs.BG1_IC,ICs.G_IC,ICs.BG_IC,ICs.B_IC,ICs.cb_IC,ICs.cf_IC];
%  pause

% pause;
% rates.kpi
% rates.kbcat
% rates.kacat

 t_start = 0; % Default start time is 0
 t_final = tend; % Default final time is 1
   
% rates.fracgp

% rates.ThromScale=0;

rates.gamma;
[sol]=RHSweiselFPBtest(t_final,t_start,init_cond,pw,pb,pe,rates.fracgp, [rates.alpha, rates.beta, rates.gamma],rates.ThromScale);
sol.init_cond=init_cond; 
sol.alpha = rates.alpha; 
sol.beta = rates.beta;
sol.gamma = rates.gamma;



sol.kpg=rates.kpg; 
sol.kfg=rates.kfg; 

T= sol.T + sol.Efab + sol.Efb + sol.BE1 + sol.BE + sol.BG + sol.BG1 + sol.B;
F= (sol.fab +sol.Efab + sol.Efb + sol.fb) + sol.f  ...
    + 2*sol.f2 + 3*sol.f3 + 4*sol.f4 + 5*sol.f5 + 6*sol.f6 + 7*sol.f7 ...
    + 8*sol.f8 + 9*sol.f9 + 10*sol.f10 ...
    + sol.cfb + sol.cfn + sol.cfr;
% F= 1/sol.whole*(sol.fab +sol.Efab + sol.Efb + sol.fb) + sol.f  ...
%     + 2*sol.f2 + 3*sol.f3 + 4*sol.f4 + 5*sol.f5 + 6*sol.f6 + 7*sol.f7 ...
%     + 8*sol.f8 + 9*sol.f9 + 10*sol.f10 ...
%     + sol.cfb + sol.cfn + sol.cfr;

 
%  max(abs(ICs.T_IC(1)-T))
%  
%  max(abs(sol.fab(1)-F))
 
%  ICs.fab_IC(1)
%  sol.whole
% max(F)
% min(F)

end


