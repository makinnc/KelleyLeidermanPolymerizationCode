function [sol]=RHSweiselFPBtest(t_final,t_start,init_cond,pw,pb,pe,fracgp,fpbscales,ThromScale)
% Solves a system of ODEs from t=t_start to t=t_final 
% If no start time is given, then t_start = 0 
% If no start or final time is given, then t_start = 0, t_final = 1 
%
%
% This file was created by issuing command: 
% 	python make_ode_m_file.py reactionstest.csv reactionstest.m
%



%      t_start = 0; % Default start time is 0
%      t_final = 1; % Default final time is 1


% fracgp

options = odeset('RelTol',2.22045e-11,'AbsTol',1e-13,'NonNegative', [1:32]);

fpbscales;
%------------------------------ Main Solve ---------------------------%
[time,y] = ode15s(@(t,y)RHS(t,y,pw,pb,pe,fracgp,fpbscales,ThromScale), [t_start:.1:t_final], init_cond, options);
%---------------------------------------------------------------------%

fab = y(:,1); 
fb = y(:,2);
Efab=y(:,3);
Efb=y(:,4);

% f2+f -> f3
% f2+fb-> f3
% fi+fb-> fi+1
% fb+fb =>f2
% 
% fracfb= fb/(fb+f)
% 
% fn*(1-fracfb) -> fr
f = y(:,5); 
f2 = y(:,6); 
f3 = y(:,7); 
f4 = y(:,8); 
f5 = y(:,9); 
f6 = y(:,10); 
f7 = y(:,11); 
f8 = y(:,12); 
f9 = y(:,13); 
f10 = y(:,14); 
fn = y(:,15); 
fr = y(:,16); 
ftotn=y(:,17);
cfb = y(:,18);
cfn=y(:,19);
cfr=y(:,20);
T = y(:,21); 
E1 = y(:,22); 
BE1 = y(:,23); 
E=y(:,24);
BE=y(:,25);
G1=y(:,26);
BG1=y(:,27);
G=y(:,28);
BG=y(:,29);
B=y(:,30);

cb=y(:,31);
cf=y(:,32);




% average number fibrin/protofibril
nalt=(cfn+cfb)./fn;
% average number of protofibril/fiber
malt=ftotn./fr;
% average length of fibers
lalt=cfr./ftotn;

% average number fibrin/protofibril
% ind = find(cfb < 10^(-15) & cfn < 10^(-15) & fn < 10^(-15)); 
n=(cfb+cfn)./fn;
% n(ind)=0; 
% average number of protofibril/fiber
m=ftotn./fr;
% average length of fibers
l=cfr./ftotn;

% n(1)=0;
% m(1)=0;
% l(1)=0;
sol.fb=fb;
sol.Efab=Efab;
sol.Efb=Efb;

sol.fab=fab;
sol.f=f;
sol.f2=f2;
sol.f3=f3;
sol.f4=f4;
sol.f5=f5;
sol.f6=f6;
sol.f7=f7;
sol.f8=f8;
sol.f9=f9;
sol.f10=f10;
sol.fn=fn;
sol.fr=fr;
sol.ftotn=ftotn;
sol.cfb= cfb; 
sol.cfn=cfn;
sol.cfr=cfr;
sol.n=n;
sol.m=m;
sol.l=l;


sol.T=T;
sol.E1=E1;
sol.BE1=BE1;
sol.E=E;
sol.BE=BE;
sol.G1=G1;
sol.BG1=BG1;
sol.G=G;
sol.BG=BG;
sol.B=B;

sol.cb=cb; 
sol.cf=cf; 

sol.time=time;

sol.kacat = pe(1); 
sol.kbcat = pe(4); 
sol.alpha = fpbscales(1);
sol.beta  = fpbscales(2);
sol.gamma = fpbscales(3); 

sol.pe = pe; 

% sol.whole=fpbscales(4); 
end




function dy = RHS(t,y,pw,pb,pe,fracgp,fpbscales,ThromScale)

% y(y<10^(-16))= 0 ;

dy = zeros(32,1);

alpha=fpbscales(1); 
beta=fpbscales(2); 
gamma=fpbscales(3);

% whole=fpbscales(4);
% gamma=1; 

% y=max(0,y); 
% Rename Variables 

fab   = y(1); 
fb = y(2);
Efab = y(3);
Efb = y (4);
f   = y(5); 
f2   = y(6); 
f3   = y(7); 
f4   = y(8); 
f5   = y(9); 
f6   = y(10); 
f7   = y(11); 
f8   = y(12); 
f9   = y(13); 
f10   = y(14); 
fn   = y(15); 
fr   = y(16); 
ftotn = y(17);
cfb = y(18);
cfn = y(19);
cfr = y(20);
T   = y(21); 
E1   = y(22); 
BE1   = y(23); 
E   = y(24); 
BE   = y(25); 
G1   = y(26); 
BG1   = y(27); 
G   = y(28); 
BG   = y(29);
B   = y(30); 

cb = y(31); 
cf = y(32); 




scl = pb(10);
Sgpi = pb(11);

[Pgg, Pgp, Pgpgp]= ProbGamP(fracgp/2);

% Scale inhibiting fibrin end to end binding multiplies kpi, kpg and all
% variants.
kpscalegpi = Pgg + Sgpi*(Pgp +Pgpgp);


% Rename Kinetic Parameters 
ka = kpscalegpi * pw(1);  
kpi = kpscalegpi * pw(2);  
kpi21 = kpscalegpi * pw(3);  
kpi22 = kpscalegpi * pw(4);  
kpi23 = kpscalegpi * pw(5);  
kpi24 = kpscalegpi * pw(6);  
kpi25 = kpscalegpi * pw(7);  
kpi26 = kpscalegpi * pw(8);  
kpi27 = kpscalegpi * pw(9);  
kpi28 = kpscalegpi * pw(10);  
kpi29 = kpscalegpi * pw(11);  
kpi31 = kpscalegpi * pw(12);  
kpi33 = kpscalegpi * pw(13);  
kpi34 = kpscalegpi * pw(14);  
kpi35 = kpscalegpi * pw(15);  
kpi36 = kpscalegpi * pw(16);  
kpi37 = kpscalegpi * pw(17);  
kpi38 = kpscalegpi * pw(18);  
kpi41 = kpscalegpi * pw(19);  
kpi44 = kpscalegpi * pw(20);  
kpi45 = kpscalegpi * pw(21);  
kpi46 = kpscalegpi * pw(22);  
kpi47 = kpscalegpi * pw(23);  
kpi51 = kpscalegpi * pw(24);  
kpi55 = kpscalegpi * pw(25);  
kpi56 = kpscalegpi * pw(26);  
kpi61 = kpscalegpi * pw(27);  
kpi71 = kpscalegpi * pw(28);  
kpi81 = kpscalegpi * pw(29);  
kpi91 = kpscalegpi * pw(30);  
kpi101 = kpscalegpi * pw(31);  
kpg = kpscalegpi * pw(32);  
kfi = kpscalegpi * pw(33);  
kfg = kpscalegpi * pw(34); 

kpi210 = kpscalegpi * pw(35) ;
kpi39 = kpscalegpi * pw(36); 
kpi310 = kpscalegpi * pw(37); 
kpi48 = kpscalegpi * pw(38); 
kpi49 = kpscalegpi * pw(39); 
kpi410 = kpscalegpi * pw(40); 
kpi57 = kpscalegpi * pw(41); 
kpi58 = kpscalegpi * pw(42); 
kpi59 = kpscalegpi * pw(43); 
kpi510 = kpscalegpi * pw(44);
kpi66 = kpscalegpi * pw(45);
kpi67 = kpscalegpi * pw(46); 
kpi68 = kpscalegpi * pw(47); 
kpi69 = kpscalegpi * pw(48); 
kpi610 = kpscalegpi * pw(49); 
kpi77 = kpscalegpi * pw(50); 
kpi78 = kpscalegpi * pw(51); 
kpi79 = kpscalegpi * pw(52); 
kpi710 = kpscalegpi * pw(53);
kpi88 = kpscalegpi * pw(54); 
kpi89 = kpscalegpi * pw(55); 
kpi810 = kpscalegpi * pw(56);
kpi99 = kpscalegpi * pw(57); 
kpi910 = kpscalegpi * pw(58); 
kpi1010 = kpscalegpi * pw(59);

kpg2 = kpscalegpi * pw(60);
kpg3 = kpscalegpi * pw(61);
kpg4 = kpscalegpi * pw(62);
kpg5 = kpscalegpi * pw(63);
kpg6 = kpscalegpi * pw(64);
kpg7 = kpscalegpi * pw(65);
kpg8 = kpscalegpi * pw(66);
kpg9 = kpscalegpi * pw(67);
kpg10 = kpscalegpi * pw(68);


% kpi56=0;
% kpi47=0;
% kpi38=0;
% kpi29=0;


% kpi22=0;
% kpi23=0;
% kpi24=0;
% kpi25=0;
% kpi26=0;
% kpi27=0;
% kpi28=0;


% kpi33=0;
% kpi34=0;
% kpi35=0;
% kpi36=0;
% kpi37=0;


% kpi44=0;
% kpi45=0;
% kpi46=0;



% kpi55=0;





% pe=[rates.kacat,rates.kap,rates.kam,rates.kbcat,rates.kbp,rates.kbm];

% pe=[kacat,kap,kam,kbcat,kbp,kbm];
kacat=pe(1);
kap=pe(2);
kam=pe(3);
kbcat=pe(4);
kbp=pe(5);
kbm=pe(6);


kone = pb(1);  
koffe = pb(2);  
kong = pb(3);  
koffg = pb(4);  
 
kong2 = pb(5);  
koffg2 = pb(6);  
kone2 = pb(7);  
koffe2 = pb(8);



TESTSCALE = pb(9); 
scl = pb(10);
% Sgpi = Pb(11);

% [Pgg, Pgp, Pgpgp]= ProbGamP(fracgp);
% 
% % Scale inhibiting fibrin end to end binding multiplies kpi, kpg and all
% % variants.
% kpscalegpi = Pgg + Sgpi*(Pgp +Pgpgp); 
% TESTSCALE
% end

% kone = 0;  
% koffe =0;  
% kong = 0;  
% koffg = 0;  
%  
% kong2 = 0;  
% koffg2 = 0;  
% kone2 = 0;  
% koffe2 = 0; 

fabalt=fab; 
% fab=fab; 
% 
% % falt=f*2;
% % f=f; 
% 
% whole=2; 
fbalt=fb;
% fb=1/whole*fb; 

% ODEs from reaction equations 
% whole = 1; 

%fibrinogen
% fab
 dy(1)  =- kap*T*fab + kam*Efab;% -  ka * fab;
 
 %fibrin 1
 % fb
 dy(2)  = - kbp*T*fb+ kbm*Efb + kacat*Efab...
     + ( -  2 * kpi * fb^2  -  kpi * f *fb  ...
     -  kpi21 * fb * f2  -  kpi31 * fb * f3  -  kpi41 * fb * f4  -  kpi51 * fb * f5  -  kpi61 * fb * f6...
     -  kpi71 * fb * f7  -  kpi81 * fb * f8  -  kpi91 * fb * f9  -  kpi101 * fb * f10 ...
     -  kpg * fb * fn)...
     ;%+ka*fab;
 
 
%   dy(2)  = - kbp*T*fbalt+ kbm*Efb + kacat*Efab...
%      + whole*( -  2 * kpi * fb^2  -  kpi * f *fb  ...
%      -  kpi21 * fb * f2  -  kpi31 * fb * f3  -  kpi41 * fb * f4  -  kpi51 * fb * f5  -  kpi61 * fb * f6  -  kpi71 * fb * f7  -  kpi81 * fb * f8  -  kpi91 * fb * f9  -  kpi101 * fb * f10 ...
%      -  kpg * fb * fn);
 % Efab
  dy(3)  = + kap*T*fab - kam*Efab - kacat*Efab;% -  ka * fa;

  
  % Efb
 dy(4)  = + kbp*T*fb - kbm*Efb - alpha*kbcat*Efb;% -  ka * fa;

% f
 dy(5)  =  alpha*kbcat*Efb ...
     -  2 * kpi * f^2  -  kpi * f *fb  ...
     -  kpi21 * f * f2  -  kpi31 * f * f3  -  kpi41 * f * f4  -  kpi51 * f * f5  -  kpi61 * f * f6  -  kpi71 * f * f7  -  kpi81 * f * f8  -  kpi91 * f * f9  -  kpi101 * f * f10 ...
     -  kpg * f * fn;
% dy(5)  =  1/whole*kbcat*Efb ...
%      -  2 * kpi * f^2  -  kpi * f *fb  ...
%      -  kpi21 * f * f2  -  kpi31 * f * f3  -  kpi41 * f * f4  -  kpi51 * f * f5  -  kpi61 * f * f6  -  kpi71 * f * f7  -  kpi81 * f * f8  -  kpi91 * f * f9  -  kpi101 * f * f10 ...
%      -  kpg * f * fn; 
% f2
 dy(6)  =  +  kpi * f^2  + kpi * f * fb + kpi * fb^2 ...
     -  kpi21 * (f+fb) * f2...
     -  2 * kpi22 * f2^2  -  kpi23 * f2 * f3  -  kpi24 * f2 * f4  -  kpi25 * f2 * f5  -  kpi26 * f2 * f6  -  kpi27 * f2 * f7  -  kpi28 * f2 * f8  -  kpi29 * f2 * f9 - kpi210*f2*f10...
     -  kpg2 * f2 * fn;

% f3
 dy(7)  =  +  kpi21 * (f + fb) * f2 ...
     -  kpi31 * (f+fb) * f3 ...
     -  kpi23 * f2 * f3    - 2 * kpi33 * f3^2  -  kpi34 * f3 * f4  -  kpi35 * f3 * f5  -  kpi36 * f3 * f6  -  kpi37 * f3 * f7  -  kpi38 * f3 * f8...
     - kpi39*f3*f9 - kpi310*f3*f10...
     -  kpg3 * f3 * fn;

% f4
 dy(8)  = +  kpi31 * (f+fb) * f3 +  kpi22 * f2^2 ...
      -  kpi41 * (f+fb) * f4  ...
     -  kpi24 * f2 * f4   -  kpi34 * f3 * f4  -  2 * kpi44 * f4^2  -  kpi45 * f4 * f5  -  kpi46 * f4 * f6  -  kpi47 * f4 * f7...
     -  kpi48*f4*f8 - kpi49*f4*f9 - kpi410*f4*f10...
     -  kpg4 * f4 * fn;

% f5
 dy(9)  =  +  kpi41 * (f +fb) * f4 +  kpi23 * f2 * f3 ...
     -  kpi51 * (f+fb) * f5 ...
     -  kpi25 * f2 * f5  -  kpi35 * f3 * f5    -  kpi45 * f4 * f5    - 2 * kpi55 * f5^2  -  kpi56 * f5 * f6...
     -  kpi57*f5*f7 -  kpi58*f5*f8 - kpi59*f5*f9 - kpi510*f5*f10...
     -  kpg5 * f5 * fn;

% f6
 dy(10)  = +  kpi51 * (f+fb) * f5   +  kpi24 * f2 * f4 +  kpi33 * f3^2 ...
     -  kpi61 * (f+fb) * f6 ...
     -  kpi26 * f2 * f6   -  kpi36 * f3 * f6  -  kpi46 * f4 * f6    -  kpi56 * f5 * f6 ...
     -  2*kpi66*f6*f6 -  kpi67*f6*f7 -  kpi68*f6*f8 - kpi69*f6*f9 - kpi610*f6*f10...
     -  kpg6 * f6 * fn;

% f7
 dy(11)  =  +  kpi61 * (f+fb) * f6 +  kpi25 * f2 * f5 +  kpi34 * f3 * f4 ...
      -  kpi71 * (f+fb) * f7...
     -  kpi27 * f2 * f7   -  kpi37 * f3 * f7  -  kpi47 * f4 * f7 ...
     -  kpi67*f6*f7 -  kpi57*f5*f7...
     -  2*kpi77*f7*f7  -  kpi78*f7*f8 - kpi79*f7*f9 - kpi710*f7*f10 ...
     -  kpg7 * f7 * fn;
 
% f8
 dy(12)  =  +  kpi71 * (f+fb) * f7  +  kpi26 * f2 * f6 +  kpi35 * f3 * f5 +  kpi44 * f4^2 ...
   -  kpi81 * (f+fb) * f8 ...
   -  kpi28 * f2 * f8    -  kpi38 * f3 * f8 ...
   -  kpi78*f7*f8 -  kpi68*f6*f8 -  kpi58*f5*f8 - kpi48*f4*f8...
   -  2*kpi88*f8*f8 - kpi89*f8*f9 - kpi810*f8*f10 ...
    -  kpg8 * f8 * fn;

% f9
 dy(13)  =  +  kpi81 * (f+fb) * f8 +  kpi27 * f2 * f7    +  kpi36 * f3 * f6  +  kpi45 * f4 * f5  ...
     -  kpi91 * (f+fb) * f9 ...
     -  kpi29 * f2 * f9 ...
     -  kpi89*f8*f9 -  kpi79*f7*f9 -  kpi69*f6*f9 -  kpi59*f5*f9 - kpi49*f4*f9 - kpi39*f3*f9...
     -  2*kpi99*f9*f9 - kpi910*f9*f10...
     -  kpg9 * f9 * fn;

% f10
 dy(14)  =  +  kpi91 * (f+fb) * f9 +  kpi28 * f2 * f8  +  kpi37 * f3 * f7  +  kpi46 * f4 * f6  +  kpi55 * f5^2  ...
     -  kpi101 * (f+fb) * f10...
     -  kpi910*f9*f10 -  kpi810*f8*f10 -  kpi710*f7*f10 -  kpi610*f6*f10 -  kpi510*f5*f10 - kpi410*f4*f10 - kpi310*f3*f10 - kpi210*f2*f10...
     -  2*kpi1010*f10*f10...
     -  kpg10 * f10 * fn ;

% %all additional olifomer interactions -- no repeated terms
%      +  kpi210*f2*f10 + kpi39*f3*f9 + kpi310*f3*f10 +  kpi48*f4*f8 + kpi49*f4*f9 + kpi410*f4*f10 ...
%      +  kpi57*f5*f7 +  kpi58*f5*f8 + kpi59*f5*f9 + kpi510*f5*f10 ...
%      +  2*kpi66*f6*f6 +  kpi67*f6*f7 +  kpi68*f6*f8 + kpi69*f6*f9 + kpi610*f6*f10...
%      +  2*kpi77*f7*f7  +  kpi78*f7*f8 + kpi79*f7*f9 + kpi710*f7*f10...
%      +  2*kpi88*f8*f8 + kpi89*f8*f9 + kpi710*f7*f10 +  2*kpi99*f9*f9 + kpi910*f9*f10 + 2*kpi1010*f10*f10 ...
%  %repeats    
%      + kpi67*f6*f7 +  kpi57*f5*f7...
%      +  kpi68*f6*f8 +  kpi58*f5*f8 + kpi48*f4*f8...
%      +  kpi69*f6*f9 +  kpi59*f5*f9 + kpi49*f4*f9 + kpi39*f3*f9...
%      +  kpi610*f6*f10 +  kpi510*f5*f10 + kpi410*f4*f10 + kpi310*f3*f10 + kpi210*f2*f10;
     
 
%Additional protofibril growth terms...
%      +  kpg2 * f2 * fn+  kpg3 * f3 * fn +  kpg4 * f4 * fn +  kpg5 * f5 * fn +  kpg6 * f6 * fn...
%      +  kpg7 * f7 * fn +  kpg8 * f8 * fn +  kpg9 * f9 * fn +  kpg10 * f10 * fn;
     
if kbcat == 0
      fracpol = 1; 
      FPBOFF = 1; 
      TESTSCALE=0; 
      TESTENHANCE =1;
else
%       TESTSCALE=0.5;
% TESTSCALE=1;%0.5;
%       fracpol= max(0,(cfn + TESTSCALE*cfb)/(cfb+cfn));  
fracpol= max(0,(cfn + TESTSCALE*cfb)/(cfb+cfn));  
      FPBOFF = 0 ;
      TESTENHANCE =1; 
end

% cfn
% cfb
fnp = fracpol * fn;  
 
% fn
 dy(15)  =  +  kpi29 * f2 * f9  +  kpi38 * f3 * f8  +  kpi47 * f4 * f7  +  kpi56 * f5 * f6...
     +  kpi101 * (f+fb) * f10 ...
     -  2 * kfi * fnp^2  -  kfg * fnp * fr ...
     +  kpi210*f2*f10 + kpi39*f3*f9 + kpi310*f3*f10 +  kpi48*f4*f8 + kpi49*f4*f9 + kpi410*f4*f10 ...
     +  kpi57*f5*f7 +  kpi58*f5*f8 + kpi59*f5*f9 + kpi510*f5*f10 ...
     +  kpi66*f6*f6 +  kpi67*f6*f7 +  kpi68*f6*f8 + kpi69*f6*f9 + kpi610*f6*f10...
     +  kpi77*f7*f7  +  kpi78*f7*f8 + kpi79*f7*f9 + kpi710*f7*f10...
     +  kpi88*f8*f8 + kpi89*f8*f9 + kpi810*f8*f10 +  kpi99*f9*f9 + kpi910*f9*f10 + kpi1010*f10*f10;

 
%  kfa=  0*kfi; 
% fr
 dy(16)  =  +  kfi * fnp^2;%  - kfa*fr*fr;
 
%ftotn
 dy(17) = 2*kfi*fnp^2 + kfg*fr*fnp;% +kfa*fr*fr;
 
 
 
%  gamma=1; 
  fBEBE1B=kbcat*(alpha*BE1 + beta*BE + gamma*B);
%  fBEBE1B=0;%kbcat*fracpol*(alpha*BE1 + beta*BE + gamma*B);
 
 %fraction of protofibril/oligomer fibrin that is in protofibrils
 %specifically
%  fracprot = max(0,(cfb+cfn)/(cfb+cfn+ cf+cb));
%  pause 
 

frac_cb=max(0,(cb)/(cfb+cfn+ cf+cb +f));
frac_cfb=max(0,(cfb)/(cfb+cfn+ cf+cb +f));


 %fraction of oligomoer fibrin with fpb
 fracb = max(0,cb/(cf+cb));
 
 
 
% 
 %cb
 dy(31) =+ 2*kpi*fb^2 + kpi*fb*f ...
     +  kpi21 * fb * f2  +  kpi31 * fb * f3  +  kpi41 * fb * f4  +  kpi51 * fb * f5  +  kpi61 * fb * f6  +  kpi71 * fb * f7  +  kpi81 * fb * f8  +  kpi91 * fb * f9 ...
     - 10*fracb*kpi101*f10*(f+fb) ...
     - 11*fracb*(kpi29*f2*f9+kpi38*f3*f8+kpi47*f4*f7+kpi56*f5*f6) ...
     - fBEBE1B*(frac_cb) ...
     - fracb * ( + 12*kpi210*f2*f10...
     +  12*kpi39*f3*f9 +  13*kpi310*f3*f10 ...
     +  12*kpi48*f4*f8 +  13*kpi49*f4*f9 +  14*kpi410*f4*f10 ...
     +  12*kpi57*f5*f7 +  13*kpi58*f5*f8 +  14*kpi59*f5*f9 + 15*kpi510*f5*f10 ...
     +  12*kpi66*f6*f6 +  13*kpi67*f6*f7 +  14*kpi68*f6*f8 + 15*kpi69*f6*f9 + 16*kpi610*f6*f10...
     +  14*kpi77*f7*f7 +  15*kpi78*f7*f8 + 16*kpi79*f7*f9 + 17*kpi710*f7*f10...
     +  16*kpi88*f8*f8 + 17*kpi89*f8*f9 + 18*kpi810*f8*f10...
     +  18*kpi99*f9*f9 + 19*kpi910*f9*f10...
     + 20*kpi1010*f10*f10) ...
     - fracb * (+  2*kpg2 * f2 * fn+  3*kpg3 * f3 * fn +  4*kpg4 * f4 * fn +  5*kpg5 * f5 * fn +  6*kpg6 * f6 * fn...
     +  7* kpg7 * f7 * fn +  8*kpg8 * f8 * fn +  9*kpg9 * f9 * fn +  10*kpg10 * f10 * fn);
%      - fBEBE1B*(1-fracprot);

%cf
 dy(32) =+ 2*kpi*f^2 + kpi*fb*f ...
     +  kpi21 * f * f2  +  kpi31 * f * f3  +  kpi41 * f * f4  +  kpi51 * f * f5  +  kpi61 * f * f6  +  kpi71 * f * f7  +  kpi81 * f * f8  +  kpi91 * f * f9 ...
     - 10*(1-fracb)*kpi101*f10*(f+fb) ...
     - 11*(1-fracb)*(kpi29*f2*f9+kpi38*f3*f8+kpi47*f4*f7+kpi56*f5*f6) ...
     + fBEBE1B*(frac_cb)...
     - (1-fracb) * ( + 12*kpi210*f2*f10...
     +  12*kpi39*f3*f9 +  13*kpi310*f3*f10...
     +  12*kpi48*f4*f8 +  13*kpi49*f4*f9 +  14*kpi410*f4*f10 ...
     +  12*kpi57*f5*f7 +  13*kpi58*f5*f8 +  14*kpi59*f5*f9 + 15*kpi510*f5*f10 ...
     +  12*kpi66*f6*f6 +  13*kpi67*f6*f7 +  14*kpi68*f6*f8 + 15*kpi69*f6*f9 + 16*kpi610*f6*f10...
     +  14*kpi77*f7*f7 +  15*kpi78*f7*f8 + 16*kpi79*f7*f9 + 17*kpi710*f7*f10...
     +  16*kpi88*f8*f8 + 17*kpi89*f8*f9 + 18*kpi810*f8*f10 ...
     +  18*kpi99*f9*f9 + 19*kpi910*f9*f10 ...
     +  20*kpi1010*f10*f10) ...
     - (1-fracb) * (+  2*kpg2 * f2 * fn+  3*kpg3 * f3 * fn +  4*kpg4 * f4 * fn +  5*kpg5 * f5 * fn +  6*kpg6 * f6 * fn...
     +  7* kpg7 * f7 * fn +  8*kpg8 * f8 * fn +  9*kpg9 * f9 * fn +  10*kpg10 * f10 * fn);
 
 
 
  %cfb
 dy(18) = +kpi101*fb*f10 ...
     + 10*(fracb)*kpi101*f10*(f+fb) ...
     +11*fracb*(kpi29*f2*f9+kpi38*f3*f8+kpi47*f4*f7+kpi56*f5*f6) ...
     + kpg*fn*fb ...
     - fBEBE1B*(frac_cfb) ...
     + TESTSCALE*(- 2*kfi*fnp*cfb -kfg*fr*cfb)...
     + fracb * ( + 12*kpi210*f2*f10...
     + 12*kpi39*f3*f9 + 13*kpi310*f3*f10...
     +  12*kpi48*f4*f8 + 13*kpi49*f4*f9 + 14*kpi410*f4*f10 ...
     +  12*kpi57*f5*f7 +  13*kpi58*f5*f8 + 14*kpi59*f5*f9 + 15*kpi510*f5*f10 ...
     +  12*kpi66*f6*f6 +  13*kpi67*f6*f7 +  14*kpi68*f6*f8 + 15*kpi69*f6*f9 + 16*kpi610*f6*f10...
     +  14*kpi77*f7*f7  +  15*kpi78*f7*f8 + 16*kpi79*f7*f9 + 17*kpi710*f7*f10...
     +  16*kpi88*f8*f8 + 17*kpi89*f8*f9 + 18*kpi810*f8*f10 ...
     +  18*kpi99*f9*f9 + 19*kpi910*f9*f10...
     +  20*kpi1010*f10*f10) ...
     + fracb * (+  2*kpg2 * f2 * fn +  3*kpg3 * f3 * fn +  4*kpg4 * f4 * fn +  5*kpg5 * f5 * fn +  6*kpg6 * f6 * fn...
     +  7* kpg7 * f7 * fn +  8*kpg8 * f8 * fn +  9*kpg9 * f9 * fn +  10*kpg10 * f10 * fn)...
     + (- 2*kfi*fnp*cfb -kfg*fr*cfb )*FPBOFF; %this term should only be used if kbcat = 0... no fpB case 

%cfn
 dy(19) =+kpi101*f*f10 ...
     + 10*(1-fracb)*kpi101*f10*(f+fb) ...
     +11*(1-fracb)*(kpi29*f2*f9+kpi38*f3*f8+kpi47*f4*f7+kpi56*f5*f6) ...
     + kpg*fn*f ...
     + fBEBE1B*(frac_cfb)...
     +TESTENHANCE*(- 2*kfi*fnp*cfn -kfg*fr*cfn) ...
     + (1-fracb) * ( + 12*kpi210*f2*f10 ...
     + 12*kpi39*f3*f9 + 13*kpi310*f3*f10 ...
     +  12*kpi48*f4*f8 + 13*kpi49*f4*f9 + 14*kpi410*f4*f10 ...
     +  12*kpi57*f5*f7 +  13*kpi58*f5*f8 + 14*kpi59*f5*f9 + 15*kpi510*f5*f10 ...
     +  12*kpi66*f6*f6 +  13*kpi67*f6*f7 +  14*kpi68*f6*f8 + 15*kpi69*f6*f9 + 16*kpi610*f6*f10...
     +  14*kpi77*f7*f7  +  15*kpi78*f7*f8 + 16*kpi79*f7*f9 + 17*kpi710*f7*f10...
     +  16*kpi88*f8*f8 + 17*kpi89*f8*f9 + 18*kpi810*f8*f10...
     +  18*kpi99*f9*f9 + 19*kpi910*f9*f10 ...
     + 20*kpi1010*f10*f10) ...
     + (1-fracb) * (+  2*kpg2 * f2 * fn+  3*kpg3 * f3 * fn +  4*kpg4 * f4 * fn +  5*kpg5 * f5 * fn +  6*kpg6 * f6 * fn...
     +  7* kpg7 * f7 * fn +  8*kpg8 * f8 * fn +  9*kpg9 * f9 * fn +  10*kpg10 * f10 * fn);
 
% %  
 
 

%cfr
 dy(20) = + TESTENHANCE*(2*kfi*fnp*cfn +kfg*fr*cfn)...
     + TESTSCALE*( 2*kfi*fnp*cfb +kfg*fr*cfb)...
     +( 2*kfi*fnp*cfb +kfg*fr*cfb )*FPBOFF;





% dfp=- 2 * kpi * f^2  -  kpi21 * f * f2  -  kpi31 * f * f3  -  kpi41 * f * f4  -  kpi51 * f * f5  -  kpi61 * f * f6  -  kpi71 * f * f7  -  kpi81 * f * f8  -  kpi91 * f * f9  -  kpi101 * f * f10  -  kpg * f * fn;
% 
% % 
% % dE= +(2-fracgp)*kbcat*Efb;%2*totalnotpolym+1.7*totalpolym ;
% % dEp=+fracgp*kbcat*Efb;%0.3*totalpolym; 
% % dG=+fracgp*kbcat*Efb;%0.3*totalpolym;
% 
% % dE= +2*kbcat*Efb+fracgp*dfp;%2*totalnotpolym+1.7*totalpolym ;
% % dEp=-fracgp*dfp;%+fracgp*kbcat*Efb;%0.3*totalpolym; 
% % dG=+fracgp*kbcat*Efb;%-fracgp*dfp%+fracgp*kbcat*Efb;%0.3*totalpolym;
% % dGp=-fracgp*kbcat*Efb;%+fracgp*dfp
% % %
% dE= +2*kacat*Efab+fracgp*dfp;%2*totalnotpolym+1.7*totalpolym ;
% dEp=-fracgp*dfp;%+fracgp*kbcat*Efb;%0.3*totalpolym; 
% dG=+fracgp*kacat*Efab;%-fracgp*dfp%+fracgp*kbcat*Efb;%0.3*totalpolym;
% dGp=-fracgp*kacat*Efab;%+fracgp*dfp



fracBE1 = max(0,BE1/(E1+BE1));
fracBE = max(0,BE/(E+BE));
fracBG1 = max(0,BG1/(G1+BG1));
fracBG =max(0,BG/(G+BG));

% 
% ( -  2 * kpi * fb^2  -  kpi * f *fb  ...
%      -  kpi21 * fb * f2  -  kpi31 * fb * f3  -  kpi41 * fb * f4  -  kpi51 * fb * f5  -  kpi61 * fb * f6  -  kpi71 * fb * f7  -  kpi81 * fb * f8  -  kpi91 * fb * f9  -  kpi101 * fb * f10 ...
%      -  kpg * fb * fn)

dfbp=  abs(-  2 * kpi * fb^2  -  kpi * f *fb  ...
     -  kpi21 * fb * f2  -  kpi31 * fb * f3  -  kpi41 * fb * f4  -  kpi51 * fb * f5  -  kpi61 * fb * f6  -  kpi71 * fb * f7  -  kpi81 * fb * f8  -  kpi91 * fb * f9  -  kpi101 * fb * f10 ...
     -  kpg * fb * fn);

dfp = abs(-  2 * kpi * f^2  -  kpi * f *fb  ...
     -  kpi21 * f * f2  -  kpi31 * f * f3  -  kpi41 * f * f4  -  kpi51 * f * f5  -  kpi61 * f * f6  -  kpi71 * f * f7  -  kpi81 * f * f8  -  kpi91 * f * f9  -  kpi101 * f * f10 ...
     -  kpg * f * fn);

 
dE1 = + 2*kbcat*Efb + dfbp*(2-fracgp) - dfp*(fracgp)*(1-fracBE1);
dBE1 =  - dfp*(fracgp)*fracBE1;

dE =  + dfbp*(fracgp) + dfp*(fracgp)*(1-fracBE1);
dBE = dfp*(fracgp)*fracBE1; 

dG1 = (- dfp*(fracgp) -dfbp*(fracgp))*(1-fracBG1);
dBG1 = (- dfp*(fracgp) -dfbp*(fracgp))*(fracBG1);

dG = (+ dfp*(fracgp) +dfbp*(fracgp))*(1-fracBG1);
dBG = (+ dfp*(fracgp) +dfbp*(fracgp))*(fracBG1);
% dfbp
% dfp
% dBG1
% dG1
% fracBG1
% scl =.5; 

%  +(+ 2*kpi * fb^2 + kpi * f *fb  +  kpi21 * fb * f2  +  kpi31 * fb * f3 ...
%     +  kpi41 * fb * f4  +  kpi51 * fb * f5  +  kpi61 * fb * f6  ...
%     + kpi71 * fb * f7  +  kpi81 * fb * f8  +  kpi91 * fb * f9  +  kpi101 * fb * f10)*(1-fracgp) * (fracBE1)...
% T
 dy(21)  =  -  kone * T * E1  +  koffe * BE1  ...
     -  kone * T * E  +  koffe * BE ...
      -  kong * T * G1  +  koffg * BG1...
     -  kong * T * G  +  koffg * BG...
     - kap*T*fab+ kam*Efab +kacat*Efab...
     - kbp*T*fb+ kbm*Efb + alpha*kbcat*Efb...
     + kbcat*(alpha*BE1 + beta*BE)*(frac_cfb + frac_cb)+ kbcat*(scl * gamma*B)*(frac_cfb + frac_cb)*ThromScale;
%      + kbcat*(alpha*BE1 + beta*BE )*(frac_cfb + frac_cb)*ThromScale;
 
%  kbcat*alpha*BE1 
%  kbcat*beta*BE 
%  kbcat*gamma*B
    
% size(fracgp)


% E1
 dy(22)  =  -  kone * T * E1  +  koffe * BE1 +  dE1 + (frac_cfb + frac_cb)*kbcat*alpha*BE1;

% BE1
 dy(23)  =  +  kone * T * E1  -  koffe * BE1 + dBE1 - kbcat*alpha*BE1*(frac_cfb + frac_cb);

% E
 dy(24)  =  -  kone * T * E  +  koffe * BE -  kone2 * BG  +  koffe2 * B + dE +  kbcat*beta*BE*(frac_cfb + frac_cb) + kbcat*gamma*B*(frac_cfb + frac_cb)*ThromScale;

% BE
 dy(25)  =  +  kone * T * E  -  koffe * BE  -  kong2 * BE  +  koffg2 * B + dBE - kbcat*beta*BE*(frac_cfb + frac_cb) ;

  % G1
 dy(26)  =  -  kong * T * G1  +  koffg * BG1 + dG1;
% BG1
 dy(27)  =  +  kong * T * G1  -  koffg * BG1 + dBG1;

 % G
 dy(28)  =  -  kong * T * G  +  koffg * BG -  kong2 * BE  +  koffg2 * B +  dG + scl*kbcat*gamma*B*(frac_cfb + frac_cb)*ThromScale;
 

% BG
 dy(29)  =  +  kong * T * G  -  koffg * BG  -  kone2 * BG  +  koffe2 * B+ dBG + (1-scl)*kbcat*gamma*B*(frac_cfb + frac_cb)*ThromScale;

% B
 dy(30)  =  +  kong2 * BE  -  koffg2 * B  +  kone2 * BG  -  koffe2 * B -  kbcat*gamma*B*(frac_cfb + frac_cb)*ThromScale;
 
 

 
%  cfb
%  cfn
%  cb
%  cf
% %  

%  pause
%  fracprot
%  fracb
%  
%   fb
%   Efab
%   Efb
%  f
%  f2
%  f3
%  f4
%  f5
%  f6
%  f7
%  f8
%  f9
%  f10
%   fn
%  dfbp
%  dfp
%  pause
% % diff= abs(sum(dy(6:14))-sum(dy(31:32)))
% dy(6:14)
% F=sum(dy(6:14).*(2:10)')
% F2=sum(dy(31:32))
%  if t>1.2
% %       fracpol
% %  fracprot
% %  fracb
% %      pause  
% 
% % pause
%  end
 
 

end


function [Pgg, Pgp, Pgpgp]= ProbGamP(fracgp)
Pgg = (1-fracgp)^2;
Pgp = 2 * (1-fracgp)*fracgp; 
Pgpgp =  fracgp^2; 



end