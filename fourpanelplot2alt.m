function fourpanelplot2(sol,fignum,FS,legnd)


pf=[0, 30,60,90,120,150,180,210,240];
diameter=[54, 75.5, 92, 105.5, 118, 128.7,139.1,148.7];
numb_innac=[0, 4,17,35,54,75,96,120,141];
frac_innac=numb_innac./pf;
frac_innac(1)=0;

time=sol.time;
t_final=time(end);
% 
% 
% 
% 
% 


 totalpolym= sol.cfn+sol.cfr+2*sol.f2+3*sol.f3+4*sol.f4+5*sol.f5 ...
       +6*sol.f6+7*sol.f7+8*sol.f8+9*sol.f9+10*sol.f10;
totalnotpolym= sol.f;%sol.fa+sol.f;
%     
%   totalthrombin=sol.BE+sol.BEp+sol.BG+sol.B+ sol.Efb+sol.Efab+sol.T
% totalfibrin= sol.cfn+sol.cfr+2*sol.f2+3*sol.f3+4*sol.f4+5*sol.f5 ...
%        +6*sol.f6+7*sol.f7+8*sol.f8+9*sol.f9+10*sol.f10+sol.f+sol.fab+sol.fb
% TotalE=sol.E+sol.BE
% TotalEp=sol.Ep+sol.BEp
% TotalG=sol.G+sol.BG
% TotalB=sol.B
% pause 



% E=2*totalnotpolym+1.7*totalpolym;% + sol.Efab+sol.Efb;
E1=1.7*totalpolym+2*totalnotpolym;
E=0.3*totalpolym; 
G=0.3*totalpolym;
 
fracE1=1.7*(sol.cfr)./E1;
fracE=0.3*(sol.cfr)./E;
fracG=0.3*(sol.cfr)./G;
% fracE=(sol.cfr)./E;
% fracEp=(sol.cfr)./Ep;
% fracG=(sol.cfr)./G;

fracEp(1)=0;
fracG(1)=0;
totalboundthrombin=sol.BE1+sol.BE+sol.BG+sol.B+sol.BG1+ sol.Efb+sol.Efab;
trappedthrombin=sol.BE1.*fracE1+sol.BE.*fracE+sol.BG.*fracG+sol.B.*fracG;
percenttrappedthrombin=(trappedthrombin)./totalboundthrombin;

sol.percentboundthrombin=percenttrappedthrombin;
save('examplesol.mat','sol')

% pttm=percenttrappedthrombin(percenttrappedthrombin<.3);
% pttp=percenttrappedthrombin(percenttrappedthrombin>.3);
% 
% tind=length(pttm)+1;
% ttrap=time(tind);




m=sol.m;
max(m);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(fignum(1))

% subplot(2,2,2)
% hold on
% % set(gca,'ColorOrderIndex',1)
% plot(time,100*percenttrappedthrombin)
% 
% %  plot(time(1:tind),ones(1,length(time(1:tind)))*percenttrappedthrombin(tind)*100,'r')
% %  plot(ttrap*ones(100,1),linspace(0,percenttrappedthrombin(tind)*100,100),'r')
%  xlim([0 t_final]);
% title('Percent of fiber bound thrombin');
% ylabel('% thrombin')
% xlabel('time (s)')
% set(gca,'FontSize',FS)
% % legend('0','0.05','0.1','0.15','0.2','0.3','1')
% % legend('100 nM thrombin','10 nM thrombin','1 nM thormbin','0.1 nM thormbin')
% 

 switch legnd 
     case 'BS'
%           legend('0 (\gamma'' BS''s)/(E domain BS''s)','0.3 (\gamma'' BS''s)/(E domain BS''s)','1(\gamma'' BS''s)/(E domain BS''s)','2 (\gamma'' BS''s)/(E domain BS''s)')
subplot(1,2,2)
     case 'Thrombin'
%           legend('100 nM thrombin','10 nM thrombin','1 nM thrombin','0.1 nM thrombin')
subplot(1,2,1)
 end

% subplot(2,2,1)
hold on
plot(time,sol.m);

%  plot(ttrap*ones(100,1),linspace(0,m(tind),100),'r')
%  plot(time(1:tind),ones(1,length(time(1:tind)))*m(tind),'r')
% ylim([0 (max(sol.m)+5)]);
 xlim([0 t_final]);
title('# of bundled protofibrils per fiber')
ylabel('# protofibrils')
xlabel('time (s)')
set(gca,'FontSize',FS)
% legend('0','0.05','0.1','0.15','0.2','0.3','1')
% legend('100 nM thrombin','10 nM thrombin','1 nM thormbin','0.1 nM thormbin')
% 
% vq=interp1(x,v,xq);


switch legnd 
     case 'BS'
          legend('2:1 \gamma'' per E-domain','1:1 \gamma'' per E-domain','0.3:1 \gamma'' per E-domain','0:1 \gamma'' per E-domain')

     case 'Thrombin'
          legend('100 nM thrombin','10 nM thrombin','1 nM thrombin','0.1 nM thrombin')
 end
fracvsm=interp1(pf,frac_innac,sol.m);

% subplot(2,2,3)
% hold on
% %  plot(pf,100*frac_innac)
%   plot(time,sol.fr)
% % ylim([0 (max(sol.m)+5)]);
% %  xlim([0 t_final]);
% % title('Percent innacesible fibrin binding sites')
% % ylabel('# fibers')
% % xlabel('Number fibers')
% set(gca,'FontSize',FS)
% legend('100 nM thrombin','10 nM thrombin','1 nM thrombin','0.1 nM thrombin')
% 
% title('Concentration of fibrin fibers')
% ylabel('fibers (\muM)')
% xlabel('time (s)')
%  xlim([0 t_final]);


% subplot(2,2,3)
% hold on
%  plot(pf,100*frac_innac)
% % ylim([0 (max(sol.m)+5)]);
% %  xlim([0 t_final]);
% title('Percent innacesible fibrin binding sites')
% ylabel('% (protofibrils)')
% xlabel('Number protofibrils')
% set(gca,'FontSize',FS)
% % legend('0','0.05','0.1','0.15','0.2','0.3','0.75','1')
% % legend('10 nM thrombin','1 nM thrombin','0.1 nM thormbin','0.01 nM thormbin')

figure(fignum(2))

 switch legnd 
     case 'BS'
%           legend('0 (\gamma'' BS''s)/(E domain BS''s)','0.3 (\gamma'' BS''s)/(E domain BS''s)','1(\gamma'' BS''s)/(E domain BS''s)','2 (\gamma'' BS''s)/(E domain BS''s)')
subplot(1,2,2)
     case 'Thrombin'
%           legend('100 nM thrombin','10 nM thrombin','1 nM thrombin','0.1 nM thrombin')
subplot(1,2,1)
 end

% subplot(2,2,4)
hold on
 plot(time, 100*percenttrappedthrombin.*fracvsm)
% ylim([0 (max(sol.m)+5)]);
 xlim([0 t_final]);
title('Percent trapped thrombin')
ylabel('% thrombin')
xlabel('time (s)')
set(gca,'FontSize',FS)
% legend('100 nM thrombin','10 nM thrombin','1 nM thormbin','0.1 nM thormbin')


switch legnd 
     case 'BS'
%           legend('0 (\gamma'' BS''s)/(E domain BS''s)','0.3 (\gamma'' BS''s)/(E domain BS''s)','1(\gamma'' BS''s)/(E domain BS''s)','2 (\gamma'' BS''s)/(E domain BS''s)')
%           legend('2 (\gamma'' BS''s)/(E domain BS''s)','1(\gamma'' BS''s)/(E domain BS''s)','0.3 (\gamma'' BS''s)/(E domain BS''s)','0 (\gamma'' BS''s)/(E domain BS''s)')
legend('2:1 \gamma'' per E-domain','1:1 \gamma'' per E-domain','0.3:1 \gamma'' per E-domain','0:1 \gamma'' per E-domain')
     case 'Thrombin'
          legend('100 nM thrombin','10 nM thrombin','1 nM thrombin','0.1 nM thrombin')
 end
