% plot fig 2 violins, need to import G1 times first
[C]=cbrewer('div','RdBu' , 3);

binsG0=[1:1:100];  %since there are 3 frames in an hour, this bin makes sense / surpressing zero helps in later fitting
binsSM=[1:1:50];
vio_cntr=0.2;
vio_room=0.05;
g0_centers=[0.5:1:99.5];

% pos1 = [0.1 0.55 0.3 0.4];
% pos2 = [0.4 0.55 0.3 0.4];
% pos3 = [0.7 0.55 0.3 0.4];

g1_emg_fits = zeros(6,3);
ft_lognc = fittype('(1./(x*sig*sqrt(2*pi)))*exp(-(log(x)-mu)^2./(2*sig^2))'...
    ,'options',fo_lognc);
%%
fig2 = figure; hold on;
fig2.Color = 'w';
fig2.Position = [100 100 800 400];

clear gca
sb1=subplot(1,3,1); hold on
box on
set(gca,'linewidth',1.5,'fontsize',15)
sb1.XLim = [-0.15 0.15]; sb1.XTick = [-0.15 -0.1 -.05 0 0.05 0.1 0.15];
Xlabels = {[],0.1,[],0,[],0.1,[]};
sb1.XTickLabel = Xlabels;
Ylabels = {0,[],20,[],40,[],60,[],80,[],100};
sb1.YTick = [0:10:100];
sb1.YTickLabel = Ylabels;
sb1.YGrid = 'on';
ax = gca;
gca.Xaxis.LineWidth = 2; gca.Yaxis.LineWidth = 2;
[a, b]=histcounts(Ds8_02kPa_DMS0_G0(2,:),binsG0,'Normalization','probability');
a(100)=0; a(1)=0;%have to pad to make even to use stairs fun?
[xb, yb] = stairs(b,a);
%log-normal fit
[logn_curve_g0,gofln] = fit(b',a',ft_lognc);
logn_sol=logn_curve_g0(binsG0);
%exp-gauss fit
[exgfit, egci] = exgaussErr(Ds8_02kPa_DMS0_G0(2,:),binsG0);
egtau(1,1)=exgfit(1); egtau(1,2)=egci(1,1); egtau(1,3)=egci(1,2);
egmu(1,1)=exgfit(2); egmu(1,2)=egci(2,1); egmu(1,3)=egci(2,2);
egsig(1,1)=exgfit(3); egsig(1,2)=egci(3,1); egsig(1,3)=egci(3,2);
egp_sol=exgausspdf(egmu(1,1),egsig(1,1),egtau(1,1),b);
%gamma fit
gam=gamfit(Ds8_02kPa_DMS0_G0(2,:));
gam_sol=gampdf(0:1:100,gam(1),gam(2));
%plot solutions
fill(-yb,xb,C(3,:)...
    ,'FaceAlpha',0.8,'EdgeColor','k');
%plot(-logn_sol,binsG0,'k');
plot(-logn_sol,binsG0,'Color',[0.15 0.23 0.37],'LineWidth',2);
%plot(-gam_sol,0:1:100,'r');

[a, b]=histcounts(Ds8_02kPa_ERL_G0(2,:),binsG0,'Normalization','probability');
a(100)=0; a(1)=0;%have to pad to make even to use stairs fun?
[xb, yb] = stairs(b,a);
%log-normal fit
[logn_curve_g0,gofln] = fit(b',a',ft_lognc);
logn_sol=logn_curve_g0(binsG0);
%exp-gauss fit
[exgfit, egci] = exgaussErr(Ds8_02kPa_ERL_G0(2,:),binsG0);
egtau(2,1)=exgfit(1); egtau(2,2)=egci(1,1); egtau(2,3)=egci(1,2);
egmu(2,1)=exgfit(2); egmu(2,2)=egci(2,1); egmu(2,3)=egci(2,2);
egsig(2,1)=exgfit(3); egsig(2,2)=egci(3,1); egsig(2,3)=egci(3,2);
egp_sol=exgausspdf(egmu(2,1),egsig(2,1),egtau(2,1),b);
%gamma fit
gam=gamfit(Ds8_02kPa_ERL_G0(2,:)); 
gam_sol=gampdf(0:1:100,gam(1),gam(2));
fill(yb,xb,[0.6350, 0.0780, 0.1840]...
    ,'FaceAlpha',0.6,'EdgeColor','k');
%plot(logn_sol,binsG0,'k');
plot(logn_sol,binsG0,'Color',[0.6350, 0.0780, 0.1840],'LineWidth',2);
%plot(gam_sol,0:1:100,'r');

%%
clear gca
sb2=subplot(1,3,2); hold on
box on
set(gca,'linewidth',1.5,'fontsize',15,'YTickLabel',[])
%sb2.YAxis.Visible = 'off';
sb2.XLim = [-0.15 0.15]; sb2.XTick = [-0.15 -0.1 -.05 0 0.05 0.1 0.15];
sb2.XTickLabel = Xlabels;
sb2.YGrid = 'on';
sb2.YTick = [0:10:100];
[a, b]=histcounts(Ds8_4kPa_DMS0_G0(2,:),binsG0,'Normalization','probability');
a(100)=0; a(1)=0;%have to pad to make even to use stairs fun?
[xb, yb] = stairs(b,a);
%log-normal fit
[logn_curve_g0,gofln] = fit(b',a',ft_lognc);
logn_sol=logn_curve_g0(binsG0);
%exp-gauss fit
[exgfit, egci] = exgaussErr(Ds8_4kPa_DMS0_G0(2,:),binsG0);
egtau(3,1)=exgfit(1); egtau(3,2)=egci(1,1); egtau(3,3)=egci(1,2);
egmu(3,1)=exgfit(2); egmu(3,2)=egci(2,1); egmu(3,3)=egci(2,2);
egsig(3,1)=exgfit(3); egsig(3,2)=egci(3,1); egsig(3,3)=egci(3,2);
egp_sol=exgausspdf(egmu(3,1),egsig(3,1),egtau(3,1),b);
%gamma fit
% gam=gamfit(Ds8_4kPa_DMS0_G0(2,:));
% gam_sol=gampdf(0:1:100,gam(1),gam(2));
%plot solutions
fill(-yb,xb,C(3,:)...
    ,'FaceAlpha',0.8,'EdgeColor','k');
%plot(-logn_sol,binsG0,'k');
plot(-logn_sol,binsG0,'Color',[0.15 0.23 0.37],'LineWidth',2);
%plot(-gam_sol,0:1:100,'r');

[a, b]=histcounts(Ds8_4kPa_ERL_G0(2,:),binsG0,'Normalization','probability');
a(100)=0; a(1)=0;%have to pad to make even to use stairs fun?
[xb, yb] = stairs(b,a);
%log-normal fit
[logn_curve_g0,gofln] = fit(b',a',ft_lognc);
logn_sol=logn_curve_g0(binsG0);
%exp-gauss fit
[exgfit, egci] = exgaussErr(Ds8_4kPa_ERL_G0(2,:),binsG0);
egtau(4,1)=exgfit(1); egtau(4,2)=egci(1,1); egtau(4,3)=egci(1,2);
egmu(4,1)=exgfit(2); egmu(4,2)=egci(2,1); egmu(4,3)=egci(2,2);
egsig(4,1)=exgfit(3); egsig(4,2)=egci(3,1); egsig(4,3)=egci(3,2);
egp_sol=exgausspdf(egmu(4,1),egsig(4,1),egtau(4,1),b);
%gamma fit
% gam=gamfit(Ds8_4kPa_ERL_G0(2,:));
% gam_sol=gampdf(0:1:100,gam(1),gam(2));
fill(yb,xb,[0.6350, 0.0780, 0.1840]...
    ,'FaceAlpha',0.6,'EdgeColor','k');
%plot(logn_sol,binsG0,'k');
plot(logn_sol,binsG0,'Color',[0.6350, 0.0780, 0.1840],'LineWidth',2);
%plot(gam_sol,0:1:100,'r');

%%
clear gca
sb3=subplot(1,3,3); hold on
box on
set(gca,'linewidth',1.5,'fontsize',15,'YTickLabel',[])
%sb3.YAxis.Visible = 'off';
sb3.XLim = [-0.15 0.15]; sb3.XTick = [-0.15 -0.1 -.05 0 0.05 0.1 0.15];
sb3.XTickLabel = Xlabels;
sb3.YGrid = 'on';
sb3.YTick = [0:10:100];
[a, b]=histcounts(Ds8_55kPa_DMS0_G0(2,:),binsG0,'Normalization','probability');
a(100)=0; a(1)=0;%have to pad to make even to use stairs fun?
[xb, yb] = stairs(b,a);
%log-normal fit
[logn_curve_g0,gofln] = fit(b',a',ft_lognc);
logn_sol=logn_curve_g0(binsG0);
%exp-gauss fit
[exgfit, egci] = exgaussErr(Ds8_55kPa_DMS0_G0(2,:),binsG0);
egtau(5,1)=exgfit(1); egtau(5,2)=egci(1,1); egtau(5,3)=egci(1,2);
egmu(5,1)=exgfit(2); egmu(5,2)=egci(2,1); egmu(5,3)=egci(2,2);
egsig(5,1)=exgfit(3); egsig(5,2)=egci(3,1); egsig(5,3)=egci(3,2);
egp_sol=exgausspdf(egmu(5,1),egsig(5,1),egtau(5,1),b);
%gamma fit
gam=gamfit(Ds8_55kPa_DMS0_G0(2,:));
gam_sol=gampdf(0:1:100,gam(1),gam(2));
%plot solutions
fill(-yb,xb,C(3,:)...
    ,'FaceAlpha',0.8,'EdgeColor','k');
%plot(-logn_sol,binsG0,'k');
plot(-logn_sol,binsG0,'Color',[0.15 0.23 0.37],'LineWidth',2);
%plot(-gam_sol,0:1:100,'r');

[a, b]=histcounts(Ds8_55kPa_ERL_G0(2,:),binsG0,'Normalization','probability');
a(100)=0; a(1)=0;%have to pad to make even to use stairs fun?
[xb, yb] = stairs(b,a);
%log-normal fit
[logn_curve_g0,gofln] = fit(b',a',ft_lognc);
logn_sol=logn_curve_g0(binsG0);
%exp-gauss fit
[exgfit, egci] = exgaussErr(Ds8_55kPa_ERL_G0(2,:),binsG0);
egtau(6,1)=exgfit(1); egtau(6,2)=egci(1,1); egtau(6,3)=egci(1,2);
egmu(6,1)=exgfit(2); egmu(6,2)=egci(2,1); egmu(6,3)=egci(2,2);
egsig(6,1)=exgfit(3); egsig(6,2)=egci(3,1); egsig(6,3)=egci(3,2);
egp_sol=exgausspdf(egmu(6,1),egsig(6,1),egtau(6,1),b);
%gamma fit
gam=gamfit(Ds8_55kPa_ERL_G0(2,:));
gam_sol=gampdf(0:1:100,gam(1),gam(2));
fill(yb,xb,[0.6350, 0.0780, 0.1840]...
    ,'FaceAlpha',0.6,'EdgeColor','k');
%plot(logn_sol,binsG0,'k');
plot(logn_sol,binsG0,'Color',[0.6350, 0.0780, 0.1840],'LineWidth',2);
%plot(gam_sol,0:1:100,'r');
%%
