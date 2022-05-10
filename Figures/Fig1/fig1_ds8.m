%let's make some figures!
t = (0:419)/3;
t3 = (0:362)/3;

DIP_DS8_0p2kPa_DMSO = importdata('DIP_DS8_0p2kPa_DMSO.txt');
DIP_DS8_0p2kPa_ERL = importdata('DIP_DS8_0p2kPa_ERL.txt');
DIP_DS8_4kPa_DMSO = importdata('DIP_DS8_4kPa_DMSO.txt');
DIP_DS8_4kPa_ERL = importdata('DIP_DS8_4kPa_ERL.txt');
DIP_DS8_55kPa_DMSO = importdata('DIP_DS8_55kPa_DMSO.txt');
DIP_DS8_55kPa_ERL = importdata('DIP_DS8_55kPa_ERL.txt');

fig1a=figure; hold on
fig1a.Color='w';
fig1a.Position = [100 100 800 300];
%%
%ds8 0.2
[C]=cbrewer('div','RdBu' , 3);

sb1=subplot(1,3,1);
title('0.2 kPa','fontsize',14)
sb1.Position(2)=0.2;
sb1.Position(4)=0.7;
axis([0 120 -1 4])
set(gca,'linewidth',1.5,'fontsize',12,'TickLength',[0.03,0.03],'TickDir','out')
grid on
a=gca;
set(a,'box','off','color','none')
b=axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
set(gca,'linewidth',1.5,'fontsize',12)
axes(a)
linkaxes([a b])

%dmso
dip_tmp=zeros(length(DIP_DS8_0p2kPa_DMSO),2);
for ii = 1:length(DIP_DS8_0p2kPa_DMSO)
    dip_tmp(ii,1)=mean(DIP_DS8_0p2kPa_DMSO(ii,:));
    dip_tmp(ii,2)=std(DIP_DS8_0p2kPa_DMSO(ii,:));
end
fit_tmp = fit(t(216:end)',dip_tmp(216:end,1),'poly1');

s1d=shadedErrorBar(t,smooth(smooth(dip_tmp(:,1),7),15),smooth(smooth(dip_tmp(:,2),7),15));
% Set face and edge properties
set(s1d.edge,'Visible','off')
set(s1d.mainLine,'LineWidth',2,'color',C(3,:))
s1d.patch.FaceColor = C(3,:);
hold on
tmpx = 72:1:120;
tmpy = fit_tmp(tmpx);
plot(tmpx,tmpy,'k');

%erl
dip_tmp=zeros(length(DIP_DS8_0p2kPa_ERL),2);
for ii = 1:length(DIP_DS8_0p2kPa_ERL)
    dip_tmp(ii,1)=mean(DIP_DS8_0p2kPa_ERL(ii,:));
    dip_tmp(ii,2)=std(DIP_DS8_0p2kPa_ERL(ii,:));    
end
fit_tmp = fit(t(216:end)',dip_tmp(216:end,1),'poly1');

s1e=shadedErrorBar(t,smooth(smooth(dip_tmp(:,1),7),15),smooth(smooth(dip_tmp(:,2),7),15));
% Set face and edge properties
set(s1e.edge,'Visible','off')
set(s1e.mainLine,'LineWidth',2,'LineStyle','--','color',[0.6350, 0.0780, 0.1840]	)
s1e.patch.FaceColor = [0.6350, 0.0780, 0.1840]	;

tmpx = 72:1:120;
tmpy = fit_tmp(tmpx);
tmp_plot = plot(tmpx,tmpy,'k');

legend([s1d.mainLine s1e.mainLine tmp_plot],'+DMSO','+ERL', 'DIP fit', 'Location', 'northwest') 


%%
%ds8 4

sb2=subplot(1,3,2);
title('4 kPa','fontsize',14)
sb2.Position(2)=0.2;
sb2.Position(4)=0.7;
axis([0 120 -1 4])
set(gca,'linewidth',1.5,'fontsize',12,'TickLength',[0.03,0.03],'TickDir','out')
grid on
a=gca;
set(a,'box','off','color','none')
b=axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
set(gca,'linewidth',1.5,'fontsize',12)
axes(a)
linkaxes([a b])

%dmso
dip_tmp=zeros(length(DIP_DS8_4kPa_DMSO),2);
for ii = 1:length(DIP_DS8_4kPa_DMSO)
    dip_tmp(ii,1)=mean(DIP_DS8_4kPa_DMSO(ii,:));
    dip_tmp(ii,2)=std(DIP_DS8_4kPa_DMSO(ii,:));
end
fit_tmp = fit(t(216:end)',dip_tmp(216:end,1),'poly1');

s2d=shadedErrorBar(t,smooth(smooth(dip_tmp(:,1),7),15),smooth(smooth(dip_tmp(:,2),7),15));
% Set face and edge properties
set(s2d.edge,'Visible','off')
set(s2d.mainLine,'LineWidth',2,'color',C(3,:))
s2d.patch.FaceColor = C(3,:);
hold on
tmpx = 72:1:120;
tmpy = fit_tmp(tmpx);
plot(tmpx,tmpy,'k');

%erl
dip_tmp=zeros(length(DIP_DS8_4kPa_ERL),2);
for ii = 1:length(DIP_DS8_4kPa_ERL)
    dip_tmp(ii,1)=mean(DIP_DS8_4kPa_ERL(ii,:));
    dip_tmp(ii,2)=std(DIP_DS8_4kPa_ERL(ii,:));
end
fit_tmp = fit(t(216:end)',dip_tmp(216:end,1),'poly1');

s2e=shadedErrorBar(t,smooth(smooth(dip_tmp(:,1),7),15),smooth(smooth(dip_tmp(:,2),7),15));
% Set face and edge properties
set(s2e.edge,'Visible','off')
set(s2e.mainLine,'LineWidth',2,'LineStyle','--','color',[0.6350, 0.0780, 0.1840]	)
s2e.patch.FaceColor = [0.6350, 0.0780, 0.1840]	;
plot((0:1:120),(zeros(1,121)),'color',[0,0,0]+.85,'LineWidth',2);
tmpx = 72:1:120;
tmpy = fit_tmp(tmpx);
plot(tmpx,tmpy,'k');

%%
%ds8 55

sb3=subplot(1,3,3);
title('55 kPa','fontsize',14)
sb3.Position(2)=0.2;
sb3.Position(4)=0.7;
axis([0 120 -1 4])
set(gca,'linewidth',1.5,'fontsize',12,'TickLength',[0.03,0.03],'TickDir','out')
grid on
a=gca;
set(a,'box','off','color','none')
b=axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
set(gca,'linewidth',1.5,'fontsize',12)
axes(a)
linkaxes([a b])

%dmso
dip_tmp=zeros(length(DIP_DS8_55kPa_DMSO),2);
for ii = 1:length(DIP_DS8_55kPa_DMSO)
    dip_tmp(ii,1)=mean(DIP_DS8_55kPa_DMSO(ii,:));
    dip_tmp(ii,2)=std(DIP_DS8_55kPa_DMSO(ii,:));
end
fit_tmp = fit(t3(216:end)',dip_tmp(216:end,1),'poly1');

s3d=shadedErrorBar(t3,smooth(smooth(dip_tmp(:,1),7),15),smooth(smooth(dip_tmp(:,2),7),15));
% Set face and edge properties
set(s3d.edge,'Visible','off')
set(s3d.mainLine,'LineWidth',2,'color',C(3,:))
s3d.patch.FaceColor = C(3,:);
hold on
tmpx = 72:1:120;
tmpy = fit_tmp(tmpx);
plot(tmpx,tmpy,'k');

%erl
dip_tmp=zeros(length(DIP_DS8_55kPa_ERL),2);
for ii = 1:length(DIP_DS8_55kPa_ERL)
    dip_tmp(ii,1)=mean(DIP_DS8_55kPa_ERL(ii,:));
    dip_tmp(ii,2)=std(DIP_DS8_55kPa_ERL(ii,:));
end
fit_tmp = fit(t3(216:end)',dip_tmp(216:end,1),'poly1');

s3e=shadedErrorBar(t3,smooth(smooth(dip_tmp(:,1),7),15),smooth(smooth(dip_tmp(:,2),7),15));
% Set face and edge properties
set(s3e.edge,'Visible','off')
set(s3e.mainLine,'LineWidth',2,'LineStyle','--','color',[0.6350, 0.0780, 0.1840]	)
s3e.patch.FaceColor = [0.6350, 0.0780, 0.1840]	;
plot((0:1:120),(zeros(1,121)),'color',[0,0,0]+.85,'LineWidth',2);
tmpx = 72:1:120;
tmpy = fit_tmp(tmpx);
plot(tmpx,tmpy,'k');

%%
all_ax=axes(fig1a,'visible','off'); 
all_ax.XLabel.Visible='on';
all_ax.YLabel.Visible='on';
ylabel(all_ax,'Population Doubling','fontsize',14);
xlabel(all_ax,'Time (Hours)','fontsize',14);
