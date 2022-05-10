t1=(1:186)*(32/60);
t2 = (0:405)/3;
t5 = (0:360)/3;


DIP_DS9_0p2kPa_DMSO = importdata('DIP_DS9_0p2kPa_DMSO.txt');
DIP_DS9_0p2kPa_ERL = importdata('DIP_DS9_0p2kPa_ERL.txt');
DIP_DS9_4kPa_DMSO = importdata('DIP_DS9_4kPa_DMSO.txt');
DIP_DS9_4kPa_ERL = importdata('DIP_DS9_4kPa_ERL.txt');
DIP_DS9_55kPa_DMSO = importdata('DIP_DS9_55kPa_DMSO.txt');
DIP_DS9_55kPa_ERL = importdata('DIP_DS9_55kPa_ERL.txt');

fig1a=figure; hold on
fig1a.Color='w';
fig1a.Position = [100 100 1000 300];

%%
%ds9 0.2
[C]=cbrewer('div','PRGn' , 3);

sb1=subplot(1,3,1);
title('0.2 kPa','fontsize',18)
sb1.Position(2)=0.2;
sb1.Position(4)=0.7;
axis([0 120 -1 4])
set(gca,'linewidth',1.5,'fontsize',15,'TickLength',[0.03,0.03],'TickDir','out')
grid on
a=gca;
set(a,'box','off','color','none')
b=axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
set(gca,'linewidth',1.5,'fontsize',15)
axes(a)
linkaxes([a b])

%dmso
dip_tmp=zeros(length(DIP_DS9_0p2kPa_DMSO),2);
for ii = 1:length(DIP_DS9_0p2kPa_DMSO)
    dip_tmp(ii,1)=mean(DIP_DS9_0p2kPa_DMSO(ii,:));
    dip_tmp(ii,2)=std(DIP_DS9_0p2kPa_DMSO(ii,:));
end
fit_ds9_dmso_0p2 = fit(t2(216:end)',dip_tmp(216:end,1),'poly1');

s4d=shadedErrorBar(t2,smooth(smooth(dip_tmp(:,1),7),15),smooth(smooth(dip_tmp(:,2),7),15));
% Set face and edge properties
set(s4d.edge,'Visible','off')
set(s4d.mainLine,'LineWidth',2,'color',C(3,:))
s4d.patch.FaceColor = C(3,:);
hold on
tmpx = 72:1:120;
tmpy = fit_ds9_dmso_0p2(tmpx);
plot(tmpx,tmpy,'k');

%erl
dip_tmp=zeros(length(DIP_DS9_0p2kPa_ERL),2);
for ii = 1:length(DIP_DS9_0p2kPa_ERL)
    dip_tmp(ii,1)=mean(DIP_DS9_0p2kPa_ERL(ii,:));
    dip_tmp(ii,2)=std(DIP_DS9_0p2kPa_ERL(ii,:));
end
fit_ds9_erl_0p2 = fit(t2(216:end)',dip_tmp(216:end,1),'poly1');

s4e=shadedErrorBar(t2,smooth(smooth(dip_tmp(:,1),7),15),smooth(smooth(dip_tmp(:,2),7),15));
% Set face and edge properties
set(s4e.edge,'Visible','off')
set(s4e.mainLine,'LineWidth',2,'LineStyle','--','color',[0.4940, 0.1840, 0.5560])
s4e.patch.FaceColor = [0.4940, 0.1840, 0.5560];
plot((0:1:120),(zeros(1,121)),'color',[0,0,0]+.85,'LineWidth',2);
tmpx = 72:1:120;
tmpy = fit_ds9_erl_0p2(tmpx);
tmp_plot=plot(tmpx,tmpy,'k');

legend([s4d.mainLine s4e.mainLine tmp_plot],'+DMSO','+ERL', 'DIP fit', 'Location', 'northwest') 

%%
%ds8 4

sb2=subplot(1,3,2);
title('4 kPa','fontsize',18)
sb2.Position(2)=0.2;
sb2.Position(4)=0.7;
axis([0 120 -1 4])
set(gca,'linewidth',1.5,'fontsize',15,'TickLength',[0.03,0.03],'TickDir','out')
grid on
a=gca;
set(a,'box','off','color','none')
b=axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
set(gca,'linewidth',1.5,'fontsize',15)
axes(a)
linkaxes([a b])

%dmso
dip_tmp=zeros(length(DIP_DS9_4kPa_DMSO),2);
for ii = 1:length(DIP_DS9_4kPa_DMSO)
    dip_tmp(ii,1)=mean(DIP_DS9_4kPa_DMSO(ii,:));
    dip_tmp(ii,2)=std(DIP_DS9_4kPa_DMSO(ii,:));
end
fit_ds9_dmso_4 = fit(t5(216:end)',dip_tmp(216:end,1),'poly1');

s5d=shadedErrorBar(t5,smooth(smooth(dip_tmp(:,1),7),15),smooth(smooth(dip_tmp(:,2),7),15));
% Set face and edge properties
set(s5d.edge,'Visible','off')
set(s5d.mainLine,'LineWidth',2,'color',C(3,:))
s5d.patch.FaceColor = C(3,:);
hold on
tmpx = 72:1:120;
tmpy = fit_ds9_dmso_4(tmpx);
plot(tmpx,tmpy,'k');

%erl
dip_tmp=zeros(length(DIP_DS9_4kPa_ERL),2);
for ii = 1:length(DIP_DS9_4kPa_ERL)
    dip_tmp(ii,1)=mean(DIP_DS9_4kPa_ERL(ii,:));
    dip_tmp(ii,2)=std(DIP_DS9_4kPa_ERL(ii,:));
end
fit_ds9_erl_4 = fit(t5(216:end)',dip_tmp(216:end,1),'poly1');

s5e=shadedErrorBar(t5,smooth(smooth(dip_tmp(:,1),7),15),smooth(smooth(dip_tmp(:,2),7),15));
% Set face and edge properties
set(s5e.edge,'Visible','off')
set(s5e.mainLine,'LineWidth',2,'LineStyle','--','color',[0.4940, 0.1840, 0.5560])
s5e.patch.FaceColor = [0.4940, 0.1840, 0.5560];
plot((0:1:120),(zeros(1,121)),'color',[0,0,0]+.85,'LineWidth',2);
tmpx = 72:1:120;
tmpy = fit_ds9_erl_4(tmpx);
plot(tmpx,tmpy,'k');


%%
%ds8 55
sb3=subplot(1,3,3);
title('55 kPa','fontsize',18)
sb3.Position(2)=0.2;
sb3.Position(4)=0.7;
axis([0 120 -1 4])
set(gca,'linewidth',1.5,'fontsize',15,'TickLength',[0.03,0.03],'TickDir','out')
grid on
a=gca;
set(a,'box','off','color','none')
b=axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
set(gca,'linewidth',1.5,'fontsize',15)
axes(a)
linkaxes([a b])

%dmso
dip_tmp=zeros(length(DIP_DS9_55kPa_DMSO),2);
for ii = 1:length(DIP_DS9_55kPa_DMSO)
    dip_tmp(ii,1)=mean(DIP_DS9_55kPa_DMSO(ii,:));
    dip_tmp(ii,2)=std(DIP_DS9_55kPa_DMSO(ii,:));
end
fit_ds9_dmso_55 = fit(t2(216:end)',dip_tmp(216:end,1),'poly1');

s6d=shadedErrorBar(t2,smooth(smooth(dip_tmp(:,1),7),15),smooth(smooth(dip_tmp(:,2),7),15));
% Set face and edge properties
set(s6d.edge,'Visible','off')
set(s6d.mainLine,'LineWidth',2,'color',C(3,:))
s6d.patch.FaceColor = C(3,:);
hold on
tmpx = 72:1:120;
tmpy = fit_ds9_dmso_55(tmpx);
plot(tmpx,tmpy,'k');

%erl 
dip_tmp=zeros(length(DIP_DS9_55kPa_ERL),2);
for ii = 1:length(DIP_DS9_55kPa_ERL)
    dip_tmp(ii,1)=mean(DIP_DS9_55kPa_ERL(ii,:));
    dip_tmp(ii,2)=std(DIP_DS9_55kPa_ERL(ii,:));
end
fit_ds9_erl_55 = fit(t2(216:end)',dip_tmp(216:end,1),'poly1');

s6e=shadedErrorBar(t2,smooth(smooth(dip_tmp(:,1),7),15),smooth(smooth(dip_tmp(:,2),7),15));
% Set face and edge properties
set(s6e.edge,'Visible','off')
set(s6e.mainLine,'LineWidth',2,'LineStyle','--','color',[0.4940, 0.1840, 0.5560])
s6e.patch.FaceColor = [0.4940, 0.1840, 0.5560];
plot((0:1:120),(zeros(1,121)),'color',[0,0,0]+.85,'LineWidth',2);
tmpx = 72:1:120;
tmpy = fit_ds9_erl_55(tmpx);
plot(tmpx,tmpy,'k');

%%
all_ax=axes(fig1a,'visible','off'); 
all_ax.XLabel.Visible='on';
all_ax.YLabel.Visible='on';
ylabel(all_ax,'Population Doubling','fontsize',18);
xlabel(all_ax,'Time (Hours)','fontsize',18);