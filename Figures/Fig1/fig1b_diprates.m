load('dip_rates.mat')

fig1b=figure;hold on;
fig1b.Color='w';
set(gca,'XScale','log','TickLength',[0.03,0.03])
axis([0.1 100 -0.02 0.05])
DIP_fits = zeros(12,2);

%%
[C]=cbrewer('div','RdBu' , 3);

for ii = 1:12
    DIP_fits(ii,1)=mean(all_DIPs{ii});
    DIP_fits(ii,2)=std(all_DIPs{ii})/length(all_DIPs{ii});
end   
    


%%
box on
set(gca,'linewidth',1.5,'fontsize',15)
x_pts = [0.2,4,55];
ds8_dmso_dips = [DIP_fits(1,1),DIP_fits(3,1),DIP_fits(5,1)];
ds8_dmso_err = [DIP_fits(1,2),DIP_fits(3,2),DIP_fits(5,2)];
p1 = errorbar(x_pts,ds8_dmso_dips,ds8_dmso_err,'k','LineWidth',1.5,'Marker','s','MarkerEdgeColor',C(3,:),...
    'MarkerFaceColor',C(3,:),'MarkerSize',10,'LineStyle','none');
plot(x_pts,ds8_dmso_dips,'Color',C(3,:),'LineWidth',2)

ds8_erl_dips = [DIP_fits(2,1),DIP_fits(4,1),DIP_fits(6,1)];
ds8_erl_err = [DIP_fits(2,2),DIP_fits(4,2),DIP_fits(6,2)];
p2=errorbar(x_pts,ds8_erl_dips,ds8_erl_err,'k','LineWidth',1.5,'Marker','o','MarkerEdgeColor',[0.6350, 0.0780, 0.1840],...
    'MarkerFaceColor',[0.6350, 0.0780, 0.1840],'MarkerSize',10,'LineStyle','none');
plot(x_pts,ds8_erl_dips,'Color',[0.6350, 0.0780, 0.1840],'LineWidth',2)

[C]=cbrewer('div','PRGn' , 3);
ds9_dmso_dips = [DIP_fits(7,1),DIP_fits(9,1),DIP_fits(11,1)];
ds9_dmso_err = [DIP_fits(7,2),DIP_fits(9,2),DIP_fits(11,2)];
p3=errorbar(x_pts,ds9_dmso_dips,ds9_dmso_err,'k','LineWidth',1.5,'Marker','d','MarkerEdgeColor',C(3,:),...
    'MarkerFaceColor',C(3,:),'MarkerSize',10,'LineStyle','none');
plot(x_pts,ds9_dmso_dips,'Color',C(3,:),'LineWidth',2)

ds9_erl_dips = [DIP_fits(8,1),DIP_fits(10,1),DIP_fits(12,1)];
ds9_erl_err = [DIP_fits(8,2),DIP_fits(10,2),DIP_fits(12,2)];
p4=errorbar(x_pts,ds9_erl_dips,ds9_erl_err,'k','LineWidth',1.5,'Marker','^','MarkerEdgeColor',[0.4940, 0.1840, 0.5560],...
    'MarkerFaceColor',[0.4940, 0.1840, 0.5560],'MarkerSize',10,'LineStyle','none');
plot(x_pts,ds9_erl_dips,'Color',[0.4940, 0.1840, 0.5560],'LineWidth',2)


legend([p1, p2, p3, p4],'DS8 +DMSO','DS8 +ERL','DS9 +DMSO','DS9 +ERL', 'location','northwest') 
ylabel('DIP Rate','fontsize',18)
xlabel('Stiffness (kPa)', 'fontsize', 18)