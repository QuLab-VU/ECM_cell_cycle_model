FAKs = 0:0.01:1;
o_p21 = [];
p21 = [];

%p.v_s1p21 = 0.7;
% p.F_n = 8
% p.K_p21fak = 5
p.K_p21fak2 = 0.72;
%p.v_s1p21_2 = 1;

for ii = 1:length(FAKs)
    FAK=FAKs(ii);
    o_p21(ii) = p.s*((p.v_s1p21/(1+p.K_p21fak*FAK^p.F_n)));
    p21(ii) = p.s*((p.v_s1p21^p.F_n/(p.K_p21fak2^p.F_n+FAK^p.F_n)));
end

figure; hold on
plot(cells{1,3},cells{1,2}(:,19),'DisplayName','BAX')
plot(cells{1,3},cells{1,2}(:,21),'DisplayName','BH3')
plot(cells{1,3},cells{1,2}(:,24),'DisplayName','BCL2')
plot(cells{1,3},cells{1,2}(:,23),'DisplayName','BAXm')


d_avg = [];
for ii=1:100
    tmp=[];
    for jj=1:20
        tmp=[tmp;d{ii,jj}];
    end
    d_avg(ii,:)=mean(tmp);
end



figure; hold on
t = 0:0.001*2000:1000;
lyap_erls = zeros(1,100);
for ii = 1:100
    subplot(10,10,ii);
    hold on
    [TF,S1,S2] = ischange(smooth(smooth(d_avg(ii,:))),'linear','Threshold',75);
    if max(TF) == 0
        P = polyfit(t,d_avg(ii,:),1);
        yfit = polyval(P,t);
        plot(t,d_avg(ii,:))
        plot(t,yfit)
    else
        ch_ind = find(d_avg(ii,:)<-6);
        P = polyfit(t(1:ch_ind(1)),d_avg(ii,1:ch_ind),1);
        yfit = polyval(P,t(1:ch_ind(1)));
        plot(t,d_avg(ii,:))
        plot(t(1:ch_ind(1)),yfit)
    end
    lyap_erls(ii) = P(1);
end

d_avg = [];
for ii=1:100
    tmp=[];
    for jj=1:20
        tmp=[tmp;d{ii,jj}];
    end
    d_avg(ii,:)=mean(tmp);
end

%figure; hold on
t = 0:0.001*2000:1000;
lyap_erls = zeros(1,100);
for ii = 1:100
    % subplot(10,10,ii);
    % hold on
    % [TF,S1,S2] = ischange(smooth(smooth(d_avg(ii,:))),'linear','Threshold',75);

    ch_ind = find(d_avg(ii,:)<-6);
    if isempty(ch_ind)
        P = polyfit(t,d_avg(ii,:),1);
        yfit = polyval(P,t);
    else
        P = polyfit(t(1:ch_ind(1)),d_avg(ii,1:ch_ind),1);
        yfit = polyval(P,t(1:ch_ind(1)));
    end
    % plot(t,d_avg(ii,:))
    % plot(t(1:ch_ind(1)),yfit)

    lyap_erls(ii) = P(1);
end

d_avg = [];
for ii=1:100
    tmp=[];
    for jj=1:20
        tmp=[tmp;d{ii,jj}];
    end
    d_avg(ii,:)=mean(tmp);
end

%figure; hold on
t = 0:0.001*2000:1000;
lyap_dmsos = zeros(1,100);
for ii = 1:100
    % subplot(10,10,ii);
    % hold on
    % [TF,S1,S2] = ischange(smooth(smooth(d_avg(ii,:))),'linear','Threshold',75);

    ch_ind = find(d_avg(ii,:)<-6);
    if isempty(ch_ind)
        P = polyfit(t,d_avg(ii,:),1);
        yfit = polyval(P,t);
    else
        P = polyfit(t(1:ch_ind(1)),d_avg(ii,1:ch_ind),1);
        yfit = polyval(P,t(1:ch_ind(1)));
    end
    % plot(t,d_avg(ii,:))
    % plot(t(1:ch_ind(1)),yfit)

    lyap_dmsos(ii) = P(1);
end