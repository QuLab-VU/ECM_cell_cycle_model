figure; hold on
t = 0:0.001*2000:1000;
lyap_erls = zeros(1,100);
for ii = 1:100
    subplot(10,10,ii);
    hold on
    [TF,S1,S2] = ischange(smooth(smooth(d{ii})),'linear','Threshold',75);
    if max(TF) == 0
        P = polyfit(t,d{ii},1);
        yfit = polyval(P,t);
        plot(t,d{ii})
        plot(t,yfit)
    else
        ch_ind = find(TF == 1);
        P = polyfit(t(1:ch_ind(1)),d{ii}(1:ch_ind(1)),1);
        yfit = polyval(P,t(1:ch_ind(1)));
        plot(t,d{ii})
        plot(t(1:ch_ind(1)),yfit)
    end
    lyap_erls(ii) = P(1);
end


figure; hold on
t = 0:0.001*2000:1000;
lyap_erls = zeros(1,100);
for ii = 1:100
    subplot(10,10,ii);
    if mod(ii,10) ~= 1
        set(gca,'YTickLabel',[])
    end
    if ii < 91
        set(gca,'XTickLabel',[])
    end
    ylim([-10 -2])    
    hold on
    [TF,S1,S2] = ischange(smooth(smooth(d{ii})),'linear','Threshold',75);
    if ii < 9
        P = polyfit(t,d{ii},1);
        yfit = polyval(P,t);
        plot(t,d{ii})
        plot(t,yfit)
    else
        ch_ind = find(d{ii}<-8.5);
        P = polyfit(t(1:ch_ind(1)),d{ii}(1:ch_ind(1)),1);
        yfit = polyval(P,t(1:ch_ind(1)));
        plot(t,d{ii})
        plot(t(1:ch_ind(1)),yfit)
    end
    lyap_erls(ii) = P(1);
end

figure;plot(lyap_erls)
ylim([-0.09 0.01])  

figure; hold on
t = 0:0.001*2000:1000;
lyap_dmso = zeros(1,100);
for ii = 1:100
    subplot(10,10,ii);
    if mod(ii,10) ~= 1
        set(gca,'YTickLabel',[])
    end
    if ii < 91
        set(gca,'XTickLabel',[])
    end
    ylim([-10 -2])    
    hold on
    [TF,S1,S2] = ischange(smooth(smooth(d{ii})),'linear','Threshold',75);
    % if ii < 9
    %     P = polyfit(t,d{ii},1);
    %     yfit = polyval(P,t);
    %     plot(t,d{ii})
    %     plot(t,yfit)
    % else
        ch_ind = find(d{ii}<-8.5);
        if isempty(ch_ind)
            P = polyfit(t,d{ii},1);
            yfit = polyval(P,t);
            plot(t,d{ii})
            plot(t,yfit)
        else
            P = polyfit(t(1:ch_ind(1)),d{ii}(1:ch_ind(1)),1);
            yfit = polyval(P,t(1:ch_ind(1)));
            plot(t,d{ii})
            plot(t(1:ch_ind(1)),yfit)
        end
    % end
    lyap_dmso(ii) = P(1);
end

figure;plot(lyap_dmso)
ylim([-0.09 0.01])  




figure; hold on
t = 0:0.001*2000:1000;
lyap_erls1 = zeros(1,100);
for ii = 1:100
    subplot(10,10,ii);
    hold on
    [TF,S1,S2] = ischange(smooth(smooth(d1{ii})),'linear','Threshold',75);
    if ii < 9
        P = polyfit(t,d1{ii},1);
        yfit = polyval(P,t);
        plot(t,d1{ii})
        plot(t,yfit)
    else
        ch_ind = find(d1{ii}<-7.5);
        P = polyfit(t(1:ch_ind(1)),d1{ii}(1:ch_ind(1)),1);
        yfit = polyval(P,t(1:ch_ind(1)));
        plot(t,d1{ii})
        plot(t(1:ch_ind(1)),yfit)
    end
    lyap_erls1(ii) = P(1);
end

figure;plot(smooth(lyap_erls1))



%%%%%%%%%%%%%%%%%%%%%%%%
lyap_erls = [];
for n = 1:6
    figure; hold on
    t = 0:0.001*2000:1000;
    
    for ii = 51:92
        subplot(10,10,ii);
        if mod(ii,10) ~= 1
            set(gca,'YTickLabel',[])
        end
        if ii < 91
            set(gca,'XTickLabel',[])
        end
        ylim([-10 -2])    
        hold on
        [TF,S1,S2] = ischange(smooth(smooth(d{ii,n})),'linear','Threshold',75);
        if ii < 9
            P = polyfit(t,d{ii,n},1);
            yfit = polyval(P,t);
            plot(t,d{ii,n})
            plot(t,yfit)
        else
            ch_ind = find(d{ii,n}<-6);
            P = polyfit(t(1:ch_ind(1)),d{ii,n}(1:ch_ind(1)),1);
            yfit = polyval(P,t(1:ch_ind(1)));
            plot(t,d{ii,n})
            plot(t(1:ch_ind(1)),yfit)
        end
        lyap_erls(ii,n) = P(1);
    end
end

figure; hold on
for ii = 1:6
    plot(lyap_erls(:,ii))
end

legend('loc1','loc2','loc3','loc4','loc5','loc6')

lyap_avg = mean(lyap_erls,dim,2);

%%%%%%%%%%%%%%%%%
lyap_dmso = [];
for n = 1:6
    figure; hold on
    t = 0:0.001*2000:1000;
    
    for ii = 1:11
        subplot(10,10,ii);
        if mod(ii,10) ~= 1
            set(gca,'YTickLabel',[])
        end
        if ii < 91
            set(gca,'XTickLabel',[])
        end
        ylim([-10 -2])    
        hold on
        [TF,S1,S2] = ischange(smooth(smooth(d{ii,n})),'linear','Threshold',75);
        if ii < 9
            P = polyfit(t,d{ii,n},1);
            yfit = polyval(P,t);
            plot(t,d{ii,n})
            plot(t,yfit)
        else
            ch_ind = find(d{ii,n}<-6);
            P = polyfit(t(1:ch_ind(1)),d{ii,n}(1:ch_ind(1)),1);
            yfit = polyval(P,t(1:ch_ind(1)));
            plot(t,d{ii,n})
            plot(t(1:ch_ind(1)),yfit)
        end
        lyap_dmso(ii,n) = P(1);
    end
end

figure; hold on
for ii = 1:6
    plot(lyap_dmso(:,ii))
end

legend('loc1','loc2','loc3','loc4','loc5','loc6')

lyap_avg = mean(lyap_dmso,dim,2);
