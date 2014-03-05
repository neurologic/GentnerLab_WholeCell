function [out_struct,hfig] = MetaResponseAnal_SpikeCluster(input_struct)  

allfields = fieldnames(input_struct);
for ifield = 1:size(allfields,1)
    s = [allfields{ifield} ' = input_struct.' allfields{ifield} ';'];
    eval(s)
end

spknorm = [];
spkind = 1;
for ispk = 1:size(shape,1)
    
    prespk = shape(ispk,225:276);
    
    if ~isempty(find(prespk))
    threshind = find(prespk == initVm(ispk));
%     if isempty(threshind)
%         threshind = crossing(prespk,[],initVm(itrial)) + 1;
%     end
    startind = (225+threshind-1);
    spktrunk = shape(ispk,startind:(startind + 70));
    spktrunk = spktrunk - initVm(ispk);
    spknorm(ispk,:) = spktrunk;
    
    %spike height (Vm at threshold to peak Vm)
    spkheight(ispk) = max(spktrunk);
    
    %get width at half-max height
    halfH = max(spktrunk)/2;
    fitfunc = spktrunk - halfH;
    xtime = [1:size(spktrunk,2)]*dt;
    [ind,t,s] = crossing(fitfunc,xtime,0,'linear');
    if max(size(t)) == 1
        spkwid(ispk) = NaN;
    end
    if max(size(t)) > 1
    spkwid(ispk) = min(diff(t)) * 1000;
    end
    %get max dv/dt during rise time in volts/sec
    spkdvdt(ispk) = round(max(diff(spktrunk/1000'))/dt);
%     spkdvdt(itrial) = (spkheight(itrial)/1000) / xtime(min(find(spktrunk == max(spktrunk)))) - xtime(1);

spkind = spkind +1;
    end
end

X = [spkwid',spkheight',spkdvdt'];
[idx,ctrs] = kmeans(X,2);

usegrp = find(ctrs(:,2) == max(ctrs(:,2)));
out_struct.useSpks = spknorm(idx == usegrp,:);
out_struct.spkwid = spkwid(1,idx == usegrp);
out_struct.spkheight = spkheight(1,idx == usegrp);
out_struct.spkdvdt = spkdvdt(1,idx == usegrp);

hfig(1) = figure
scatter3(X(idx==1,1),X(idx==1,2),X(idx==1,3),12,'r','fill')
hold on
scatter3(X(idx==2,1),X(idx==2,2),X(idx==2,3),12,'b','fill')
plot3(ctrs(:,1),ctrs(:,2),ctrs(:,3),'kx',...
     'MarkerSize',12,'LineWidth',2)
plot3(ctrs(:,1),ctrs(:,2),ctrs(:,3),'ko',...
     'MarkerSize',12,'LineWidth',2)
xlabel('spkwid msec')
ylabel('spkheight')
zlabel('spkdvdt V/sec')
title(['n = ' num2str(size(idx,1))]);

hfig(2) = figure;
spks1 = spknorm(idx ==1,:);
line(xtime,spks1','color',[0.5 0.5 0.5])
line(xtime,mean(spks1)','color','k','LineWidth',3)
title(['n = ' num2str(size(find(idx == 1),1))])

hfig(3) = figure;
spks2 = spknorm(idx ==2, :);
line(xtime,spks2','color',[0.5 0.5 0.5])
line(xtime,mean(spks2)','color','k','LineWidth',3)
title(['n = ' num2str(size(find(idx == 2),1))])

hfig(4) = figure;
hold on
subplot(3,1,1)
hist(spkwid)
title('spike width (msec)')
subplot(3,1,2)
hist(spkheight)
title('spike height (mV)')
subplot(3,1,3)
hist(spkdvdt)
title('max dvdt in initial rise (V/sec)')
set(hfig(4),'Position',[98   131   463   663])


