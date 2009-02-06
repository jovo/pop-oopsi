clf, clc, clear
fig=figure(5); 
col = [0 0 0; .75 .75 .75; 0 0 0; 1 0.5 0; 1 0 1];          % define colors for mean
ls{1} = '-'; ls{2} =  '-'; ls{3} = '-.';
fs = 14;
fs2 = 12;

Ts=[29999 30000 4];

fnames{1}=['SimConnector', num2str(Ts(1)), '_10_1.mat'];
fnames{2}=['SimConnector', num2str(Ts(2)), '_10_1.mat'];
fnames{3}=['SimConnector', num2str(Ts(3)), '_10_1.mat'];

ylabs{1}='~2000 sp/neuron';
ylabs{2}='~2000 sp/neuron';
ylabs{3}='~600 sp/neuron';

ticks{2}=1:5;
ticks{1}=2:2:10;
ticks{3}=4:4:20;

load(fnames{2})
tstart=2/Sim.dt;
tdt=1/Sim.dt;
tend = tstart+4*tdt;

Ts([1 3])=[];
nrows=numel(Ts)+1;
subplot(nrows,3,[1 2 3]), hold on, if ~exist('R'), R=S; end
for i=1:2
    load(fnames{i})
%     if i==1, R=S; end
    plot(z1(R(i).F(tstart:tend))+1,'Color',col(i,:),'LineStyle',ls{i});
    bar(R(i).n(tstart:tend),'EdgeColor',col(i,:),'FaceColor',col(i,:))
end

axis('tight'),
xticks  = tstart:tdt:tend;               % XTick positions
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks-tstart,'XTickLabel',(xticks-tstart)*Sim.dt)
ylabel([{'Fluorescence'};{'and spikes'}],'FontSize',fs)
xlabel('Time (sec)','FontSize',fs2)

clims=[inf -inf];
for i=1:numel(Ts)
    load(fnames{i})
    clims(1)=min([clims(1) P.omega(:)' Phat{1}.omega(:)' Phat{2}.omega(:)']);
    clims(2)=max([clims(2) P.omega(:)' Phat{1}.omega(:)' Phat{2}.omega(:)']);
end


for i=1:numel(Ts)
%     load(fnames{i}); if i==1; R=S;  end
    subplot(nrows,3,3*i+1), imagesc(P.omega,clims), colormap(gray),
    set(gca,'XTick',[],'YTick',[])
    if i==nrows-1, title('True matrix','FontSize',fs2), end
    ylab=ylabel(ylabs{i},'FontSize',fs);
    if i==1
        ylabel('Presynaptic','FontSize',fs), 
        xlabel('Postsynaptic','FontSize',fs), 
    end
    set(gca,'XTick',ticks{i},'YTick',ticks{i})

    subplot(nrows,3,3*i+2), imagesc(Phat{1}.omega,clims), %colorbar
    set(gca,'XTick',[],'YTick',[])
    if i==nrows-1, title('Matrix from spikes','FontSize',fs2), end
    % ylabel('Presynaptic'), xlabel('Postsynaptic')

    subplot(nrows,3,3*i+3), imagesc(Phat{2}.omega,clims), %colorbar
    set(gca,'XTick',[],'YTick',[])
    if i==nrows-1, 
        title('Matrix from fluorescence','FontSize',fs2), 
    end
    
    
    % ylabel('Presynaptic'), xlabel('Postsynaptic')
    nmean=0;
    for j=1:Sim.Nc,
        nmean=nmean+sum(R(j).n);
    end
    nmean=nmean/Sim.Nc;
    disp(['nmean of row ' num2str(i) ' = ' num2str(nmean)])
end

%% print fig
wh=[7 4.5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','SimConnector2')
print('-dpdf','SimConnector2')
