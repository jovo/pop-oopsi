function PlotMatrix2(Cell,Sim,P,Phat)

%%
fig=figure(5); clf,
col = jet(Sim.Nc);          % define colors for mean

tstart=2/Sim.dt;
tdt=1/Sim.dt;
tend = min(tstart+8*tdt,Sim.T);
subplot(2,3,[1 2 3]), hold on
for i=1:Sim.Nc
    plot(((Cell{i}.F)./max(Cell{i}.F))+1,'Color',col(i,:));
end
axis('tight'),
xticks  = tstart:tdt:tend;               % XTick positions
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',xticks*Sim.dt)
ylabel('Fluorescence')
xlabel('Time (sec)')

clims(1)=min([P.omega(:)' Phat{1}.omega(:)' Phat{2}.omega(:)']);
clims(2)=max([P.omega(:)' Phat{1}.omega(:)' Phat{2}.omega(:)']);
subplot(234), imagesc(P.omega,clims), colormap(gray),
set(gca,'XTick',[1:Sim.Nc],'YTick',[1:Sim.Nc]) %colorbar
title('True matrix')
ylabel('Presynaptic'), xlabel('Postsynaptic')

subplot(235), imagesc(Phat{1}.omega,clims), %colorbar
set(gca,'XTick',[1:Sim.Nc],'YTick',[1:Sim.Nc]) %colorbar
title('Matrix from spikes')

subplot(236), imagesc(Phat{2}.omega,clims), %colorbar
set(gca,'XTick',[1:Sim.Nc],'YTick',[1:Sim.Nc]) %colorbar
title('Matrix from fluorescence')

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','SimConnector')
print('-dpdf','SimConnector')

end