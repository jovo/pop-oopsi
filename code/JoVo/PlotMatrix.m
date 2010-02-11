function PlotMatrix(Cell,V,P,Phat)

fig=figure(5); clf,
col = [1 0 0; 0 0 1; 0 0.5 0; 1 0.5 0; 1 0 1];          % define colors for mean

tstart=2/V.dt;
tdt=1/V.dt;
tend = min(tstart+8*tdt,V.T);
subplot(2,3,[1 2 3]), hold on
for i=1:2
    plot(((Cell{i}.F(tstart:tend))./max(Cell{i}.F(tstart:tend)))+i,'Color',col(i,:));
end
axis('tight'),
xticks  = tstart:tdt:tend;               % XTick positions
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',xticks*V.dt)
ylabel('Fluorescence')
xlabel('Time (sec)')

if isempty(Phat{1})
    clims(1)=min([P.omega(:)' Phat{2}.omega(:)']); %Phat{2}.omega(:)']);
    clims(2)=max([P.omega(:)' Phat{2}.omega(:)']);% Phat{2}.omega(:)']);
else
    clims(1)=min([P.omega(:)' Phat{1}.omega(:)' Phat{2}.omega(:)']);
    clims(2)=max([P.omega(:)' Phat{1}.omega(:)' Phat{2}.omega(:)']);
end
subplot(234), imagesc(P.omega,clims), colormap(gray),
set(gca,'XTick',[1:V.Ncells],'YTick',[1:V.Ncells]) %colorbar
title('True matrix')
ylabel('Presynaptic'), xlabel('Postsynaptic')

subplot(235),
if ~isempty(Phat{1})
imagesc(Phat{1}.omega,clims), %colorbar
end
set(gca,'XTick',[1:V.Ncells],'YTick',[1:V.Ncells]) %colorbar
title('Matrix from spikes')
% ylabel('Presynaptic'), xlabel('Postsynaptic')

subplot(236), imagesc(Phat{2}.omega,clims), %colorbar
set(gca,'XTick',[1:V.Ncells],'YTick',[1:V.Ncells]) %colorbar
title('Matrix from fluorescence')
% ylabel('Presynaptic'), xlabel('Postsynaptic')

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','VConnector')
print('-dpdf','VConnector')

end