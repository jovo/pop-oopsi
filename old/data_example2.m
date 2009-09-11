% load('~/Research/oopsi/smc-oopsi/data/BurstData.mat')

%%
Algs=7;
fig=figure(3); clf,
nrows=2;
ncols=2;
gray  = [.5 .5 .5];              % define gray color
col     = [1 0 0; 0 .5 0];          % define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;   % define colors for std
inter = 'latex';                     % interpreter for axis labels
xlims = [5 Sim.T-Sim.freq+1];       % xmin and xmax for current plot
fs=12;                              % font size
ms=4;                              % marker size for real spike
sw=1.5;                             % spike width
lw=2;                               % line width
I{7}.name=[{'Nonlinear Observation'}; {'PFS Spike Inference'}];
I{7}.name=[{'$n$'}; ];
spt=find(D{j}.n==1); spt=spt+1; Sim.n=zeros(size(Sim.n)); Sim.n(spt)=1;
% make xticks
Nsec = floor(D{j}.T_o*D{j}.dt_o);
secs = zeros(1,Nsec-2);
for i=1:Nsec-2
    secs(i) = find(D{j}.FrameStartTime>i,1);
end


col   = [1 0 0; 0.2 0.2 1];            % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;     % define colors for std


% plot real data
i=1; 
subplot(nrows,ncols,i), hold on
plot(tvec_o,z1(D{j}.F(tvec_o)),'-k','LineWidth',1,'MarkerSize',ms);
% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
% ylab=ylabel([{'in vitro'}; {'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
ylab=ylabel([{'$F$'}; ],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',secs,'XTickLabel',[])
axis([xlims 0 1.1*max(D{j}.F(1:xlims(2)))])

subplot(nrows,ncols,i+1), hold on

S=Hill_v1(I{7}.P,I{7}.M.Cbar);
noisyF = I{7}.P.alpha*S+I{7}.P.beta + sqrt(I{7}.P.gamma*S*1000+I{7}.P.zeta).*randn(size(F))';
plot(tvec_o,z1(noisyF(tvec_o)),'-k','LineWidth',1,'MarkerSize',ms);
% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
% ylab=ylabel([{'in vitro'}; {'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
ylab=ylabel([{'$F$'}; ],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',secs,'XTickLabel',[])
axis([xlims 0 1.1*max(D{j}.F(1:xlims(2)))])


i=i+2;
subplot(nrows,ncols,i), hold on,
BarVar=I{m}.M.nbar+I{m}.M.nvar; BarVar(BarVar>1)=1;
spts=find(BarVar>1e-3);
stem(spts,BarVar(spts),'Marker','none','LineWidth',sw,'Color',gray);
spts=find(I{m}.M.nbar>1e-3);
stem(spts,I{m}.M.nbar(spts),'Marker','none','LineWidth',sw,'Color','k')
stem(spt,1.1*Sim.n(spt),'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
axis([xlims 0 1])


hold off,
ylab=ylabel(I{m}.name,'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',secs,'XTickLabel',[])

% % plot calcium
% i=i+1;
% subplot(nrows,ncols,i+2), hold on
% C = I{m}.M.Cbar/I{m}.P.k_d;
% hfill=fill([1:Sim.T Sim.T:-1:1],[I{m}.M.Cptiles(1,:) I{m}.M.Cptiles(end,Sim.T:-1:1)]/I{m}.P.k_d,ccol(2,:));
% set(hfill,'edgecolor','k')
% plot(C(1:end-1),'Color','k','LineWidth',2)
% set(gca,'YTick',1,'YTickLabel',[])
% ylab=ylabel([{'Nonlinear Observation'}; {'PFS [Ca^{2+}] Inference'}],'Interpreter',inter,'FontSize',fs);
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% axis([xlims(1) xlims(2)-2 min(C) 1.1])


set(gca,'XTick',secs,'XTickLabel',round((secs-xlims(1))*Sim.dt),'FontSize',fs)
xlabel('Time (sec)','FontSize',fs)

