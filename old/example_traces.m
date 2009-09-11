clear, clc
load('~/Research/oopsi/meta-oopsi/data/tom/20081126_13_05_43_orientation_Bruno_reg1_ori1_135umdepth.mat')

G{1}=Im.F(3,:); 
G{2}=Im.F(4,:); 

% F{1}=detrend(double(G{1})); 
% F{2}=detrend(double(G{2})); 
% 
% F{1}=F{1}-min(F{1})+mean(G{1});
% F{2}=F{2}-min(F{2})+mean(G{1});
F=G;

save('~/Research/oopsi/meta-oopsi/data/tom/example_traces','F','G')
%%
fname='~/Research/oopsi/meta-oopsi/data/tom/example_traces';
run_datafit(fname);

%%
load('~/Research/oopsi/meta-oopsi/data/tom/example_traces')

load('~/Research/oopsi/meta-oopsi/data/tom/example_trace1','I')
E{1}=I;

load('~/Research/oopsi/meta-oopsi/data/tom/example_trace2','I')
E{2}=I;
I=E;
clear E


%%
nrows=2;
ncols=1;
fs=12;
figure(1), clf

a=Im.MeanFrame;
% MeanFrame = (a-min(a(:)))/(max(a(:))-min(a(:)));
% roi_edges = Im.roi_edges;
% roi_edges(roi_edges>1)=1;
% 
% seg_frame = MeanFrame+roi_edges;
% subplot(nrows,ncols,[1 4])
% subimage(seg_frame)
% set(gca,'XTick',[0 30 60 90 120])
% colormap(gray)
% title('Mean Frame','FontSize',fs)

tvec=[1000:1900]+1000;
xtick=[0:300:tvec(end)];
fs=12;
xticklabel=xtick/30;
subplot(nrows,ncols,1)
plot(z1(detrend(double(Im.F(3,tvec)))),'k')
axis([0 max(tvec)-min(tvec) 0 1])
title('High SNR','FontSize',fs)
set(gca,'XTick',xtick,'XTickLabel',[])
set(gca,'YTick',[])
ylab=ylabel([{'$\mathbf{F}$'}],'FontSize',fs,'Interpreter','latex');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% set(gca,'XTick',xtick,'XTickLabel',xticklabel,'FontSize',fs)
% xlabel('Time (sec)')

subplot(nrows,ncols,2)
plot(z1(detrend(double(Im.F(4,tvec)))),'k')
axis([0 max(tvec)-min(tvec) 0 1])
title('Low SNR','FontSize',fs)
set(gca,'YTick',[])
set(gca,'XTick',xtick,'XTickLabel',[])
set(gca,'XTick',xtick,'XTickLabel',xticklabel,'FontSize',fs)
xlabel('Time (sec)')
ylab=ylabel([{'$\mathbf{F}$'}],'FontSize',fs,'Interpreter','latex');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% subplot(nrows,ncols,3)
% plot(I{1}.M.nbar(tvec),'k')
% axis([0 max(tvec)-min(tvec) 0 1])
% set(gca,'XTick',xtick,'XTickLabel',[])
% set(gca,'YTick',[])
% ylab=ylabel([{'$\widehat{\mathbf{n}}$'}],'FontSize',fs,'Interpreter','latex');
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% 
% subplot(nrows,ncols,4)
% plot(I{2}.M.nbar(tvec),'k')
% axis([0 max(tvec)-min(tvec) 0 1])
% set(gca,'XTick',xtick,'XTickLabel',xticklabel,'FontSize',fs)
% set(gca,'YTick',[])
% xlabel('Time (sec)')

%%

wh=[7 3];
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
FigName = '~/Research/oopsi/pop-oopsi/figs/example_traces';
print('-depsc',FigName)
print('-dpdf',FigName)

