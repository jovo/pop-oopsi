clear, clc
load('~/Research/oopsi/meta-oopsi/data/tom/20081126_13_05_43_orientation_Bruno_reg1_ori1_135umdepth.mat')
%%
nrows=2;
ncols=2;
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

tvec=1000:2800;
xtick=[0:900:tvec(end)];
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

subplot(nrows,ncols,2)
plot(z1(detrend(double(Im.F(4,tvec)))),'k')
axis([0 max(tvec)-min(tvec) 0 1])
title('Low SNR','FontSize',fs)
set(gca,'YTick',[])
set(gca,'XTick',xtick,'XTickLabel',[])

subplot(nrows,ncols,3)
plot(z1(detrend(double(Im.F(3,tvec)))),'k')
axis([0 max(tvec)-min(tvec) 0 1])
set(gca,'XTick',xtick,'XTickLabel',[])
set(gca,'YTick',[])
ylab=ylabel([{'$\widehat{\mathbf{n}}$'}],'FontSize',fs,'Interpreter','latex');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'XTick',xtick,'XTickLabel',xticklabel,'FontSize',fs)

subplot(nrows,ncols,4)
plot(z1(detrend(double(Im.F(4,tvec)))),'k')
axis([0 max(tvec)-min(tvec) 0 1])
set(gca,'XTick',xtick,'XTickLabel',xticklabel,'FontSize',fs)
set(gca,'YTick',[])

%%

wh=[4 4];
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
FigName = '~/Research/oopsi/pop-oopsi/figs/example_traces';
print('-depsc',FigName)
print('-dpdf',FigName)

