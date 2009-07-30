clear, clc
load('~/Research/oopsi/meta-oopsi/data/tom/20081126_13_05_43_orientation_Bruno_reg1_ori1_135umdepth.mat')
%%

a=Im.MeanFrame;
MeanFrame = (a-min(a(:)))/(max(a(:))-min(a(:)));
roi_edges = Im.roi_edges;
roi_edges(roi_edges>1)=1;

seg_frame = MeanFrame+roi_edges;
figure(1), clf
subplot(2,2,[1 3])
imagesc(seg_frame)
colormap(gray)
title('Mean Frame','FontSize',fs)

tvec=1000:2800;
xtick=[0:900:tvec(end)];
fs=12;
xticklabel=xtick/30;
subplot(2,2,2)
plot(detrend(double(Im.F(3,tvec))),'k')
axis('tight')
title('High SNR','FontSize',fs)
set(gca,'XTick',xtick,'XTickLabel',[])
set(gca,'YTick',[])
ylab=ylabel([{'F   '}; {'(a.u.)'}],'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

subplot(2,2,4)
plot(detrend(double(Im.F(4,tvec))),'k')
axis('tight')
title('Low SNR','FontSize',fs)
set(gca,'XTick',xtick,'XTickLabel',xticklabel,'FontSize',fs)
set(gca,'YTick',[])
ylab=ylabel([{'F   '}; {'(a.u.)'}]);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

wh = [7 4];
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
FigName = ['~/Research/oopsi/pop-oopsi/figs/data_example'];
print('-depsc',FigName)
print('-dpdf',FigName)
