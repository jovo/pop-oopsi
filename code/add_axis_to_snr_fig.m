clc
open('../figs/FigureA6_recvar_SNR.fig')

h=get(gcf,'Children');

h1=get(h(3),'Children');

x1=get(h1(1),'XData'); y1=get(h1(1),'YData');

x2=get(h1(2),'XData'); y2=get(h1(2),'YData');

x3=get(h1(3),'XData'); y3=get(h1(3),'YData');

X1=x1;
YMatrix1=[y1; y2; y3];
%%
% Create figure
figure1 = figure('XVisual',...
    '0x24 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)');
% figure1=figure(1); clf
fs=18;
fs2=14;
% Create axes
axes1 = axes('Parent',figure1,'FontSize',fs,'FontName','Arial','XTick',0:10:80);
box('on');
hold('all');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',axes1,'MarkerSize',16,'Marker','.',...
    'LineWidth',2,...
    'Color',[0 0 0]);
set(plot1(1),'DisplayName','66Hz');
set(plot1(2),'LineStyle','--','DisplayName','33Hz');
set(plot1(3),'LineStyle',':','DisplayName','15Hz');

% Create xlabel
xlabel('Photon budget, Kph/neuron/frame','FontSize',fs,'FontName','Arial');

% Create ylabel
ylabel('$r^2$','FontSize',fs,'FontName','Arial','interpreter','latex');

axis([0 80 0 0.7])

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthWest');

%%
hold on
plot(2651.7/1000*ones(100,1),linspace(0,100,100),'k')
plot(18973/1000*ones(100,1),linspace(0,100,100),'k')
%%

ax2 = axes('Position',get(gca,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','r','YColor','k',...
    'XTick',linspace(0,1,9),...
    'XTickLabel',[0 4.8 5.9 6.7 7.5 7.8 8.1 8.4 8.7],...  %[4.8 5.9 6.7 7.5 7.8 8.1 8.4 8.7]
    'YTickLabel',[],...
    'TickLength',[0 0],...
    'fontsize',fs); %,...

xlabel(ax2,'Effective Signal to Noise Ratio','color','r','fontsize',fs)

wh=[11 8];
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
FigName = '~/Research/oopsi/pop-oopsi/figs/FigureA6_recvar_SNR';
print('-depsc',FigName)
print('-dpdf',FigName)