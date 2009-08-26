%%
% clc
% open('../figs/FigureA6_recvar_SNR.fig')
% 
% h=get(gcf,'Children');
% 
% h1=get(h(2),'Children');
% 
% x1=get(h1(1),'XData'); y1=get(h1(1),'YData');
% 
% x2=get(h1(2),'XData'); y2=get(h1(2),'YData');
% 
% x3=get(h1(3),'XData'); y3=get(h1(3),'YData');
% 
% X1=x1;
% YMatrix1=[y1; y2; y3];

%%

pb = [1  5  10  20  40  80];
eSNR=[2.2  4.8   6.4   8.3   10.3   12];
figure(2), clf, plot(pb,eSNR,'.-'), grid on
%%
close all
open('../figs/FigureA6_recvar_SNR.fig')
h=get(gcf,'Children');

fs=10;

set(h(2),'fontsize',fs)
set(h(1),'fontsize',fs)
xlabel('Photon budget (Kph/neuron/frame)','fontsize',fs)
xlabel('Photon budget (Kph/neuron/frame)','fontsize',fs)
ylabel('$r^2$','interpreter','latex','fontsize',fs)

set(gca,'OuterPosition',[0 0 1 .93])

ticklength=get(gca,'ticklength');
ticklength(2)=0;
% set(gca,'TickLength',[0 0]);

grid on

hold on
plot(2*ones(10,1),linspace(0,.7,10),'k','linewidth',1.5)
plot(27*ones(10,1),linspace(0,.7,10),'k','linewidth',1.5)

r=[3 2.23;
  6 8.75;
  9 27;
  10 37;
  11 56.5'
  12 80];

ax2 = axes('Position',get(gca,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k',...
    'XTick',r(:,2)/max(r(:,2)),...
    'XTickLabel',r(:,1),...  %[4.8 5.9 6.7 7.5 7.8 8.1 8.4 8.7]
    'YTick',linspace(0,1,8),...
    'YTickLabel',[],...
    'TickLength',ticklength,...
    'fontsize',fs); %,...
  
xlabel(ax2,'Effective Signal to Noise Ratio (unitless)','color','k','fontsize',fs)

wh=[6 4];
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
FigName = '~/Research/oopsi/pop-oopsi/figs/FigureA6_recvar_SNRb';
print('-depsc',FigName)
print('-dpdf',FigName)


%%
% % Create figure
% figure1 = figure('XVisual',...
%     '0x24 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)');
% % figure1=figure(1); clf
% fs2=14;
% % Create axes
% axes1 = axes('Parent',figure1,'FontSize',fs,'FontName','Arial','XTick',0:10:80);
% box('on');
% hold('all');
% 
% % Create multiple lines using matrix input to plot
% plot1 = plot(X1,YMatrix1,'Parent',axes1,'MarkerSize',16,'Marker','.',...
%     'LineWidth',2,...
%     'Color',[0 0 0]);
% set(plot1(1),'DisplayName','66Hz');
% set(plot1(2),'LineStyle','--','DisplayName','33Hz');
% set(plot1(3),'LineStyle',':','DisplayName','15Hz');
% 
% % Create xlabel
% xlabel('Photon budget, Kph/neuron/frame','FontSize',fs,'FontName','Arial');
% 
% % Create ylabel
% ylabel('$r^2$','FontSize',fs,'FontName','Arial','interpreter','latex');
% 
% axis([0 80 0 0.7])
% 
% % Create legend
% legend1 = legend(axes1,'show');
% set(legend1,'Location','EastOutside');
% set(gca,'OuterPosition',[0 0 1 1])
% 
% hold on
% plot(2651.7/1000*ones(100,1),linspace(0,100,100),'k')
% plot(18973/1000*ones(100,1),linspace(0,100,100),'k')
% 
% ax2 = axes('Position',get(gca,'Position'),...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none',...
%     'XColor','r','YColor','k',...
%     'XTick',linspace(0,1,9),...
%     'XTickLabel',[0 4.8 5.9 6.7 7.5 7.8 8.1 8.4 8.7],...  %[4.8 5.9 6.7 7.5 7.8 8.1 8.4 8.7]
%     'YTickLabel',[],...
%     'TickLength',[0 0],...
%     'fontsize',fs); %,...
% 
% xlabel(ax2,'Effective Signal to Noise Ratio','color','r','fontsize',fs)
% %%
% 
% h=get(gcf,'Children');
% outerposition=[0 0 .9 1];
% set(h(1),'OuterPosition',[0 0 .9 1])
% % set(h(2),'OuterPosition',[0 0 1 1])
% set(h(3),'OuterPosition',[0 0 .8 1])
% 
% %%
% wh=[6 4]*1.4;
% set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
% FigName = '~/Research/oopsi/pop-oopsi/figs/FigureA6_recvar_SNR';
% print('-depsc',FigName)
% print('-dpdf',FigName)