function omegahat = plot_poopsi(V,smc,omega)

% given new inference for each neuron, update parameters

for i=1:V.Ncells
    EE{i}    = smc{i}.P;
end
    
omegahat = GetMatrix(V.Ncells,EE);

clims(1)=min([omega(:)' omegahat(:)']);
clims(2)=max([omega(:)' omegahat(:)']);

figure(1), clf
subplot(121), 
imagesc(omega,clims), colormap(gray),
set(gca,'XTick',[1:V.Ncells],'YTick',[1:V.Ncells]) %colorbar
title('True matrix')
ylabel('Presynaptic'), xlabel('Postsynaptic')

subplot(122), 
imagesc(omegahat,clims), %colorbar
set(gca,'XTick',[1:V.Ncells],'YTick',[1:V.Ncells]) %colorbar
title('Matrix from fluorescence')

% print fig
wh=[7 5];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
figname='../../figs/sim1';
print('-depsc',figname)
print('-dpdf',figname)
saveas(gcf,figname)

end

function omega = GetMatrix(Nc,P)

omega=zeros(Nc);
for i=1:Nc
    omega(i,i)=P{i}.omega;
    Pre=1:Nc;                                       % generate list of presynaptic neurons
    Pre(Pre==i)=[];                                 % remove self
    k=0;                                            % counter of dimension
    for j=Pre
        k=k+1;                                      % generate input to neuron based on posterior mean spike train from neuron j
        omega(i,j)=P{i}.k(k+1);
    end
end

end