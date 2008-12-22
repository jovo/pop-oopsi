function PlotPop_v1_0(Sim,P,omega,omegahat)

Phat.omega=zeros(Sim.Nc);
for i=1:Sim.Nc
    Phat.omega(i,i)=P{i}.omega;
    Pre=1:Sim.Nc;                                   % generate list of presynaptic neurons
    Pre(Pre==i)=[];                                 % remove self
    k=0;                                            % counter of dimension
    for j=Pre
        k=k+1;                                      % generate input to neuron based on posterior mean spike train from neuron j
        Phat.omega(i,j)=P{i}.k(k+1);
    end
end
[omega round(Phat.omega*10)/10]

% figure(2), clf, nrows=Sim.Nc;
% for i=1:Sim.Nc
%     subplot(nrows,1,i), hold on,
%     stem(S(i).n,'LineStyle','none','Color','k'),
%     stem(I{i}.M.nbar,'Marker','none','Color',col(1,:))
%     ylabel(num2str(i)) % axis('tight')
% end

% %     if i==1
% figure(3), nrows=Sim.Nc;
% subplot(nrows,1,1), hold on, plot(S(1).h+1,'k'), plot(I{1}.S.h(1,:,1)+1)
% for j=2:Sim.Nc
%     subplot(nrows,1,j), hold on, plot(S(j).h+1,'k'), plot(Tim.x(j,:)+1)
% end
% %     end


figure(4), clf,
clims(1)=min([omega(:)' Phat.omega(:)' omegahat(:)']);
clims(2)=max([omega(:)' Phat.omega(:)' omegahat(:)']);
subplot(131), imagesc(omega,clims), colormap(gray), %colorbar
subplot(132), imagesc(omegahat,clims), %colorbar
subplot(133), imagesc(Phat.omega,clims), %colorbar

drawnow
end