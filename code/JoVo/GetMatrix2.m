function omega = GetMatrix2(Nc,P)

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
