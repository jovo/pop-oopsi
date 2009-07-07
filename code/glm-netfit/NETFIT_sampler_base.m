%SCRIPT TO DRAW SAMPLE OF SPIKE HISTORIES FOR ESTIMATING W
% --- factorized/independent approximation, using PF

T=length(F{1}); 
for k=nrange{id_proc}
  H_loc=zeros(spkM,T);
  n_loc=false(spkM,T);
  
  n_loc(:,min(T,ceil(n_GT{k}/netSim.K)))=1;
  H_loc(:,min(T,ceil(n_GT{k}/netSim.K)))=1;
  
  n{k}=n_loc;
  H{k}=H_loc;
end
  
save(fname,'n','H','spkM');

%clean the mess
clear T H_loc n_loc