%SCRIPT TO DRAW SAMPLE OF SPIKE HISTORIES FOR ESTIMATING W
% --- factorized/independent approximation, using PF

spkM=1;
T_loc=length(F{1}); 
for k=1:N
  out_H=zeros(1,T_loc);
  out_n=false(1,T_loc);
  
  out_n(1,ceil(n_GT{k}/netSim.K))=1;
  out_H(1,ceil(n_GT{k}/netSim.K))=1;
  
  n{k}=out_n;
  H{k}=out_H;
end
  