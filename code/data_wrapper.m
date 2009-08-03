clear, clc
load('~/Research/oopsi/meta-oopsi/data/vincent/timecourses.mat')
%%
G=tcs.raw;
G=detrend(G);

%%
siz=size(G);
for k=1:siz(2)
    H{k}=G(:,k);
    if min(H{k})<0
        H{k}=H{k}-min(H{k});
    end
end
%%
fname='~/Research/oopsi/meta-oopsi/data/vincent/raw';

%%
k=5;
F{1}=H{k};
save([fname, '_' num2str(k)] ,'F')

%%
run_datafit('~/Research/oopsi/meta-oopsi/data/vincent/raw_5.mat')