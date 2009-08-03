clear, clc
load('~/Research/oopsi/meta-oopsi/data/vincent/timecourses.mat')
%%
G=tcs.raw;
G=detrend(G);

%%
siz=size(G);
for k=1:siz(2)
%     G=
end
% run_datafit(fname)