% load frame-0717-66FR-result_25_600
% load frame-0717-66FR-sample-25

%%
clc
maxlag=10;
cc=zeros(25,25);
xcov2=zeros(25,25);
for i=1:25
    n1=double(n{i});
    for j=1:25
        n2=double(n{j});
        xxcov=xcov(n1,n2,maxlag,'coeff');
        if j<i
            cc_temp=corrcoef(n1,n2);
            cc(i,j)=cc_temp(2);  
        end
        xcov2(i,j)=xxcov(maxlag+1);
    end
end
for i=1:25, xcov2(i,i)=0; end

%%
figure(1), clf
scatter(xcov2(:), result.GT_W(:))
hold all,
scatter(cc(:), result.GT_W(:))
legend('xcov2','cc')
print('corr_v_w','-dpdf')