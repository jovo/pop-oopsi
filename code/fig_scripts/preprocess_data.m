clear, clc, %close all
datapath = '/Users/joshyv/Research/oopsi/meta-oopsi/data/';
dataset = 3;
switch dataset
    case 1
        Im.path = [datapath 'rafa/tanya/051809/s1m1/'];
        Im.fname= 's1m1';
        Im.Fura = 1;
    case 2
        Im.path = [datapath 'mrsic-flogel/'];
        Im.fname= '20081126_13_05_43_orientation_Bruno_reg1_ori1_135umdepth';
        Im.Fura = 0;
    case 3
        Im.path = [datapath 'mrsic-flogel/'];
        Im.fname= '20081126_13_15_31_natural_reg1_nat_135umdepth';
        Im.Fura = 0;
end

% set switches of things to do
LoadTif     = 0;
GetROI      = 1;
GetF        = 1;

cd /Users/joshyv/Research/oopsi/pop-oopsi/code/
Im.tifname     = [Im.path Im.fname '.tif'];
Im.matname     = [Im.path Im.fname '.mat'];


%% get image data

if LoadTif == 1                                     % get whole movie
    MovInf  = imfinfo(Im.tifname);                     % get number of frames
    Im.T  = numel(MovInf);                  % only alternate frames have functional data
    Im.h  = MovInf(1).Height;
    Im.w  = MovInf(1).Width;
    Im.Np = Im.w*Im.h;
    Nframes  = 1000;

    DataMat = sparse(Im.w*Im.h,Nframes);% initialize mat to store movie
    i=1;
    for t=round(linspace(1,Im.T,Nframes))
        X = imread(Im.tifname,t);
        DataMat(:,i)=(X(:));
        i=i+1;
        if mod(i,10)==0, display(t), end
    end
    Im.MeanFrame=reshape(mean(DataMat,2),Im.h,Im.w);
    save(Im.matname,'Im')
end

%% select rois

if GetROI == 1
    if ~isfield(Im,'MeanFrame'), load(Im.matname); end
    figure(1); clf,

    Im.seg_frame=z1(Im.MeanFrame);
    imagesc(Im.seg_frame)
    %     set(gcf,'Position',[[Im.w 1/2.5]*2.5 [Im.w Im.h]*4])


    % manually determine roi radius
    title('select roi radius, double click when complete')
    ellipse0=imellipse; wait(ellipse0);
    vertices0=getVertices(ellipse0);
    xdat=vertices0(:,1);
    ydat=vertices0(:,2);
    x0 = min(xdat) + .5*(max(xdat)-min(xdat));
    y0 = min(ydat) + .5*(max(ydat)-min(ydat));
    a = max(xdat)-x0;
    b = max(ydat)-y0;
    Im.radius0=mean([a b]);

    % define the null roi
    [pixmatx pixmaty] = meshgrid(1:Im.w,1:Im.h);
    inellipse0 = (((pixmatx-x0).^2 + (pixmaty-y0).^2 )<= Im.radius0^2);

    % manually select roi centers
    i=0;
    button=0;
    Im.rois=0*inellipse0;
    Im.roi_edges=Im.rois;

    imagesc(Im.seg_frame)
    title('select roi centers, hit space when complete')
    while button ~= 32
        i=i+1;
        [x y button] = ginput(1);
        if button == 32, break, end
        Im.roi{i}.x0        = round(x);
        Im.roi{i}.y0        = round(y);
        Im.roi{i}.roi       = (((pixmatx-Im.roi{i}.x0).^2 + (pixmaty-Im.roi{i}.y0).^2 )<= Im.radius0^2);
        Im.roi{i}.edge      = edge(uint8(Im.roi{i}.roi));
        [ssub_i ssub_j]     = find(Im.roi{i}.roi);
        Im.roi{i}.sub       = [ssub_i ssub_j]; %[sub_i-x0+Im.roi{i}.x0, sub_j-y0+Im.roi{i}.y0];
        Im.roi{i}.ind       = sub2ind([Im.h Im.w],Im.roi{i}.sub(:,1),Im.roi{i}.sub(:,2));
        Im.roi{i}.Np        = numel(Im.roi{i}.ind);
        Im.roi{i}.F         = zeros(Im.roi{i}.Np,Im.T);
        Im.rois             = Im.rois + i*Im.roi{i}.roi;
        Im.roi_edges        = Im.roi_edges + i*Im.roi{i}.edge;

        Im.seg_frame(Im.roi{i}.ind)=1;
        imagesc(Im.seg_frame)
        axis off
    end
    Im.Nrois=i-1;
    save(Im.matname,'Im')
end

%% get F

if GetF == 1
    if ~isfield(Im,'Nrois'), load(Im.matname); end
    for t=1:Im.T
        X = imread(Im.tifname,t);
        for i=1:Im.Nrois
            Im.roi{i}.F(:,t)=X(Im.roi{i}.ind);
        end
        if mod(t,50)==0, display(t), end
    end

    Im.F = zeros(Im.Nrois,Im.T);
    for i=1:Im.Nrois
        Im.F(i,:) = mean(Im.roi{i}.F);
    end
    if Im.Fura==1, Im.F=-double(Im.F); else  Im.F=double(Im.F); end
    Im.F=Im.F-repmat(min(Im.F'),Im.T,1)';
    Im.F=Im.F./repmat(max(Im.F'),Im.T,1)';
    Im.F=Im.F*2^16;
    Im.F=uint16(Im.F);

    figure(3), clf,
    subplot(121), plot(Im.F')
    subplot(122), imagesc(Im.F)
    save(Im.matname,'Im')
end

%% extract good stuff and detrend

if ~isfield(Im,'F'), load(Im.matname); end
switch dataset
    case 1
        keyboard
        Ftemp = double(Im.F);
        Ftrunc = Ftemp(:,[1:1820 2050:3740 3950:7250 7550:8250 8400:end]);
        Fdetrend = detrend(Ftrunc')';
        Ftrends = Fdetrend-Ftrunc;
        mins = min(Ftrends');
        maxs = max(Ftrends');
        trends = 0*F;
        for k=1:Im.Nrois
            trends(k,:) = linspace(mins(k),maxs(k),Im.T);
        end
        Ftemp = Ftemp - trends;
        Ftemp = Ftemp-repmat(min(Ftemp'),Im.T,1)';
        Ftemp = Ftemp./repmat(max(Ftemp'),Im.T,1)';
        Ftemp = Ftemp*2^16;
        Ftemp = uint16(Ftemp);
    case 2
        keyboard
        tvec  = 1000:5000;
        Ftemp = double(Im.F(:,tvec));
        Ftemp = detrend(Ftemp);
        Ftemp = Ftemp-repmat(min(Ftemp'),length(tvec),1)';
        Ftemp = Ftemp./repmat(max(Ftemp'),length(tvec),1)';
        Ftemp = Ftemp*2^16;
        Ftemp = uint16(Ftemp);
    case 3
        keyboard
        L       = 2^nextpow2(Im.T);
        F_FFT   = fft(Im.F(i,:),L)/Im.T;
        F_FFT(1:20) = 0;
        F_iFFT  = ifft(F_FFT,Im.T);
        imagesc(F_iFFT);
end

for k=1:Im.Nrois
    F{k}=Ftemp(k,:);
end

save(['../data/' Im.fname],'F')


%% plot ROI
% Pl = PlotParams;
nrows = 2;
ncols = 2;

% roi2=Im.roi';
% roi_edge2=Im.roi_edge';
%
% ROI_im      = Im.MeanFrame+max(Im.MeanFrame)*roi_edge2(:);
% weighted_ROI= Im.MeanFrame.*roi2(:);

figure(2); clf,
subplot(nrows,ncols,1);
imagesc(Im.seg_frame)
colormap(gray)
set(gca,'XTickLabel',[],'YTickLabel',[])

subplot(nrows,ncols,2);
imagesc(Im.F)
colormap(gray)
% colorbar
set(gca,'XTickLabel',[],'YTickLabel',[])

subplot(nrows,ncols,nrows+1:nrows*ncols);
plot(Im.F(3,:))

% print fig
wh=[7 5];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
FigName = ['../figs/' Im.fname];
print('-depsc',FigName)
print('-dpdf',FigName)
