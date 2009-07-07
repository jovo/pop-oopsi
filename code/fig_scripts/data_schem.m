clear, clc, %close all
path        = '/Users/joshyv/Research/oopsi/meta-oopsi/data/mrsic-flogel/';
imname      = '(20081126_13_15_31)-_natural_reg1_nat_135umdepth';
tifname     = [path imname '.tif'];
fname       = 'data_schem';
matname     = ['/Users/joshyv/Research/oopsi/pop-oopsi/data/' fname '.mat'];

% set switches of things to do
LoadTif     = 0;
GetROI      = 1;
GetEphys    = 0;

%% get image data

if LoadTif == 1                                     % get whole movie
    MovInf  = imfinfo(tifname);                     % get number of frames
    Im.T  = 100; %numel(MovInf);                  % only alternate frames have functional data
    Im.h  = MovInf(1).Height;
    Im.w  = MovInf(1).Width;
    Im.Np = Im.w*Im.h;
    Im.fname = fname;

    Im.DataMat = sparse(Im.w*Im.h,Im.T);% initialize mat to store movie
    for j=1:Im.T
        X = imread(tifname,j);
        Im.DataMat(:,j)=(X(:));
    end
    Im.MeanFrame=mean(Im.DataMat,2);
    save(matname,'Im')
else
    load(matname)
end

%% select roi

if GetROI == 1
    figure(1); clf,
    imagesc(reshape(Im.MeanFrame,Im.h,Im.w)')
    title('select roi radius, double click when complete')
    set(gca,'BusyAction','queu','Interruptible','on');
    ellipse0=imellipse;
    wait(ellipse0);
    vertices0=getVertices(ellipse0);
    xdat=vertices0(:,1);
    ydat=vertices0(:,2);
    Im.x0 = min(xdat) + .5*(max(xdat)-min(xdat));
    Im.y0 = min(ydat) + .5*(max(ydat)-min(ydat));
    a = max(xdat)-Im.x0;
    b = max(ydat)-Im.y0;
    Im.radius0=mean([a b]);

    [pixmatx pixmaty] = meshgrid(1:Im.h,1:Im.w);
    Im.roi = (((pixmatx-Im.x0).^2 + (pixmaty-Im.y0).^2 )<= Im.radius0^2);
    Im.roi_edge = edge(uint8(Im.roi));
    save(matname,'Im')
else
    load(matname)
end

%% get spike data

if GetEphys ==1
    Ep       = paqread(paqname, 'info');
    Ep       = Ep.ObjInfo;
    channels= 1:length(Ep.Channel);
    Data    = paqread(paqname,'Channels',channels);
    for i=channels
        Ep.Channel(i).Data = Data(:,i);
    end

    j=1;
    for i=channels
        if Ep.Channel(i).ChannelName(1)=='V'
            Ep.spt{j}=GetSpikeTimes(Ep.Channel(i).Data,0.5);
            j=j+1;
        elseif strcmp(Ep.Channel(i).ChannelName,'CameraSync')
            CameraSync = Data(:,i);
        end
    end
    Ep.Nc = j-1;
    dC                  = diff(CameraSync);
    frame_onset_times   = find(dC>2);
    fake_frames         = length(frame_onset_times)-Im.T;
    frame_onset_times(1:fake_frames)=[];
    Im.dt   = median(diff(frame_onset_times))/Ep.SampleRate;
    
    Ep.n=zeros(Ep.Nc,Ep.SamplesAcquired);
    for i=1:Ep.Nc
        Ep.n(i,Ep.spt{i})=1;
    end

    Ep.n=zeros(Im.T,Ep.Nc);
    for i=1:Ep.Nc
        Ep.n(:,i) = SubSampleSpikeTrain(frame_onset_times,Ep.spt{i});
    end

    save(matname,'Ep','-append')
else
    load(matname)
end

%% plot ROI
% Pl = PlotParams;
nrows = 2;
ncols = 3;

roi2=Im.roi';
roi_edge2=Im.roi_edge';

ROI_im      = Im.MeanFrame+max(Im.MeanFrame)*roi_edge2(:);
weighted_ROI= Im.MeanFrame.*roi2(:);

figure(2); clf,
subplot(nrows,ncols,1);
imagesc(reshape(ROI_im,Im.h,Im.w)')
colormap(gray)
set(gca,'XTickLabel',[],'YTickLabel',[])

% print fig
wh=[7 5];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
FigName = 'ass'; %['../code/' fname];
print('-depsc',FigName)
print('-dpdf',FigName)
