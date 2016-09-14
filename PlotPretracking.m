function [] = PlotPretracking(data,b,tracks,FeatSize,NmPerPixel, ...
                              FileStub,PlotOption)

% PlotPretracking plots the imagedata, the bandpassed data together with
% the results of the pretracking of the particles.
%
% INPUT:
%     data: Collection of image frames.
%     b: Collection of bandpassed image frames.
%     tracks: The collection of particle tracks in the frames.
%     FeatSize: Full optical diameter of particle (pixels).
%     NmPerPixel: Actual pixel width (nm).
%     FileStub: Path to the image_stack to be analyzed. 
%     PlotOption: Option for what plotting function will be used.

if strcmp(PlotOption,'simple')
    PlotRawData(data,tracks,FeatSize,NmPerPixel,FileStub)
    
elseif strcmp(PlotOption,'bandpass')
    PlotBandpassedData(data,b,tracks,FeatSize,NmPerPixel,FileStub)
else
    disp('Warning: PlotOpt in Pretrack not recognized.') 
end  
end

function [] = PlotRawData(data,tracks,FeatSize,NmPerPixel,FileStub)
% PLOTRAWDATA plots the rawdata together with the centroids and puts out 
% frame to a movie file every 10 frames.

fig1 = figure;
whitebg(fig1);

for frame = 1:size(data,3)
    tempx = tracks(tracks(:,5)==frame,1);
    tempy = tracks(tracks(:,5)==frame,2);
    
    hold off
    imagesc(data(:,:,frame))
    colormap gray
    hold on
    scatter(tempx,tempy,'r')
    title(['\fontsize{16} Raw data [1 px =' num2str(NmPerPixel) ' nm]'])
    xlabel('x-axis (pixels)')
    ylabel('y-axis (pixels)')
    truesize
    if mod(frame,10) == 0
        f = getframe;
        imwrite(frame2im(f),[FileStub 'tracking_movie.tif'], ...
                'tiff','compression','none','writemode','append');
    end
end
close
end  
 
function [] = PlotBandpassedData(data,b,tracks,FeatSize,NmPerPixel,FileStub)
% PLOTBANDPASSEDDATA plots the rawdata and the bandpassed data together 
% with the centroids and the boundaries of the feature as defined.
% The function puts out a frame to a movie file every 10 frames.

fig1 = figure;
whitebg(fig1);
for frame = 1:size(data,3)   
    tempx = tracks(tracks(:,5)==frame,1);
    tempy = tracks(tracks(:,5)==frame,2);
    hold off
    ax1 = subplot(1,2,1);
    h1 = imagesc(data(:,:,frame));
    colormap gray
    hold(ax1,'on')
    h2 = scatter(ax1,tempx,tempy, 2*pi*(FeatSize/2)^2,'g');
    hold(ax1,'on')
    h3 = scatter(ax1,tempx,tempy,10,'r','filled');
    title(ax1, ['\fontsize{16} Raw data [1 px =' num2str(NmPerPixel) ' nm]'])
    xlabel(ax1, 'x-axis (pixels)')
    ylabel(ax1, 'y-axis (pixels)')
    ax2 = subplot(1,2,2); 
    h4 = imagesc(b(:,:,frame));
    colormap gray
    hold(ax2,'on')        
    h5 = scatter(ax2,tempx,tempy, 2*pi*(FeatSize/2)^2,'g');
    hold(ax2,'on')
    h6 = scatter(ax2,tempx,tempy,10,'r','filled');
    title(ax2, ['\fontsize{16} Bandpassed data [1 px =' num2str(NmPerPixel) ' nm]'])
    xlabel(ax2, 'x-axis (pixels)')
    ylabel(ax2, 'y-axis (pixels)')   
    truesize
    if mod(frame,10) == 0
        f = getframe(fig1);
        imwrite(frame2im(f),[FileStub 'tracking_movie.tif'], ...
                'tiff','compression','none','writemode','append');
    end
end
close
<<<<<<< HEAD
end
=======
end
>>>>>>> 4107f2ab22e3c274933e723d7e9e1c5b30102c15
