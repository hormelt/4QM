function [] = fake_data_maker(NFrames,SigPtcle,stdStep,RelativeNoiseAmplitude,Version)

% fake_data_maker from Tommy
% Optimized for making fake data for Maria Kilfoil algorithm by Justin
% Puts out a full video of all frames, fills a folder fov1 with all the
% frames seperate and puts in a time file wih
%
% INPUT:
%   NFrames: Number of frames.
%   sigPtcle: Radius of the particle.
%   stdStep: The standard step size of the particle.
%   RelativeNoiseAmplitude: Adjusts the Noise level.
%   Version: A version number to avoid conflicts.

% VARIABLES
maxInt = 255;
dt = 1;
xbound = 512;
ybound = 512;
% FILENAMES
stdStepStr = strrep(num2str(stdStep, '%.2f'), '.', '_');
RelativeNoiseAmplitudeStr = strrep(num2str(RelativeNoiseAmplitude, '%.2f'), '.', '_');
Experiment = ['S' stdStepStr 'N' RelativeNoiseAmplitudeStr 'SIM' num2str(Version,'%02d')]
time_filename = ['data/' Experiment '/fov1_times'];
track_filename = ['data/' Experiment '/fov1_track'];
video_filename = ['data/' Experiment '/' Experiment '.tif'];
mkdir(['data/' Experiment])
mkdir(['data/' Experiment '/fov1'])
% DATA CONSTRUCTION
noise_amplitude = RelativeNoiseAmplitude*maxInt;

x = 1:xbound;
y = 1:ybound;

[x, y] = meshgrid(x,y);

x0 = 51:47:xbound;
y0 = 51:47:ybound;

[x0, y0] = meshgrid(x0,y0);

x0 = x0 + 15*(rand(size(x0))-0.5);
y0 = y0 + 15*(rand(size(y0))-0.5);
nptcles = length(x0(:));
Tracks = zeros(nptcles*NFrames,4);

for frame = 1:NFrames
    x1 = x0 + stdStep*(randn(size(x0)));
    y1 = y0 + stdStep*(randn(size(y0)));
    Tracks(1+(frame-1)*nptcles:(frame)*nptcles,1) = reshape(x1,nptcles,1);
    Tracks(1+(frame-1)*nptcles:(frame)*nptcles,2) = reshape(y1,nptcles,1);
    Tracks(1+(frame-1)*nptcles:(frame)*nptcles,3) = ones(1,nptcles)*frame;
    Tracks(1+(frame-1)*nptcles:(frame)*nptcles,4) = 1:nptcles;
    data = zeros(size(x));
    
    
    for ptcle = 1:nptcles      
        data = data + exp(-((x-x1(ptcle)).^2 + (y-y1(ptcle)).^2)/(2*SigPtcle^2))*maxInt;       
    end
    
    data = data + noise_amplitude*rand(size(data,1),size(data,2));
    data = data/(maxInt+noise_amplitude)*maxInt;
    
    frame_filename = ['data/' Experiment '/fov1/fov1_' num2str(frame, '%04d') '.tif'];
    imwrite(uint8(data), frame_filename ,'tiff','compression','none','writemode', 'append')
    imwrite(uint8(data), video_filename ,'tiff','compression','none','writemode', 'append')
    
end

time = transpose(0:dt:NFrames*dt);
save(time_filename, 'time' );
save(track_filename,'Tracks')