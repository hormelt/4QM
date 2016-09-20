% fake_data_maker from Tommy
% Optimized for making fake data for Maria Kilfoil algorithm by Justin

clear

% VARIABLES
nframes = 110;
std_step = 0.50;
sig_ptcle = 5;
max_int = 255;
dt = 0.01;
relative_noise_amplitude = 0.1;

% FILENAMES
std_step_str = strrep(num2str(std_step, '%.2f'), '.', '_');
relative_noise_amplitude_str = strrep(num2str(relative_noise_amplitude, '%.2f'), '.', '_');
experiment = ['S' std_step_str 'N' relative_noise_amplitude_str 'SIM']
time_filename = ['data/' experiment '/fov1_times'];
track_filename = ['data/' experiment '/fov1_track'];
video_filename = ['data/' experiment '/fov1_full.tif'];

% DATA CONSTRUCTION
noise_amplitude = relative_noise_amplitude*max_int;

x = 1:512;
y = 1:512;

[x y] = meshgrid(x,y);

x0 = 51:47:512;
y0 = 51:47:512;

[x0 y0] = meshgrid(x0,y0);

x0 = x0 + 15*(rand(size(x0))-0.5);
y0 = y0 + 15*(rand(size(y0))-0.5);

% imagesc(sqrt(x.^2 + y.^2));
% colormap gray
% hold on
% scatter(x0(:),y0(:));

nptcles = length(x0(:));
Tracks = zeros(nptcles*nframes,4);

for frame = 1:nframes
    frame
    x1 = x0 + std_step*(randn(size(x0)));
    y1 = y0 + std_step*(randn(size(y0)));
    Tracks(1+(frame-1)*nptcles:(frame)*nptcles,1) = reshape(x1,nptcles,1);
    Tracks(1+(frame-1)*nptcles:(frame)*nptcles,2) = reshape(y1,nptcles,1);
    Tracks(1+(frame-1)*nptcles:(frame)*nptcles,3) = ones(1,nptcles)*frame;
    Tracks(1+(frame-1)*nptcles:(frame)*nptcles,4) = 1:nptcles;
    data = zeros(size(x));
    
    for ptcle = 1:nptcles
        
        data = data + exp(-((x-x1(ptcle)).^2 + (y-y1(ptcle)).^2)/(2*sig_ptcle^2))*max_int;
        %         imagesc(data)
        %         colormap gray
        %         getframe
        
    end
    
    data = data + noise_amplitude*rand(size(data,1),size(data,2));
    data = data/(max_int+noise_amplitude)*max_int;
    
    frame_filename = ['data/' experiment '/fov1/fov1_' num2str(frame, '%04d') '.tif'];
    imwrite(uint8(data), frame_filename ,'tiff','compression','none','writemode', 'append')
    imwrite(uint8(data), video_filename ,'tiff','compression','none','writemode', 'append')
    
end

time = transpose(0:dt:nframes*dt);
save(time_filename, 'time' );
save(track_filename,'Tracks')