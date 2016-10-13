% recall that column 3 is time and column 4 is particle ID.

function res = calcMSD(particle_tracks,nm_per_pixel,collective_motion_flag)

%
particle_tracks2 = zeros([size(particle_tracks,1) size(particle_tracks,2)+2]);
particle_tracks2(:,1:4)=particle_tracks;

%
% accumulate instantaneous displacements; store in columnes 5&6 of
% particle_tracks2.
%
trackids = unique(particle_tracks(:,4));
for j = 1:numel(trackids) % can't just count from min to max in case of missed tracks. I've corrected this error at least 5 times. --Grouchy TH
    
    ptcleid = trackids(j);
    ptcle_logic = particle_tracks(:,4)==ptcleid;
    
    if sum(ptcle_logic)~=0
        spt = particle_tracks(ptcle_logic,:);
        disp = diff(spt(:,1:2),1);
        disp = [[0 0]; disp];
        particle_tracks2(ptcle_logic,5:6)=disp;
    end
    
end

particle_tracks3 = particle_tracks2;

%
% subtract collective motion at each instant
%

if collective_motion_flag == 1;
    
    for frame = 2:(max(particle_tracks3(:,3)))
        %
        subtracks = particle_tracks3(particle_tracks3(:,3)==frame,:);
        
        all_disp = mean(subtracks(:,5:6),1);
        particle_tracks3(particle_tracks3(:,3)==frame,1)=particle_tracks3(particle_tracks3(:,3)==frame,1)-all_disp(1);
        particle_tracks3(particle_tracks3(:,3)==frame,2)=particle_tracks3(particle_tracks3(:,3)==frame,2)-all_disp(2);
        
    end
    
end

%
% compute MSD
%
ptcle_index = 0;
newtracks = unique(particle_tracks3(:,4));
for j = 1:numel(newtracks) % ... --TH
    
    ptcleid = newtracks(j);
    
    spt = particle_tracks3(particle_tracks3(:,4)==ptcleid,:);
    %     flucts(ptcleid) = sqrt(mean((spt(:,1)-mean(spt(:,1))).^2 + (spt(:,2)-mean(spt(:,2))).^2))*(nm_per_pixel);
    if ~isempty(spt)
        ptcle_index = ptcle_index + 1;
        for t = 0:round(0.5*length(spt))
            
            dx_temp = spt((1+t):end,1)-spt(1:(end-t),1);
            dy_temp = spt((1+t):end,2)-spt(1:(end-t),2);
            
            dr(t+1,ptcle_index) = mean(dx_temp.^2+dy_temp.^2);
        end
    end
    
    %     crr(1:size(dr,1),ptcleid) = 1 - dr(:,ptcleid)/2/(flucts(ptcleid)^2);
    
end

dr(1,:)=[];

dr(dr(:)==0)=NaN;
dr = [zeros(1,size(dr,2)); dr];

counter = dr;
counter(~isnan(dr(:)))=1;
counter(isnan(dr(:)))=0;

msd_total = nansum(dr,2)./nansum(counter,2)*(nm_per_pixel^2);
msd_std = sqrt(nansum((dr-repmat(msd_total,[1 size(dr,2)])).^2,2)./nansum(counter,2));

res = [msd_total msd_std dr*(nm_per_pixel^2)];

end