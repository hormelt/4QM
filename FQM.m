function res = FQM(subdata,fake_dx,fake_dy,calibrate,p_coef,plotopt)

% FQM returns the 3-dimensional subpixel locations of particles in subdata
% using the Four Quadrant Method.
%
% INPUTS:
%   subdata: image containing the particle to be tracked
%   fake_dx:
%   fake_dy:
%   calibrate: TRUE to have this function perform calibration
%   p_coef:
%
% OUTPUTS:
%   res: linear corrections in position in the form p1(1)*x+p1(2)...
%       p2(1)*y+p2(2) are output as res = [p1 errx p2 erry]
%


%% step through each particle and frame then calibrate or measure

cutoff = round(size(subdata,1)/2); %are we requiring even, odd input?

% Define quadrants

QLR = subdata(cutoff+1:end,cutoff+1:end,:); QUR = subdata(1:cutoff,cutoff+1:end,:);
QLL = subdata(cutoff+1:end,1:cutoff,:); QUL = subdata(1:cutoff,1:cutoff,:);

% Perform sums to obtain bias

A = squeeze(sum(sum(QUL))); B = squeeze(sum(sum(QUR)));
C = squeeze(sum(sum(QLL))); D = squeeze(sum(sum(QLR)));

if calibrate
    
    cnt = [(A+C-B-D)./(A+B+C+D) (A+B-C-D)./(A+B+C+D)];
    refShift = [fake_dx fake_dy];
    
    % now find the shift that gives the smallest error
    
    [p1,fvalx] = fminsearch(@(p1) squeeze(mean((p1(1)*(cnt(:,1)+p1(2))-refShift(:,1)).^2,1)),...
        [range(refShift(:,1))/range(cnt(:,1)),mean(refShift(:,1))]);
    errx = sqrt(fvalx);
    
    [p2,fvaly] = fminsearch(@(p2) squeeze(mean((p2(1)*(cnt(:,2)+p2(2))-refShift(:,2)).^2,1)),...
        [range(refShift(:,2))/range(cnt(:,2)),mean(refShift(:,2))]);
    erry = sqrt(fvaly);
    
    if (errx<=1e-1) && (erry<=1e-1)
        
        %csvwrite([num2str(round(rand*1000)) '.csv'],[p1(1)*(cnt(:,1)+p1(2)),refshft(:,1) p2(1)*(cnt(:,2)+p2(2)),refshft(:,2)]);
        
        res = [p1 errx p2 erry];
        
        if plotopt
            
            scatter(p1(1)*(cnt(:,1)+p1(2)),refShift(:,1),'b')
            hold on
            scatter(p2(1)*(cnt(:,2)+p2(2)),refShift(:,2),'g')
            getframe;
            
        end
        
    else
        
        res = [NaN NaN NaN NaN NaN NaN];
        
    end
    
else
    
    res = [p_coef(1)*(A+C-B-D)./(A+B+C+D)+p_coef(2) p_coef(4)*(A+B-C-D)./(A+B+C+D)+p_coef(5) [1:size(subdata,3)]'];
    
end

end