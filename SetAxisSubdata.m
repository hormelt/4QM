function Axis = SetAxisSubdata(CoarseCoord,FeatSize,DeltaFit)

% Sets the coordinates in pixels of one of the axes of a subdata frame
% around a particle of interest.
%
% INPUT:
%   CoarseCoord: Coordinate of the particle of interest.
%   FeatSize: Full optical diameter of particle (pixels).
%   DeltaFit: Narrows analysis region around particle (pixels).
%
% OUTPUT:
%   Axis: The pixel numbers of the pixels around the particle. 

Offset = ceil(FeatSize/2)+DeltaFit;      
if round(CoarseCoord) > CoarseCoord
    Axis = (round(CoarseCoord)-(Offset)): ...
           (round(CoarseCoord)+(Offset))-1;
else
    Axis = (round(CoarseCoord)-(Offset))+1: ...
           (round(CoarseCoord)+(Offset));
end