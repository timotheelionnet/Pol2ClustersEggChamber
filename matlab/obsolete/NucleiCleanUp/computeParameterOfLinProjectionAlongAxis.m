function t = computeParameterOfLinProjectionAlongAxis(pts,ptsl)

% input:
% pts is a 3 x n array of [x; y; z] coordinates for n data points.
% ptsl is a 3 x 2 array of [x; y; z] coordinates for 2 points defining an
% axis.

% output t is a 1xn vector which measures the euclidian distance along the
% axis of the projection of each original datapoint. ptsl(:,2) is used as
% the origin of the axis.



% compute parameter value of the point on the fitted line that 
% is closest from each original point
t =  sum( ...
    (pts - repmat(ptsl(:,2),1,size(pts,2))) ...
    .* repmat( ptsl(:,1) - ptsl(:,2)  ,1,size(pts,2)) ...
    ) ...
    ./repmat(norm(ptsl(:,1) -ptsl(:,2)),1,size(pts,2));