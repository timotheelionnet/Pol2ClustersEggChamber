function [ptsl,err,ptsf,dist,elong] = linearFit3D(pts)
%pts should be a 3 x n array where each col are the x y z coordinates of a point


pts0 = mean(pts,2);
A = pts-pts0;
[U,~,~] = svd(A);
d = U(:,1);
t = d'*A;
t1 = min(t);
t2 = max(t);
ptsl = pts0 + [t1,t2].*d; % size 3x2

% ptsl is 3x2 contains coordinates of 2 points of the lines
% ptsl(:,1) is the 3D coordinates of the first point
% ptsl(:,2) is the 3D coordinates of the second point
% The line equation (parametric form) is then
% xyz(t) = t*ptsl(:,1) + (1-t)*(ptsl(:,2);

% compute parameter value of the point on the fitted line that 
% is closest from each original point
tf = computeParameterOfLinProjectionAlongAxis(pts,ptsl);

% compute elongation, i.e. max spread along the line of points closest to
% the data
elong = (max(tf) - min(tf));

% compute corresponding coordinate of the nearest point
ptsf = ptsl(:,2) + tf.*(ptsl(:,1) - ptsl(:,2))./norm(ptsl(:,1) - ptsl(:,2));

dist = pts - ptsf;

% err is the sqrt of the mean squared distance from each data point to the fitted line
err = mean(    (ptsf(1,:) - pts(1,:)).^2 ...
            + (ptsf(2,:) - pts(2,:)).^2 ...
            + (ptsf(3,:) - pts(3,:)).^2  );
err = sqrt(err);

