function [ mvals,svals,xfit,yfit ] = fit_2_gaussian_bootstrap( ...
    r, itnum, histbin, use_weights, lb, ub )

% fits a distribution of datapoints to a sum of two gaussian populations:
% a1 / (sqrt(2*pi).*sig1) * exp( - (x - mu1)^2 / (2*sig1^2) ) ...
% + (1- a1) / (sqrt(2*pi).*sig2) * exp( - (x - mu2)^2 / (2*sig2^2) )

% INPUT
% r: the data (1D array of random variables from bimodal distribution)
% itnum: Number of iterations of the bootstrap algorithm (100 to 1000
% should do it).
% histbin:  bin size of the histogram used to fit the data (depends on your
% data)
% use_weights: set this flag to 1 to include statistical weights into 
    % the fit (sqrt(n_i) for each bar n_i); recommended. Set to 0 to ignore
    % statistical weights.

% lb, ub: fit parameters constraints
% lb = [a1min, mu1min, sig1min, mu2min, sig2min];
% ub = [a1max, mu1max, sig1max, mu2max, sig2max];
    % where:
    % a1 is the fraction (0 to 1) of the population in the low mean-value mode
    % mu1 is the mean value of the low mean-value mode
    % sig1 is the sigma of the low mean-value mode
    % mu2 is the mean value of the high mean-value mode
    % sig2 is the sigma of the high mean-value mode

% OUTPUT
% mvals: the median parameters of the 2-normal fit, in the following order:
% a1, the fraction (0 to 1) of the population in the low mean-value mode
% mu1, the mean value of the low mean-value mode
% sig1, the sigma of the low mean-value mode
% mu2, the mean value of the high mean-value mode
% sig2, the sigma of the high mean-value mode

% svals: the standard deviation of the parameters (i.e. error on the fit),
% same order as above

% xfit: is the x coordinate of the best fitted curve over the data range
% yfit: is the y coordinate of the best fitted curve over the data range

%% get started
n = numel(r);

% set fit convergence criteria
options = optimset('TolX',.001,'MaxIter',1000);

%% create fitting function
fit2gauss = fittype('a1./(sqrt(2*pi).*sig1).*exp(-(x-mu1).^2./(2*sig1.^2)) + (1-a1)./(sqrt(2*pi).*sig2).*exp(-(x-mu2).^2./(2*sig2.^2))',...
            'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'a1','mu1','sig1','mu2','sig2'});

fitopts = fitoptions('Method', 'LinearLeastSquares');

%% bootstrap loop

% generate the edges of each bin of the histogram
hist_edges = min(r):histbin:max(r); 

vals = [];
for i = 1:itnum
    % random sample of opriginal dataset
    data_tmp = r(randsample(n,n,1));
    
    % generate histogram
    [ntmp,xtmp] = histcounts(data_tmp,hist_edges);
    
    % normalize to 1
    ntmp = ntmp/(sum(ntmp)*histbin);
    
    % compute bin centers coordinates
    xtmp = (xtmp(1:end-1) + xtmp(2:end))/2;
    
    % generate guess values of the coefficients to initialize the fit
    m1tmp = quantile(data_tmp,0.25);
    m2tmp = quantile(data_tmp,0.75);
    s1tmp = (m2tmp - m1tmp)/2;
    s2tmp = (m2tmp - m1tmp)/2;
    
    % set options
    if use_weights == 1
        fitopts = fitoptions('Method', 'NonLinearLeastSquares',...
            'Weights',sqrt(ntmp'),...
            'Lower',lb,...
            'Upper',ub,...
            'Startpoint',[0.5,m1tmp,s1tmp,m2tmp,s2tmp]);
    else
        fitopts = fitoptions('Method', 'NonLinearLeastSquares',...
            'Lower',lb,...
            'Upper',ub,...
            'Startpoint',[0.5,m1tmp,s1tmp,m2tmp,s2tmp]);
    end
    
    % fit to 2 gaussians
    tmpfit = fit(xtmp',ntmp',fit2gauss,fitopts);
    
    % collect coefficients in vector
    tmpcoefs = coeffvalues(tmpfit);
    
    % make sure that the gaussian modes are ordered by increasing mean
    is_ordered = tmpcoefs(2) < tmpcoefs(4);
    if ~is_ordered
        tmpcoefs = [1-tmpcoefs(1) , ...
            tmpcoefs(4),...
            tmpcoefs(5),...
            tmpcoefs(2),...
            tmpcoefs(3)];         
    end
    
    % append ordered coefficients to coef list
    vals = [vals; tmpcoefs];
end

%% compute results statistics and fitted curve 
mvals = median(vals);
svals = std(vals);

% compute fitted curve over 1000 bins spanning the data range
xfit = min(hist_edges): ...
    (max(hist_edges) - min(hist_edges))/1000 : ...
    max(hist_edges);

% compute fit values
yfit = mvals(1)/(sqrt(2*pi).*mvals(3))*exp(-(xfit-mvals(2)).^2./(2*mvals(3).^2)) + ...
    (1-mvals(1))/(sqrt(2*pi).*mvals(5))*exp(-(xfit-mvals(4)).^2./(2*mvals(5).^2));

%normalize to cumulated sum = 1
yfit = yfit * histbin;

end

