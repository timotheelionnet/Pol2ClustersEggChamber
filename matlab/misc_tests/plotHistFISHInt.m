% this code takes values of FISH intensity per cell in different conditions
% (input as different variables, see next section for data format)
% fits them to a sum of lognormal distributions, plots and saves the result of the
% fit.

%% Set up here the variables, their names for plot labels, their color
% this line for demo only, comment out if using your own data
load('histFISHtestDataSet.mat');

% make sure that there are the same number of elements in each of the three
% following variables!

% variables start here:

% each entry of this cell array should be a vertical vector corresponding 
% to one of the conditions of the experiment.
condVars = {mCherry,ints11};

% each entry of this cell array should be a character array that stores the
% name of each condition to be used in plot labels.
condNames = {'mCherry RNAi','Ints11 RNAi'}; 

% each entry of this cell array stores a 3x1 array of RGB values (0-1) used
% to plot the corresponding conditions.
condColors = {  [0.5,0.5,0.5];... % gray
                [0.5,0,0]};       % red

% this is the value of the background to be subtracted from each intensity
% value.
backgroundInt = 100;

% bin size of the histogram of the FISH intensity in the final plot
histBin = 150; 

% this stores the output table of the fit results. Without a directory
% name, it will be saved in the active matlab folder (likely where your
% code is).
% saves two files, one called <outFileName>Res.txt and one called
% <outFileName>Params.txt.
outFileName = 'HistFISH';

% Note that the fit boundaries might need to be adjusted to match different
% experiments. Go to the code lines where these values are set and adjust the values:
    % histLogBin = 0.1; % bin size for the histogram of the log transformed variable
    % a1min = 0.1; % minimum fraction of low-intensity population (0-1)
    % a1max = 0.9; % maximum fraction of low-intensity population (0-1)
    % mu1min = 1; % minimum value of the mean of the low-intensity population in the log-transformed variable
    % mu1max = 7; % maximum value of the mean of the low-intensity population in the log-transformed variable
    % sig1min = 0.1; % minimum value of the std of the low-intensity population in the log-transformed variable
    % sig1max = 1.2; % maximum value of the std of the low-intensity population in the log-transformed variable
    % mu2min = 7; % same as mu1min for high-intensity population
    % mu2max = 12; % same as mu1max for high-intensity population
    % sig2min = 0.1; % same as sig1min for high-intensity population
    % sig2max = 1.2; % same as sig1max for high-intensity population

%% subtract background from intensity values
condVarsCorr = cell(size(condVars));

for i=1:numel(condVars)
    condVarsCorr{i} = condVars{i} - backgroundInt;
end    

%%  Side by side scatter plot of ctrl vs RNAi
figure('Name','Scatter plots intensities');
hold;
for i=1:numel(condVarsCorr)
    scatter(i-1+ (1:numel(condVarsCorr{i}))/(2*numel(condVarsCorr{i})),condVarsCorr{i},'o',...
        'MarkerFaceColor',condColors{i},'MarkerEdgeColor',condColors{i});
end

alpha(0.25);
xlim([-0.5,2]);
xticks([0.25,1.25]);
xticklabels(condNames);
xtickangle(45);
ylabel('H3 RNA FISH intensity per cell');

%% fit the histogram of log(intensity+1) to the sum of two gaussians

% set up the parameters for fit of the log-transformed variable to a sum of 2 gaussians (bootstrap)
itnum = 100; % number of iterations
use_weights = 1; % whether to use statistical weights in the fit

%----------------------------------------------------------------------------------
% the following parameters might need to be adjusted to match different experiments
histLogBin = 0.1; % bin size for the histogram of the log transformed variable
a1min = 0.1; % minimum fraction of low-intensity population (0-1)
a1max = 0.9; % maximum fraction of low-intensity population (0-1)
mu1min = 1; % minimum value of the mean of the low-intensity population in the log-transformed variable
mu1max = 7; % maximum value of the mean of the low-intensity population in the log-transformed variable
sig1min = 0.1; % minimum value of the std of the low-intensity population in the log-transformed variable
sig1max = 1.2; % maximum value of the std of the low-intensity population in the log-transformed variable
mu2min = 7; % same as mu1min for high-intensity population
mu2max = 12; % same as mu1max for high-intensity population
sig2min = 0.1; % same as sig1min for high-intensity population
sig2max = 1.2; % same as sig1max for high-intensity population
% end of adjustable parameters.
%------------------------------

% combine all boundary conditions into two vectors
lb = [a1min, mu1min, sig1min, mu2min, sig2min];
ub = [a1max, mu1max, sig1max, mu2max, sig2max];

% compute histograms of the log-transformed variables, and fit it to
% the sum of two gaussians
for i=1:numel(condVarsCorr)

    [nLog{i},xLog{i}] = hist( log(condVarsCorr{i}+1),...
        0:histLogBin:max(log(condVarsCorr{i}+1))+histLogBin);

    [ mVals{i},sVals{i},xFitLog{i},yFitLog{i} ] = fit_2_gaussian_bootstrap( ...
    log(condVarsCorr{i}+1), itnum, histLogBin, use_weights, lb, ub );
end

% plot log-transformed results of the fit
figure('Name','Log-transformed Histograms and fits');
hold;
for i=1:numel(condVarsCorr)
    plot(xLog{i}-1,nLog{i},'Color',condColors{i},...
        'DisplayName',[condNames{i},' log(x+1) histogram']);
    plot(xFitLog{i}-1,yFitLog{i}*sum(nLog{i}),'Color',condColors{i},...
        'LineWidth',2,'DisplayName',[condNames{i},' fit']);
end
xlabel('log(x+1) transform of H3 RNA FISH intensity per cell');
ylabel('Counts');


%% reverse-transform results of the log-transform fit back into linear space

% factor used to define the sampling interval 
    % of the fitted function for a smooth plot. 
    % E.g. if the histogram of the data is computed with a bin size histBin = 150, 
    % the fitted function will be evaluated every histBin/overSamplingFactor, 
    % i.e. 15 if overSamplingFactor = 10.
overSamplingFactor = 10; 

% complute global maximum value    
globMax = -Inf;
for i=1:numel(condVarsCorr)
    if max(condVarsCorr{i}) > globMax
        globMax = max(condVarsCorr{i});
    end
end

x = cell(numel(condVarsCorr),1);
n = cell(numel(condVarsCorr),1);
xfit = cell(numel(condVarsCorr),1);
yfit = cell(numel(condVarsCorr),1);
for i=1:numel(condVarsCorr)
    % compute the histogram of the data in linear space
    [n{i},x{i}] = hist(condVarsCorr{i},0:histBin:globMax+histBin);
    
    % compute the values of the fit function
    xfit{i} = 0:histBin/overSamplingFactor:(globMax+histBin);

    yfit{i} = (mVals{i}(1)*normpdf(log(xfit{i}+1),mVals{i}(2),mVals{i}(3)) ...
        + (1-mVals{i}(1))*normpdf(log(xfit{i}+1),mVals{i}(4),mVals{i}(5)))./(xfit{i}+1);

    yfit{i} = overSamplingFactor*yfit{i}/sum(yfit{i})*sum(n{i});
end

fh = figure('Name','Histogram + fit');
vSpace = 0.1;
vHeight = (1-vSpace*(numel(condVarsCorr)+1))/numel(condVarsCorr);
ah = cell(numel(condVarsCorr),1);
for i=1:numel(condVarsCorr)
    ah{i} = axes('Parent',fh,'Position',[0.1,1-i*(vHeight+vSpace),0.8,vHeight],...
        'Units','normalized');
    hold;
    xlabel('FISH intensity');
    ylabel('Counts');
end

for i=1:numel(condVarsCorr)
    % plot data as bar plot
    h = bar(ah{i},x{i},n{i},'FaceColor',condColors{i},...
        'EdgeColor',condColors{i},...
        'FaceAlpha',0.5,...
        'EdgeAlpha',0,...
        'LineWidth',1,'BarWidth',1);
    set(h,'DisplayName',[condNames{i},' data']);
    
    % plot fit;
    % note that the formula of the variance of the log-normal distribution if the
    % log transformed distribution has mean mu and std sigma is:
    % (exp(sigma^2) - 1)*exp(2*mu + sigma^2)
    % from https://en.wikipedia.org/wiki/Log-normal_distribution
    plot(ah{i},xfit{i},yfit{i},'Color',condColors{i},'LineWidth',2,...
        'DisplayName',[condNames{i},' fit:',newline,...
        'a_1 = ',num2str(mVals{i}(1),'%.2f'),...
        ' +/- ',num2str(sVals{i}(1),'%.1f'),...
        '; \mu_1 = ',num2str(exp(mVals{i}(2))-1,'%.0f'),...
        ' +/- ',num2str((exp(sVals{i}(2)+mVals{i}(2))-exp(-sVals{i}(2)+mVals{i}(2)))/2,'%.0f'),...
        '; \sigma_1 = ',...
        num2str( logNormStd(mVals{i}(2),mVals{i}(3)),'%.0f'),...
        ' +/- ',num2str( (logNormStd(mVals{i}(2),mVals{i}(3)+sVals{i}(3)) ...
        - logNormStd(mVals{i}(2),mVals{i}(3)-sVals{i}(3)))/2,'%.0f'),...
        ';',newline,...
        'a_2 = ',num2str(1-mVals{i}(1),'%.2f'),...
        ' +/- ',num2str(sVals{i}(1),'%.1f'),...
        '; \mu_2 = ',num2str(exp(mVals{i}(4))-1,'%.0f'),...
        ' +/- ',num2str((exp(mVals{i}(4)+ sVals{i}(4))- exp(mVals{i}(4)- sVals{i}(4)))/2,'%.0f'),...
        '; \sigma_2 = ',...
        num2str( logNormStd(mVals{i}(4),mVals{i}(5)),'%.0f'),...
        ' +/- ',num2str( (logNormStd(mVals{i}(4),mVals{i}(5)+sVals{i}(5)) ...
        - logNormStd(mVals{i}(4),mVals{i}(5)-sVals{i}(5)))/2,'%.0f')]);

    axes(ah{i});
    legend show;
end

%% output the fitting results as a table
a1 = zeros(numel(condVarsCorr),1);
a2 = zeros(numel(condVarsCorr),1);
mu1 = zeros(numel(condVarsCorr),1);
mu2 = zeros(numel(condVarsCorr),1);
sig1 = zeros(numel(condVarsCorr),1);
sig2 = zeros(numel(condVarsCorr),1);

a1Err = zeros(numel(condVarsCorr),1);
a2Err = zeros(numel(condVarsCorr),1);
mu1Err = zeros(numel(condVarsCorr),1);
mu2Err = zeros(numel(condVarsCorr),1);
sig1Err = zeros(numel(condVarsCorr),1);
sig2Err = zeros(numel(condVarsCorr),1);

for i=1:numel(condVarsCorr)
    a1(i,1) = mVals{i}(1);
    a1Err(i,1) = sVals{i}(1);
    
    mu1(i,1) = exp(mVals{i}(2))-1;
    mu1Err(i,1) = (exp(mVals{i}(2)+sVals{i}(2))-exp(mVals{i}(2)-sVals{i}(2)))/2;
    sig1(i,1) = logNormStd(mVals{i}(2),mVals{i}(3));
    sig1Err(i,1) = (logNormStd(mVals{i}(2),mVals{i}(3)+sVals{i}(3)) ...
        - logNormStd(mVals{i}(2),mVals{i}(3)-sVals{i}(3)))/2;

    a2(i,1) = 1 - mVals{i}(1);
    a2Err(i,1) = sVals{i}(1);
    
    mu2(i,1) = exp(mVals{i}(4))-1;
    mu2Err(i,1) = (exp(mVals{i}(4)+sVals{i}(4))-exp(mVals{i}(4)-sVals{i}(4)))/2;
    sig2(i,1) = logNormStd(mVals{i}(4),mVals{i}(3));
    sig2Err(i,1) = (logNormStd(mVals{i}(4),mVals{i}(5)+sVals{i}(5)) ...
        - logNormStd(mVals{i}(4),mVals{i}(5)-sVals{i}(5)))/2;

end

t =table(condNames',a1,a1Err,mu1,mu1Err,sig1,sig1Err,a2,a2Err,mu2,mu2Err,sig2,sig2Err);
t = renamevars(t,{'Var1'},{'condition'});

writetable(t,[outFileName,'Res.txt'],'Delimiter','\t');

%% output the parameters
parameterName = {   'backgroundInt';...
                    'histBin'; ...
                    'itnum';...
                    'use_weights';...
                    'histLogBin';...
                    'a1min';...
                    'a1max';...
                    'mu1min';...
                    'mu1max';...
                    'sig1min';...
                    'sig1max';...
                    'mu2min';...
                    'mu2max';...
                    'sig2min';...
                    'sig2max'};

parameterValue = [  backgroundInt;...
                    histBin; ...
                    itnum;...
                    use_weights;...
                    histLogBin;...
                    a1min;...
                    a1max;...
                    mu1min;...
                    mu1max;...
                    sig1min;...
                    sig1max;...
                    mu2min;...
                    mu2max;...
                    sig2min;...
                    sig2max];

tp = table(parameterName,parameterValue);
writetable(tp,[outFileName,'Params.txt'],'Delimiter','\t');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENTER OBSOLETE CHUNK OF CODE - IGNORE , IT WILL NOT EXECUTE
% note there is a useful function after the obsolete chunk.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the following is obsolete for reference, ignore
runObsoleteCode = 0;

%% histogram, fit to two gaussians
% the gaussian fit of the distributions fails to capture the long tails -
% not optimal
if runObsoleteCode
    


    backgroundInt = 100;
    mCherryCorr = mCherry - backgroundInt;
    ints11Corr = ints11-backgroundInt;
    
    % fit to 2 gaussian (bootstrap)
    itnum = 1000;
    histBin = 100;
    use_weights = 0;
    
    [nCh,xCh] = hist(mCherryCorr,0:histBin:10000);
    [nInts,xInts] = hist(ints11Corr,0:histBin:10000);
    
    a1min = 0.1;
    a1max = 0.9;
    mu1min = 50;
    mu1max = 1200;
    sig1min = 50;
    sig1max = 500;
    mu2min = 1200;
    mu2max = 10000; 
    sig2min = 50;
    sig2max = 10000;
    lb = [a1min, mu1min, sig1min, mu2min, sig2min];
    ub = [a1max, mu1max, sig1max, mu2max, sig2max];
    [ mvalsCh,svalsCh,xfitCh,yfitCh ] = fit_2_gaussian_bootstrap( ...
        mCherryCorr, itnum, histBin, use_weights, lb, ub );
    
    [ mvalsInts,svalsInts,xfitInts,yfitInts ] = fit_2_gaussian_bootstrap( ...
        ints11Corr, itnum, histBin, use_weights, lb, ub );
    %
    figure('Name','Histogram Ints11 RNAi vs ctrl');
    hold;
    plot(xCh,nCh,'-o','Color',[0.5,0.5,0.5],...
        'MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor',[0.5,0.5,0.5]);
    plot(xInts,nInts,'-o','Color',[0.75,0,0],...
        'MarkerFaceColor',[0.75,0,0],'MarkerEdgeColor',[0.75,0,0]);
    
    plot(xfitCh,yfitCh*sum(nCh),'Color',[0.5,0.5,0.5],'LineWidth',2);
    plot(xfitInts,yfitInts*sum(nInts),'Color',[0.75,0,0],'LineWidth', 2);
    xlim([0,10000]);
end

%% fit to sum of two lognormals (doesn't work that well)
if runObsoleteCode

    backgroundInt = 100;
    mCherryCorr = mCherry - backgroundInt;
    ints11Corr = ints11-backgroundInt;

    % this works for ints11
    x = 0:10000;
    y = 16000*lognpdf(x,6.6,0.3)+18000*lognpdf(x,7.5,0.2);
    plot(x,y);
    % this works for ctrl
    x = 0:10000;
    y = 26000*lognpdf(x,5.9,0.9)+13000*lognpdf(x,8.2,0.27);
    plot(x,y);

    figure('Name','Fit to 2 Log-normal histograms Ints11 RNAi vs ctrl');
    hold;
    histBin = 100;
    use_weights = 0;
    
    itnum = 100;
    
    a1min = 0;
    a1max = 1;
    mu1min = 5.6;
    mu1max = 7;
    sig1min = 0.2;
    sig1max = 0.9;
    mu2min = 7.5;
    mu2max = 9;
    sig2min = 0.1;
    sig2max = 0.3;
    lb = [a1min, mu1min, sig1min, mu2min, sig2min];
    ub = [a1max, mu1max, sig1max, mu2max, sig2max];
    
    [nCh,xCh] = hist(mCherryCorr,0:histBin:10000);
    [nInts,xInts] = hist(ints11Corr,0:histBin:10000);
    
    [ mvalsCh,svalsCh,xfitCh,yfitCh ] = fit_2_lognormal_bootstrap( ...
        mCherryCorr, itnum, histBin, use_weights, lb, ub );
    
    [ mvalsInts,svalsInts,xfitInts,yfitInts ] = fit_2_lognormal_bootstrap( ...
        ints11Corr, itnum, histBin, use_weights, lb, ub );
    
    plot(xCh,nCh,'Color',[0.5,0.5,0.5]);
    plot(xInts,nInts,'Color',[0.75,0,0]);
    plot(xfitCh,yfitCh*sum(nCh),'Color',[0.5,0.5,0.5],'LineWidth',2);
    plot(xfitInts,yfitInts*sum(nInts),'Color',[0.75,0,0],'LineWidth', 2);
    xlim([0,10000]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%END OF OBSOLETE CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% the formula of the variance of the log-normal distribution if the
% log transformed distribution has mean mu and std sigma is:
% (exp(sigma^2) - 1)*exp(2*mu + sigma^2)
% from https://en.wikipedia.org/wiki/Log-normal_distribution
function s = logNormStd(mu,sigma)
    s = sqrt((exp(sigma^2) - 1)*exp(2*mu + sigma^2));
end

%% 
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

