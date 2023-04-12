%% filepaths
csvFolder = "/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/testOut";
t = loadPreSegmentedNucGeomIntCSVs(csvFolder);

%% load table w annotations
% includes manual annotations that say which nuclei are correctly vs incorrectly segmented
% also computes extra variables, e.g. surfToVolRatio, DistRank (rank of the
% distance from center), CV (coeff of variation of the DAPI intensity)
pathToAnnotation = '/Users/lionnt01/Documents/data/feiyue/egg chamber image quant/testOut/MatlabTableFirstPassWAnnotation.txt';
t2 = addAnnotationsAndVarsToNucleiTable(t,pathToAnnotation);
%%
v = t2.Kurtosis;
figure; hold; plot(v(t2.isWrong==0),'o'); plot(v(t2.isWrong==1),'o');

%% cleanup nuclei based on multiple ad hoc criteria
maxVolume = 5000;
%minVolume = 200; %100 for less stringent
minVolume = 100; %this is the best, 200 for more stringent
%distRank weak separator by itself
maxSurfToVolRatio = 1; %this is the best, 0.8 for more stringent
%maxSurfToVolRatio = 0.8; %1 for less stringent
minSphericity = 0.4;
%maxLabel = 15; %17 for less stringent
maxLabel = 17; %this is the best, 15 for more stringent
maxLabel = 50; %essentially not a criterion, use for higher recall, lower precision
minMeanBreadth = 0;
maxMeanBreadth = 100;
maxSkewness = 1;
maxKurtosis = 3;
maxCV = 0.45;
minCV = 0.15;

idxClean = t2.Volume < maxVolume & t2.Volume > minVolume ...
    & t2.SurfToVolRatio < maxSurfToVolRatio ...
    & t2.Sphericity > minSphericity ...
    & t2.Label < maxLabel ...
    & t2.MeanBreadth > minMeanBreadth & t2.MeanBreadth < maxMeanBreadth ...
    & t2.Skewness < maxSkewness ...
    & t2.Kurtosis < maxKurtosis ...
    & t2.CV < maxCV & t2.CV > minCV;

t2clean = t2(idxClean,:);

n0 = sum(t2.isWrong == 0);
n1 = sum(t2.isWrong == 1);
disp(['initially ',num2str(n0),' good nuclei and ',num2str(n1),' bad ones.']);

n0 = sum(t2clean.isWrong == 0);
n1 = sum(t2clean.isWrong == 1);
disp(['After threshold based on individual metrics, ',...
    num2str(n0),' good nuclei and ',num2str(n1),' bad ones']);
disp(['       i.e. ',...
    num2str(n0/sum(t2.isWrong == 0)),' recall; ',...
    num2str(n0/(n0+n1)),' precision']);

%%  