function s = modzscore(x)
%modified z score using median
s = (x - median(x))./std(x);
end