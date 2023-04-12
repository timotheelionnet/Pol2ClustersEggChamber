function plotNucleiStatsByGroup(t,varName)
    addpath('cbrewer');
    conds = unique(t.condIdx);
    
    hSpacing = 1;
    
    figure('Name',[varName,' vs. condition (grouped nuclei)']);
    hold;
    
    % generate unique color map per nucleus
    ctr = 0;
    for i=1:numel(conds)
        samples = unique(t.sampleIdx(t.condIdx == conds(i)));
        ctr = ctr + numel(samples);
    end
    cm = cbrewer('qual','Dark2',ctr);
    
    ctr = 1;
    for i=1:numel(conds)
        samples = unique(t.sampleIdx(t.condIdx == conds(i)));
        
        for j=1:numel(samples)
            curT = t(t.sampleIdx == samples(j) & t.condIdx == conds(i),:);
            x = [];
            y = [];
            for k=0:4
                curY = curT.(varName)(curT.ecGroup == k);
                curX = (k*(numel(conds)+hSpacing)+i*hSpacing ...
                    + j/numel(samples)*hSpacing*0.4 )*ones(size(curY));
                y = [y; curY];
                x = [x; curX];  
            end
            plot(x,y,'o','MarkerEdgeColor',cm(ctr,:),...
                    'MarkerFaceColor',cm(ctr,:),...
                    'DisplayName',['Cond. ',num2str(conds(i)),...
                    '; sample ',num2str(samples(j))]);
            ctr = ctr+1;
        end
    end