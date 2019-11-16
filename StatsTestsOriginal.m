clear all 
clear all hidden 

names = {'will';'seb';'noah';'leonie';'gemma';'dom'; 'caitlin';'adam'};

allresults=[];
allfits=[];
allslopes=[];
allsubjects=[];
figure(1);clf
hold on

for subjects=1:length(names)
    
        thissubject=names{subjects};
        %find all data files for this subject
        files=dir([pwd '/SlopeDataOriginal/' thissubject '_results.mat' ]);
 
        % Load the file you want to analyse
        load([pwd '/SlopeDataOriginal/' thissubject '_results.mat']);
   
        
        allslopes=[allslopes;slopes];
        allresults=[allresults;results];
        allfits=[allfits;resultsfit];
        allsubjects=[allsubjects; ones(size(slopes,1),1).*subjects];
   
    dur = unique(allslopes(:,6)); 
    stage = unique(allslopes(:,5));

for ss = 1:length(stage)
    
    findslopes = [];
    findslopes = find(slopes(:,5) == stage(ss));
    slopedata = slopes(findslopes, 3);
    
    if ss == 1 
        
        slopestage1 = slopedata;
        
    elseif ss == 2
        
        slopestage2 = slopedata;
        
    elseif ss == 3
        
        slopestage3 = slopedata;
   
    elseif ss == 4
        
    slopestage4 = slopedata;
    
    end 
    
end

    figure(1);
    subplot(2,5,subjects);
    plot(dur, slopestage1,'rx-')
    hold on
    plot(dur, slopestage2,'bx-')
    plot(dur, slopestage3,'kx-')
    plot(dur, slopestage4,'gx-')
    title(['Subject ', num2str(subjects)])
    xlabel('Duration (ms)')
    ylabel('Slope Value')
    %legend('All Harm', 'Low Harm', 'High Harm', 'All Harm Random Phase')
    ylim([0 600])
    xticks([8 20 36 68])
    hold off 
     
          
end
 

[p,table,stats] = anovan(allslopes(:,[3]),{allslopes(:,5),allslopes(:,6),allsubjects},'model','interaction','varnames',{'stage','duration','subject'});
[p,table,stats] = anovan(allslopes(:,[3]),{allslopes(:,5),allslopes(:,6)},'model','interaction','varnames',{'stage','duration'});
figure()
multcompare(stats,'Dimension',[1])
figure()
multcompare(stats,'Dimension',[2])

for dd = 1:length(dur)

    finddata = [];
    finddata = find(allslopes(:,6) == dur(dd));
    [p,table,stats] = anovan(allslopes(finddata,[3]),allslopes(finddata,5),'model','interaction','varnames','stage');
    figure()
    multcompare(stats,'Dimension',[1]);
    hold on
    title(['Across Stage Comparison: ', num2str(dur(dd)), 'ms'])
    hold off
    
end 

% % Run anova at only one signal:noise level
% ff=find(allslopes(:,7)==40); %
% [p,table,stats] = anovan(allslopes(ff,[3]),{allslopes(ff,6)});
% figure()
% multcompare(stats,'Dimension',[1])
