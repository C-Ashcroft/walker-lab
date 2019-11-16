function ThresholdOriginal

names = {'CA1';'CA2';'CA3';'CA4'; 'CA5';'CA6';'CA7'; 'CA8'};


allsubjectspercentcorrect = [];

for subjects = 1:length(names) 
    
    thissubject=names{subjects};
    
    cd([pwd '/TrialOriginal/' thissubject '/']) %change to the directory where your data is
 
    e = dir; %make an array called e which contains one entry for each item in the directory
    ind=1;
    files=[];
    cd ..
    cd ..
    
    for ii=1:length(e) %make a directory 'files' with all data files in

        if e(ii).isdir==0 %if this is a file and not a subdirectory
            if e(ii).name(1)~='.'
                files(ind).name=e(ii).name; %then enter its name and date in to 'files'
                files(ind).date=e(ii).date;
                ind=ind+1;
            end
        end
    end

    rawdata=[];
    ind=1;
    nind=0;
    DatTab =[];
    numbercorrect = [];
    percentcorrect = [];
  
    
    for ii=1:length(files) %for each file in our array
        rawdata(ind).file= files(ii).name; %enter these into 'rawdata'
        load([pwd '/TrialOriginal/' thissubject '/' rawdata(ind).file]) %load the file
        rawdata(ind).rawdata = d;
        DatTab =[DatTab; d];
        ind=ind+1;

    end
    
%structure of data:
%  'Correct Pitch (1=high,0=low)' 'Was subject Correct (1,0)' 'Correction Trial(0=CT)' 'Target F0' 'Ref F0' 'Stage' 'Dur Ref' 'Dur Tar' 'Onset Delay' 'Level Tar' 'ISI' 'Ramp dur s' 'Sampling rate' 'delay' 'Timeout' 'Human''s response (1=high,0=low)' 'Reaction time'

    correcttrials = find(DatTab(:,2) == 1);
    numbercorrect = length(correcttrials);
    percentcorrect = (numbercorrect/length(DatTab))*100;
    allsubjectspercentcorrect = [allsubjectspercentcorrect;percentcorrect];
  
end

figure (1);clf
histogram(allsubjectspercentcorrect,'BinWidth', 3.3)
hold on 
title('Overall Subject Performance Across Trials: Experiment 1')
xlabel('Total Percent Correct (%)')
ylabel('Number of Participants') 
xlim([50 100]);
yticks([0 1 2 3])
xline(65, '--r', 'Linewidth', 1.5);
hold off

end 

