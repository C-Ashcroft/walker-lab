function PlotSingleSubject

cd SingleTrial/ %change to the directory where your data is
 
e = dir; %make an array called e which contains one entry for each item in the directory
ind=1;
files=[];

for ii=1:length(e) %make a directory 'files' with all data files in
    
    if e(ii).isdir==0 %if this is a file and not a subdirectory
        files(ind).name=e(ii).name; %then enter its name and date in to 'files'
        files(ind).date=e(ii).date;
        ind=ind+1;
    end
end

qremove = ismember({files.name},{'.','..','.DS_Store','.MATLABDriveTag'}); % locate where the unwanted names are
files(qremove) = []; % set these names to 0
rawdata=[];
ind=1;
nind=0;
DatTab =[];

for ii=1:length(files) %for each file in our array
    rawdata(ind).file= files(ii).name; %enter these into 'rawdata'
    load(rawdata(ind).file) %load the file
    rawdata(ind).rawdata = d;
    DatTab =[DatTab; d];
    ind=ind+1;

end

%structure of data:
%  'Correct Pitch (1=high,0=low)' 'Was subject Correct (1,0)' 'Correction Trial(0=CT)' 'Target F0' 'Ref F0' 'Stage' 'Dur Ref' 'Dur Tar' 'Onset Delay' 'Level Tar' 'ISI' 'Ramp dur s' 'Sampling rate' 'delay' 'Timeout' 'Human''s response (1=high,0=low)' 'Reaction time'

stage = [];
stage = unique(DatTab(:,6));

for ss = 1:length(stage)

dur = unique(DatTab(:,7));

% identify rows of data for duration 1 (5ms) to duration 8 (132ms)
for dd = 1:length(dur)  
    finddurdata = [];
    finddurdata = find(DatTab(:,7) == dur(dd) & DatTab(:,6) == stage(ss));
    
    targets = unique(DatTab(:,4));
    percentcorrect = [];
    frequency = [];
   
    for tt = 1:length(targets) % get percent correct on each target pitch
        findtar = [];
        findtar = find(DatTab(:,4) == targets(tt));
        findtrials = intersect(findtar,finddurdata);
        frequency = [frequency; targets(tt)];
        percentcorrect = [percentcorrect ; mean(DatTab(findtrials,16))*100]; 
    
    end
    
%     create tables of target frequencies and percentage correct response for each
%     stimulus duration 
       if dd == 1
		stimdur8 = table(frequency, percentcorrect);
	elseif dd == 2
		stimdur20 = table(frequency, percentcorrect);
	elseif dd == 3
		stimdur36 = table(frequency, percentcorrect);
        elseif dd == 4 
                stimdur68 = table(frequency, percentcorrect);
                
       end      
       
end

% plot target frequency against % correct for each stim duration

if stage(ss) == 1

    figure(1); clf
    plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
    hold on
    plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
    plot(stimdur36.frequency, stimdur36.percentcorrect,'yx-', 'Linewidth', 1.5);
    plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);
    
    elseif stage(ss) == 2
    
    figure(3); clf
    plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
    hold on
    plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
    plot(stimdur36.frequency, stimdur36.percentcorrect,'yx-', 'Linewidth', 1.5);
    plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);
    
    
    elseif stage(ss) == 3
    
    figure(3); clf
    plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
    hold on
    plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
    plot(stimdur36.frequency, stimdur36.percentcorrect,'yx-', 'Linewidth', 1.5);
    plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);
    
    elseif stage(ss) == 4
    
    figure(4); clf
    plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
    hold on
    plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
    plot(stimdur36.frequency, stimdur36.percentcorrect,'yx-', 'Linewidth', 1.5);
    plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);
    
    
    elseif stage(ss) == 5
    
    figure(5); clf
    plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
    hold on
    plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
    plot(stimdur36.frequency, stimdur36.percentcorrect,'yx-', 'Linewidth', 1.5);
    plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);

    elseif stage(ss) == 6
    
    figure(6); clf
    plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
    hold on
    plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
    plot(stimdur36.frequency, stimdur36.percentcorrect,'yx-', 'Linewidth', 1.5);
    plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);


    
    elseif stage(ss) == 7
    
    figure(7); clf
    plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
    hold on
    plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
    plot(stimdur36.frequency, stimdur36.percentcorrect,'yx-', 'Linewidth', 1.5);
    plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);
    
    
    elseif stage(ss) == 8
    
    figure(8); clf
    plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
    hold on
    plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
    plot(stimdur36.frequency, stimdur36.percentcorrect,'yx-', 'Linewidth', 1.5);
    plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);
    
elseif stage(ss) == 9
    
    figure(9); clf
    plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
    hold on
    plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
    plot(stimdur36.frequency, stimdur36.percentcorrect,'yx-', 'Linewidth', 1.5);
    plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);
         
end
    
xlabel('Target Frequency (Hz)')
ylabel('Percentage "High" Response')
xlim([245 255])
ylim([0 100])
yline(50,'k--');
title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject')
legend('8ms','20ms','36ms','68ms','Location', 'best')
hold off


end
end 


 



