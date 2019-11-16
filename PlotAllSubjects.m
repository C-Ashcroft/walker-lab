function PlotAllSubjects
%     Direct to data files  
cd  /Users/Caitlin/MATLAB-Drive/TrialOriginal
    
    myFolder          = dir;
    
    qremove           = ismember({myFolder.name},{'.','..','.DS_Store','.MATLABDriveTag'}); % locate where the unwanted names are
    myFolder(qremove) = []; % set these names to 0
    participant_number     = length(myFolder);
    
    %% Define necessary variables 
        
    
        allpercenthigh8 = [];
        allpercenthigh20 = [];
        allpercenthigh36 = [];
        allpercenthigh68 = [];
        allstimdur8  = [];
        allstimdur20 = [];
        allstimdur36 = [];
        allstimdur68 = [];
        stdstim8 = [];
        stdstim20 = [];
        stdstim36 = [];
        stdstim68 = [];
        DatTabTotal = [];
     
    
    %% Load data
    for k = 1:length(myFolder)
        
        folderName      = myFolder(k).name;
        
        fileName        = dir([folderName '/*.mat']);
        
        DatTab = [];
       
        
        for a = 1:length(fileName)
            
            
            data        = load([folderName filesep fileName(a).name]);
            
            onestage = data.d;
            
            DatTab = [DatTab;onestage];
            
            DatTabTotal = [DatTabTotal;DatTab];
                   

        end 

        %structure of data:
        %  'Correct Pitch (1=high,0=low)' 'Was subject Correct (1,0)' 'Correction Trial(0=CT)' 'Target F0' 'Ref F0' 'Stage' 'Dur Ref' 'Dur Tar' 'Onset Delay' 'Level Tar' 'ISI' 'Ramp dur s' 'Sampling rate' 'delay' 'Timeout' 'Human''s response (1=high,0=low)' 'Reaction time'
        
        stage = [];
        stage = unique(DatTab(:,6));
        
        for ss = 1:length(stage)
        
        dur = unique(DatTab(:,7));
        
        % identify rows of data for duration 1 (8ms) to duration 4 (68ms)
        for dd = 1:length(dur)  
            finddurdata = [];
            finddurdata = find(DatTab(:,7) == dur(dd) & DatTab(:,6) == stage(ss));
            
            targets = unique(DatTab(:,4));
            percenthigh = [];
            frequency = [];
           
            for tt = 1:length(targets) % get percent correct on each target pitch
                findtar = [];
                findtar = find(DatTab(:,4) == targets(tt));
                findtrials = intersect(findtar,finddurdata);
                frequency = [frequency; targets(tt)];
                percenthigh = [percenthigh ; mean(DatTab(findtrials,16))*100];
            
            end
            
        %     create tables of target frequencies and percentage correct response for each
        %     stimulus duration 
        
                if dd == 1
                        allstimdur8(k,:,ss) = percenthigh;
        	    elseif dd == 2
                        allstimdur20(k,:,ss) = percenthigh;
        	    elseif dd == 3
                        allstimdur36(k,:,ss) = percenthigh;
                elseif dd == 4 
                        allstimdur68(k,:,ss) = percenthigh;
                        
               end   
                    
               
        end
      
            
        end 
        
    end
    
%     calculate standard error across participants
        
        stdstim8 = std(allstimdur8)/sqrt(participant_number)
        stdstim20 = std(allstimdur20)/sqrt(participant_number)
        stdstim36 = std(allstimdur36)/sqrt(participant_number)
        stdstim68 = std(allstimdur68)/sqrt(participant_number)

% Begin calculating mean % correct for all participants at each stage 
        
stage = unique(DatTabTotal(:,6)); 
results=[];
resultsfit=[];
slopes=[];

for ss = 1:length(stage)

dur = unique(DatTabTotal(:,7));

% identify rows of data correspnding to duration 1 (8ms) to duration 4 (68ms)
for dd = 1:length(dur)  
    finddurdata = [];
    finddurdata = find(DatTabTotal(:,7) == dur(dd) & DatTabTotal(:,6) == stage(ss));
    
    targets = unique(DatTabTotal(:,4));
    percenthigh = [];
    frequency = [];
    ntrials = [];
    stagevar= [];
    durvar = [];
   
    for tt = 1:length(targets) % get percent correct on each target pitch
        findtar = [];
        findtar = find(DatTabTotal(:,4) == targets(tt));
        findtrials = intersect(findtar,finddurdata);
        frequency = [frequency; targets(tt)];
        percenthigh = [percenthigh ; mean(DatTabTotal(findtrials,16))*100];
        ntrials = [ntrials ; length(findtrials)];
        stagevar = [stagevar ; stage(ss)];
        durvar = [durvar ; dur(dd)];
    
    end
    
    %calculate curve fit
        [b,dev,stats] = glmfit(log(frequency),[round(percenthigh.*ntrials./100) ntrials],'binomial','logit'); %run a probit regression to produce a fit to the data
        y = glmval(b,log(frequency),'logit');
        xoct=[fliplr(octavesteps(DatTab(1,5),-12,30)) octavesteps(DatTab(1,5),12,30)];%get x values that are equally spaced in octave steps:
        xoct=unique(xoct);
        yoct = glmval(b,log(xoct),'logit');
        xfreq=[frequency(1):0.5:frequency(end)]';%get x values that are equally spaced in 0.5Hz steps
        yfreq = glmval(b,log(xfreq),'logit');
        midpoint=(log(0.5/(1-0.5))-b(1))/b(2); 
        bias=exp(midpoint)-DatTab(1,5); %bias
        OctSteps=octaves(xoct(7),xoct(8)); 
        mOct=gradient(yoct*100,OctSteps);
        Pslope=max(mOct); %calculate the gradient of the function
       % mFreq=gradient(yfreq*100,50);
       % PslopeF=max(mFreq);
    
    
%     create tables of target frequencies and percentage correct response for each
%     stimulus duration 
       if dd == 1
		stimdur8 = table(frequency, percenthigh);
        stimdur8_fit = table(xfreq,yfreq)
        stimdur8_slope(dd,stage(ss)) = Pslope;
       elseif dd == 2
		stimdur20 = table(frequency, percenthigh);
        stimdur20_fit = table(xfreq,yfreq)
        stimdur20_slope(dd,stage(ss)) = Pslope;
       elseif dd == 3
		stimdur36 = table(frequency, percenthigh);
        stimdur36_fit = table(xfreq,yfreq)
        stimdur36_slope(dd,stage(ss)) = Pslope;
       elseif dd == 4 
         stimdur68 = table(frequency, percenthigh);
         stimdur68_fit = table(xfreq,yfreq)
         stimdur68_slope(dd,stage(ss)) = Pslope;
                
       end   
       
           results = [results; frequency, percenthigh, ntrials, stagevar, durvar];
           resultsfit = [resultsfit; log(frequency), y, ntrials, stagevar, durvar];
           slopes = [slopes; b', Pslope, bias, stagevar(1), durvar(1)];
       
end


% plot target frequency against % correct for each stim duration with stderror

if stage(ss) == 1

    figure(1); clf
    errorbar(stimdur8.frequency, stimdur8.percenthigh, stdstim8(:,:,1),'rx-', 'Linewidth', 1.5 );
    hold on
    errorbar(stimdur20.frequency, stimdur20.percenthigh, stdstim20(:,:,1),'bx-', 'Linewidth', 1.5);
    errorbar(stimdur36.frequency, stimdur36.percenthigh, stdstim36(:,:,1),'kx-', 'Linewidth', 1.5);
    errorbar(stimdur68.frequency, stimdur68.percenthigh, stdstim68(:,:,1), 'gx-', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''All Harmonics'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off

    figure(11); clf
    plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
    hold on
    plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
    plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
    plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go: ', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''All Harmonics'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off
    
elseif stage(ss) == 2
    
    figure(2); clf
    errorbar(stimdur8.frequency, stimdur8.percenthigh, stdstim8(:,:,2),'rx-', 'Linewidth', 1.5 );
    hold on
    errorbar(stimdur20.frequency, stimdur20.percenthigh, stdstim20(:,:,2),'bx-', 'Linewidth', 1.5);
    errorbar(stimdur36.frequency, stimdur36.percenthigh, stdstim36(:,:,2),'kx-', 'Linewidth', 1.5);
    errorbar(stimdur68.frequency, stimdur68.percenthigh, stdstim68(:,:,2), 'gx-', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''Low Harmonics'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off

    figure(12); clf
    plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
    hold on
    plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
    plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
    plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go: ', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''Low Harmonics'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off
    
elseif stage(ss) == 3
    
    figure(3); clf
    errorbar(stimdur8.frequency, stimdur8.percenthigh, stdstim8(:,:,3),'rx-', 'Linewidth', 1.5 );
    hold on
    errorbar(stimdur20.frequency, stimdur20.percenthigh, stdstim20(:,:,3),'bx-', 'Linewidth', 1.5);
    errorbar(stimdur36.frequency, stimdur36.percenthigh, stdstim36(:,:,3),'kx-', 'Linewidth', 1.5);
    errorbar(stimdur68.frequency, stimdur68.percenthigh, stdstim68(:,:,3), 'gx-', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''High Harmonics'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off

    figure(13); clf
    plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
    hold on
    plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
    plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
    plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go: ', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''High Harmonics'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off

elseif stage(ss) == 4
    
    figure(4); clf
    errorbar(stimdur8.frequency, stimdur8.percenthigh, stdstim8(:,:,4),'rx-', 'Linewidth', 1.5 );
    hold on
    errorbar(stimdur20.frequency, stimdur20.percenthigh, stdstim20(:,:,4),'bx-', 'Linewidth', 1.5);
    errorbar(stimdur36.frequency, stimdur36.percenthigh, stdstim36(:,:,4),'kx-', 'Linewidth', 1.5);
    errorbar(stimdur68.frequency, stimdur68.percenthigh, stdstim68(:,:,4), 'gx-', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''All Harmonics, Random Phase'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off

    figure(14); clf
    plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
    hold on
    plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
    plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
    plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go: ', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''All Harmonics, Random Phase'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off
elseif stage(ss) == 5
    
    figure(5); clf
    errorbar(stimdur8.frequency, stimdur8.percenthigh, stdstim8(:,:,5),'rx-', 'Linewidth', 1.5 );
    hold on
    errorbar(stimdur20.frequency, stimdur20.percenthigh, stdstim20(:,:,5),'bx-', 'Linewidth', 1.5);
    errorbar(stimdur36.frequency, stimdur36.percenthigh, stdstim36(:,:,5),'kx-', 'Linewidth', 1.5);
    errorbar(stimdur68.frequency, stimdur68.percenthigh, stdstim68(:,:,5), 'gx-', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''High Harmonics, Random Phase'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off

    figure(15); clf
    plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
    hold on
    plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
    plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
    plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go: ', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''High Harmonics, Random Phase'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off
    
    if stage(ss) == 6
    
    figure(6); clf
    errorbar(stimdur8.frequency, stimdur8.percenthigh, stdstim8(:,:,1),'rx-', 'Linewidth', 1.5 );
    hold on
    errorbar(stimdur20.frequency, stimdur20.percenthigh, stdstim20(:,:,1),'bx-', 'Linewidth', 1.5);
    errorbar(stimdur36.frequency, stimdur36.percenthigh, stdstim36(:,:,1),'kx-', 'Linewidth', 1.5);
    errorbar(stimdur68.frequency, stimdur68.percenthigh, stdstim68(:,:,1), 'gx-', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''High Harmonics, Random Phase, No Spectral Cues'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off

    figure(16); clf
    plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
    hold on
    plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
    plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
    plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go: ', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''High Harmonics, Random Phase, No Spectral Cues'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off
    
elseif stage(ss) == 7
    
    figure(7); clf
    errorbar(stimdur8.frequency, stimdur8.percenthigh, stdstim8(:,:,2),'rx-', 'Linewidth', 1.5 );
    hold on
    errorbar(stimdur20.frequency, stimdur20.percenthigh, stdstim20(:,:,2),'bx-', 'Linewidth', 1.5);
    errorbar(stimdur36.frequency, stimdur36.percenthigh, stdstim36(:,:,2),'kx-', 'Linewidth', 1.5);
    errorbar(stimdur68.frequency, stimdur68.percenthigh, stdstim68(:,:,2), 'gx-', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''All Harmonics, No Spectral Cues'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off

    figure(17); clf
    plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
    hold on
    plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
    plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
    plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go: ', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''All Harmonics, No Spectral Cues'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off
    
elseif stage(ss) == 8
    
    figure(8); clf
    errorbar(stimdur8.frequency, stimdur8.percenthigh, stdstim8(:,:,3),'rx-', 'Linewidth', 1.5 );
    hold on
    errorbar(stimdur20.frequency, stimdur20.percenthigh, stdstim20(:,:,3),'bx-', 'Linewidth', 1.5);
    errorbar(stimdur36.frequency, stimdur36.percenthigh, stdstim36(:,:,3),'kx-', 'Linewidth', 1.5);
    errorbar(stimdur68.frequency, stimdur68.percenthigh, stdstim68(:,:,3), 'gx-', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''Low Harmonics, No Spectral Cues'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off

    figure(18); clf
    plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
    hold on
    plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
    plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
    plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go: ', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''Low Harmonics, No Spectral Cues'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off
 
elseif stage(ss) == 9
    
    figure(9); clf
    errorbar(stimdur8.frequency, stimdur8.percenthigh, stdstim8(:,:,4),'rx-', 'Linewidth', 1.5 );
    hold on
    errorbar(stimdur20.frequency, stimdur20.percenthigh, stdstim20(:,:,4),'bx-', 'Linewidth', 1.5);
    errorbar(stimdur36.frequency, stimdur36.percenthigh, stdstim36(:,:,4),'kx-', 'Linewidth', 1.5);
    errorbar(stimdur68.frequency, stimdur68.percenthigh, stdstim68(:,:,4), 'gx-', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''High Harmonics, No Spectral Cues'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off

    figure(19); clf
    plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
    hold on
    plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
    plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
    plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go: ', 'Linewidth', 1.5);
    title('Pitch Discrimination Performance: ''High Harmonics, No Spectral Cues'' Condition')
    xlabel('Target Frequency (Hz)')
    ylabel('Percentage "High" Response')
    xlim([245 255])
    ylim([0 100])
    plot([245 255],[50 50],'k--');
    plot([250 250],[0 100],'k--');
    legend('8ms', '20ms', '36ms', '68ms')
    hold off
    
end

end

save(['_differentstagesresults.mat'],'results','resultsfit','slopes')
end



 
  


