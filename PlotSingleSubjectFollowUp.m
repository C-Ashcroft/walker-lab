function results = PlotSingleSubjectFollowUp(name)

name 
cd([pwd '/TrialFollowUp/' name '/']) %change to the directory where your data is
 
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
slopes = [];

for ii=1:length(files) %for each file in our array
    rawdata(ind).file= files(ii).name; %enter these into 'rawdata'
    load([pwd '/TrialFollowUp/' name '/' rawdata(ind).file]) %load the file
    rawdata(ind).rawdata = d;
    DatTab =[DatTab; d];
    ind=ind+1;

end

%structure of data:
%  'Correct Pitch (1=high,0=low)' 'Was subject Correct (1,0)' 'Correction Trial(0=CT)' 'Target F0' 'Ref F0' 'Stage' 'Dur Ref' 'Dur Tar' 'Onset Delay' 'Level Tar' 'ISI' 'Ramp dur s' 'Sampling rate' 'delay' 'Timeout' 'Human''s response (1=high,0=low)' 'Reaction time'

stage = [];
stage = unique(DatTab(:,6));
results=[];
resultsfit=[];
slopes=[];

for ss = 1:length(stage)

    dur = unique(DatTab(:,7));

    % identify rows of data for duration 1 (8ms) to duration 4 (68ms)
    for dd = 1:length(dur)  
        finddurdata = [];
        finddurdata = find(DatTab(:,7) == dur(dd) & DatTab(:,6) == stage(ss));

        targets = unique(DatTab(:,4));
        percentcorrect = [];
        frequency = [];
        ntrials = [];
        stagevar= [];
        durvar = [];

        for tt = 1:length(targets) % get percent correct on each target pitch
            findtar = [];
            findtar = find(DatTab(:,4) == targets(tt));
            findtrials = intersect(findtar,finddurdata);
            frequency = [frequency; targets(tt)];
            percentcorrect = [percentcorrect ; mean(DatTab(findtrials,16))*100]; 
            ntrials = [ntrials ; length(findtrials)];
            stagevar = [stagevar ; stage(ss)];
            durvar = [durvar ; dur(dd)];
        end
    
        
        %calculate curve fit
        [b,dev,stats] = glmfit(log(frequency),[round(percentcorrect.*ntrials./100) ntrials],'binomial','logit'); %run a probit regression to produce a fit to the data
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
            stimdur8 = table(frequency, percentcorrect, ntrials, stagevar)
            stimdur8_fit = table(xfreq,yfreq)
            stimdur8_slope(dd,stage(ss)) = Pslope;
        elseif dd == 2
            stimdur20 = table(frequency, percentcorrect, ntrials, stagevar)
            stimdur20_fit = table(xfreq,yfreq)
            stimdur20_slope(dd,stage(ss)) = Pslope;
        elseif dd == 3
            stimdur36 = table(frequency, percentcorrect, ntrials, stagevar)
            stimdur36_fit = table(xfreq,yfreq)
            stimdur36_slope(dd,stage(ss)) = Pslope;
        elseif dd == 4 
            stimdur68 = table(frequency, percentcorrect, ntrials, stagevar)
            stimdur68_fit = table(xfreq,yfreq)
            stimdur68_slope(dd,stage(ss)) = Pslope;   
       end
       
       results = [results; frequency, percentcorrect, ntrials, stagevar, durvar];
       resultsfit = [resultsfit; log(frequency), y, ntrials, stagevar, durvar];
       slopes = [slopes; b', Pslope, bias, stagevar(1), durvar(1)];
       
    end

    % plot target frequency against % correct for each stim duration

    if stage(ss) == 1

        figure(5); clf
        plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
        hold on
        plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
        plot(stimdur36.frequency, stimdur36.percentcorrect,'kx-', 'Linewidth', 1.5);
        plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, All Harmonics')
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off
        
        figure(15); clf
        plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
        hold on
        plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
        plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
        plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go:', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, All Harmonics')
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off
        
    elseif stage(ss) == 2

        figure(6); clf
        plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
        hold on
        plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
        plot(stimdur36.frequency, stimdur36.percentcorrect,'kx-', 'Linewidth', 1.5);
        plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, Low Harmonics')
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off
        
        figure(16); clf
        plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
        hold on
        plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
        plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
        plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go:', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, Low Harmonics')
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off


    elseif stage(ss) == 3

        figure(7); clf
        plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
        hold on
        plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
        plot(stimdur36.frequency, stimdur36.percentcorrect,'kx-', 'Linewidth', 1.5);
        plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, High Harmonics')
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off
        
        figure(17); clf
        plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
        hold on
        plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
        plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
        plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go:', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, High Harmonics')
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off

    elseif stage(ss) == 4

        figure(8); clf
        plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
        hold on
        plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
        plot(stimdur36.frequency, stimdur36.percentcorrect,'kx-', 'Linewidth', 1.5);
        plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, All Harmonics, Random Phase')
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off
        
        figure(18); clf
        plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
        hold on
        plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
        plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
        plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go:', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, All Harmonics, Random Phase')
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off
    
        elseif stage(ss) == 5

        figure(9); clf
        plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
        hold on
        plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
        plot(stimdur36.frequency, stimdur36.percentcorrect,'kx-', 'Linewidth', 1.5);
        plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, High Harmonics, Random Phase')
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off
        
        figure(19); clf
        plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
        hold on
        plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
        plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
        plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go:', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, High Harmonics, Random Phase')
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off


    elseif stage(ss) == 6

        figure(10); clf
        plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
        hold on
        plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
        plot(stimdur36.frequency, stimdur36.percentcorrect,'kx-', 'Linewidth', 1.5);
        plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, High Harmonics, Random Phase (Rand)' )
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off
        
        figure(20); clf
        plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
        hold on
        plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
        plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
        plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go:', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, High Harmonics, Random Phase (Rand)' )
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off

    elseif stage(ss) == 7

        figure(11); clf
        plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
        hold on
        plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
        plot(stimdur36.frequency, stimdur36.percentcorrect,'kx-', 'Linewidth', 1.5);
        plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, All Harmonics (Rand)' )
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off
        
        figure(21); clf
        plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
        hold on
        plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
        plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
        plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go:', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, All Harmonics (Rand)' )
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off

      elseif stage(ss) == 8

        figure(12); clf
        plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
        hold on
        plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
        plot(stimdur36.frequency, stimdur36.percentcorrect,'kx-', 'Linewidth', 1.5);
        plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, Low Harmonics (Rand)' )
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off
        
        figure(22); clf
        plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
        hold on
        plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
        plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
        plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go:', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, Low Harmonics (Rand)' )
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off


    elseif stage(ss) == 9

        figure(13); clf
        plot(stimdur8.frequency, stimdur8.percentcorrect,'rx-', 'Linewidth', 1.5 );
        hold on
        plot(stimdur20.frequency, stimdur20.percentcorrect,'bx-', 'Linewidth', 1.5);
        plot(stimdur36.frequency, stimdur36.percentcorrect,'kx-', 'Linewidth', 1.5);
        plot(stimdur68.frequency, stimdur68.percentcorrect,'gx-', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, High Harmonics (Rand)' )
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off
        
        figure(23); clf
        plot(stimdur8_fit.xfreq,stimdur8_fit.yfreq*100,'ro:', 'Linewidth', 1.5);
        hold on
        plot(stimdur20_fit.xfreq,stimdur20_fit.yfreq*100,'bo:', 'Linewidth', 1.5);     
        plot(stimdur36_fit.xfreq,stimdur36_fit.yfreq*100,'ko: ', 'Linewidth', 1.5);
        plot(stimdur68_fit.xfreq,stimdur68_fit.yfreq*100,'go:', 'Linewidth', 1.5);
        title('Effect of Varying Stimulus Duration of Pitch Discrimination: Single Subject, High Harmonics (Rand)' )
        xlabel('Target Frequency (Hz)')
        ylabel('Percentage "High" Response')
        xlim([245 255])
        ylim([0 100])
        plot([245 255],[50 50],'k--');
        plot([250 250],[0 100],'k--');
        legend('8ms','20ms','36ms','68ms')
        hold off

    end

   

end

save([pwd '/SlopeDataFollowUp/' name '_results.mat'],'results','resultsfit','slopes')

end 


 



