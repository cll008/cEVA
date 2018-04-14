function savedFrameNums = cEVA(videoFile, streams, chansToPlotIdx, scrollWindowInSec, rewindFrames, forwardFrames)
% Key-press controlled frame-by-frame video viewer using default Matlab VideoReader. 
% Extract event frames using single frame scrolling and pseudo-frame-based rewind.
%
% Usage:
%   >>savedFrameNums = frameExtrtact(videoFile, rewindFrames, forwardFrames);
%
% Inputs: 
%   videoFile         - full filepath to video file
%   streams 	      - stream structure from loaded .xdf file 
%   chansRToPloIdx    - indices of channels to plot on EEG scroll plot
%   scrollWindowInSec - time length (s) of EEG to plot
%   rewindFrames      - number of frames to rewind
%	forwardFrames     - number of frames to skip forward
%
% Outputs:
%   savedFrameNums - vector containing extracted frame numbers
%
% See also:
%   VideoReader
%
% Author: 
%   Clement Lee & Makoto Miyakoshi, Swartz Center for Computational Neuroscience, 
%   Institute for Neural Computation, UC San Diego
%   
%
% History: 
%  11 Apr 2018 v1.1 by CL. Added input and support for [chansToPlotIdx]
%  09 Apr 2018 v1.0 by CL & MM. Added EEG scroll at UCLA.
%  02 Feb 2018 v0.1 by Clement. Created.
for streamIdx = 1:length(streams)
    if strcmp(streams{streamIdx}.info.name,'EGI NetAmp 0')
        EEGStream = streams{streamIdx};
    elseif strcmp(streams{streamIdx}.info.name,'VideoCaptureWrite')
        videoStream = streams{streamIdx};
    end
end

opts.datascale = 100; 

% Detect the outlier to exclude from the plotting.
plotdata = bsxfun(@minus, EEGStream.time_series, mean(EEGStream.time_series,2));
stdList = std(EEGStream.time_series,[],2);
channelToExcludeIdx = find(stdList>5*median(stdList));
plotdata(channelToExcludeIdx,:) = zeros(length(channelToExcludeIdx),size(plotdata,2));

% chansToPlotIdx
plotdata(setdiff([1:size(plotdata,1)],chansToPlotIdx),:) = []; % rm chans
nbchan = size(plotdata,1);

% Build fake EEG structure to make eegfiltnew() function work.
dummyEEG = eeg_emptyset;
dummyEEG.data = plotdata;
dummyEEG.xmin = 0;
dummyEEG.xmax = str2double(EEGStream.info.last_timestamp) - str2double(EEGStream.info.first_timestamp);
dummyEEG.pnts = str2double(EEGStream.info.sample_count);
dummyEEG.srate = str2double(EEGStream.info.nominal_srate);
dummyEEG = pop_eegfiltnew(dummyEEG, 1, 0);

plotdata = dummyEEG.data;
plotoffsets = (0:size(plotdata,1)-1)'*opts.datascale;
plotdata = bsxfun(@plus, plotdata, plotoffsets);

vidObj = VideoReader(videoFile);
videoFrameIdx = 0; frameTimings = [];
firstVideoFrame = videoStream.time_series(1); 
EEGStreamTimes = interp1(EEGStream.time_stamps, EEGStream.time_stamps, videoStream.time_stamps, 'nearest');
[~, vidCorrEEGFrameIdx] = intersect(EEGStream.time_stamps, EEGStreamTimes);

tic;
t2old=0;
fprintf('Reading the full video file to preconstruct frameTimings...\n')
while hasFrame(vidObj)
   	videoFrameIdx = videoFrameIdx+1;
	vidFrame = readFrame(vidObj);
    frameTimings(length(frameTimings)+1) = vidObj.CurrentTime;
    
    if mod(videoFrameIdx,100) == 0
        t2=toc;
        disp(horzcat('Frame number reached: ',num2str(videoFrameIdx), ' of ',num2str(vidObj.Duration*vidObj.FrameRate,'%.0f'),'. Estimated time remaining: ',num2str(datestr(seconds((((vidObj.Duration*vidObj.FrameRate)-videoFrameIdx)/100)*(t2-t2old)),'HH:MM:SS')),'. Time diff: ',num2str(t2-t2old)))
        t2old=t2;
    end
        
end
fprintf('Finished reading video.\n')

fprintf(['Press "q" to start and advance frame by frame\n'...
	'Press "g" to go to specific frame\n'...
    'Press "e" to enter frame number into savedFrameNums output vector\n'...
    'Press "d" to delete last saved frame number\n'...
    'Press "f" to fastforward by "forwardFrames" amount of frames \n'...
    'Press "r" to rewind by "rewindFrames" amount of frames \n'...
	'Press "c" to check frames saved in "savedFrameNums"\n'])

videoFrameIdx = firstVideoFrame; 
vidObj.CurrentTime = 0;
savedFrameNums = [];
lineHandles = [];
figure;
videoAxes = subplot(2,1,1);
lineAxes = subplot(2,1,2);
while hasFrame(vidObj)
    k = waitforbuttonpress;
    if k == 1 %checks if waitforbuttonpress was a key
        key = get(gcf,'CurrentCharacter'); % which key was pressed (if any)?
        if strcmp(key, 'q') % next frame
            videoFrameIdx = videoFrameIdx+1;
            vidObj.CurrentTime = frameTimings(videoFrameIdx);
			vidFrame = readFrame(vidObj);
			
		elseif strcmp(key, 'g') % go
			videoFrameIdx = input('Enter frame number to skip to > ');
			vidObj.CurrentTime = frameTimings(videoFrameIdx);
			vidFrame = readFrame(vidObj);
			
        elseif strcmp(key,'e') % enter frame to vector
            if ~ismember(videoFrameIdx, savedFrameNums) % only add 1 marker per frame
                savedFrameNums(length(savedFrameNums)+1) = videoFrameIdx;
                fprintf('Added frame number %d\n',videoFrameIdx)
            else 
                fprintf('Duplicate frame number not added. Press "c" key to double check.\n')
            end
        elseif strcmp(key,'d') % delete last frame saved
            savedFrameNums(end) = [];
            fprintf('Removed last frame added\n')
            
		elseif strcmp(key,'f') % fastforward
			videoFrameIdx = videoFrameIdx + forwardFrames;
			vidObj.CurrentTime = frameTimings(videoFrameIdx);
			vidFrame = readFrame(vidObj);
            
        elseif (strcmp(key,'r') && (videoFrameIdx - rewindFrames) > 0) % rewind
            videoFrameIdx = videoFrameIdx - rewindFrames;
            vidObj.CurrentTime = frameTimings(videoFrameIdx);
            vidFrame = readFrame(vidObj);
			
		elseif strcmp(key,'c') % print frames logged
            if isempty(savedFrameNums)
                fprintf('No frames added yet')
            else
                fprintf('%d\n', savedFrameNums)
            end
        
		else
            disp('Please check key bindings and instructions')
        end
		image(vidFrame, 'Parent', videoAxes);
        set(videoAxes,'xtick',[])
        set(videoAxes,'ytick',[])
        title(['Frame ' num2str(videoFrameIdx) ' of ' num2str(length(frameTimings))])

        % draw EEG scroll
        % EEGStreamTimes is in seconds
        plottime = linspace((EEGStreamTimes(videoFrameIdx-firstVideoFrame)-EEGStreamTimes(1)),...
            EEGStreamTimes(videoFrameIdx-firstVideoFrame)-EEGStreamTimes(1)+scrollWindowInSec, scrollWindowInSec*1000);
        currentEdgeInt1 = vidCorrEEGFrameIdx(videoFrameIdx-firstVideoFrame+1);
        currentEdgeInt2 = vidCorrEEGFrameIdx(videoFrameIdx-firstVideoFrame+1) + scrollWindowInSec*1000;
        
        if isempty(lineHandles)
            lineHandles = plot(lineAxes,plottime,plotdata(:,currentEdgeInt1:currentEdgeInt2-1));
            xlabel(lineAxes,'Time (s)','FontSize',12);
        else
            for k=1:min(length(lineHandles),size(plotdata,1))
                set(lineHandles(k),'Xdata',plottime, 'Ydata',plotdata(k,currentEdgeInt1:currentEdgeInt2-1)); end
        end
        
        % update the axis limit and tickmarks
        winBeg  = (EEGStreamTimes(videoFrameIdx-firstVideoFrame)-EEGStreamTimes(1)); % In second.
        winCent = (EEGStreamTimes(videoFrameIdx-firstVideoFrame)-EEGStreamTimes(1)) + scrollWindowInSec/2; % In second.
        winEnd  = (EEGStreamTimes(videoFrameIdx-firstVideoFrame)-EEGStreamTimes(1)) + scrollWindowInSec*1; % In second.

        axis(lineAxes,[winBeg, winEnd, -opts.datascale, nbchan*opts.datascale + opts.datascale]);
        set(lineAxes, 'XTick',[winBeg winCent winEnd], ...
            'XTickLabel',...
           {winBeg winCent winEnd},...
           'YTick',[]);
      drawnow
        
    end
end