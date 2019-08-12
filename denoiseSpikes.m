function [datafilt, spikeTimes, guessData, falsePosRate, detectionRate, templates, low_spk] = denoiseSpikes(data, windowLength, sampleRate, doPlot, doClip)
if nargin<3
    sampleRate = 500;
end
if nargin<4
    doPlot = false;
end
if nargin<5
    doClip = 150; %a limit on the number of spikes used to make the template
end

%highpass filter and threshold
[bb,aa] = butter(1, 1/(sampleRate/2), 'high'); % 1Hz filter
dataHP = filtfilt(bb,aa,data);

pks = findpeaks(dataHP);
[thresh, ~, ~, low_spk] = getThresh(pks,0.25, doClip);
[~,locs] = findpeaks(dataHP', 'MinPeakHeight',thresh);
    
%peak-triggered average
window = -windowLength:windowLength;
locs = locs(locs>(-window(1)+1) & locs<(length(data)-window(end)));
PTD = data(locs+repmat(window, length(locs),1));
PTA = mean(PTD,1); %peak-triggered average

%matched filter
datafilt = whitenedMatchedFilter(data,locs, window);

%spikes detected after filter
pks2 = findpeaks(datafilt);
[thresh2, falsePosRate, detectionRate, ~] = getThresh(pks2);
[~,spikeTimes] = findpeaks(datafilt, 'MinPeakHeight',thresh2);

guessData = zeros(size(data));
guessData(spikeTimes) = 1;
guessData = conv(guessData,PTA, 'same');

%filtering shrinks the data;
%rescale so that the mean value at the peaks is same as in the input
datafilt = datafilt.*(mean(data(spikeTimes))./mean(datafilt(spikeTimes)));

%output templates
templates = PTA;

%plot
if doPlot
    figure,
    subplot(2,1,1)
    hist(pks, 500), hold on, plot([thresh thresh], get(gca, 'ylim'), 'r'); xlabel('raw data')
    subplot(2,1,2)
    hist(pks2, 500), hold on, plot([thresh2 thresh2], get(gca, 'ylim'), 'r'); xlabel('after matched filter')
    
    figure('name', 'Peak-triggered average'), plot(PTD', 'color', [0.5 0.5 0.5]), hold on, plot(PTA, 'k', 'linewidth', 2)
    
    figure('name', 'Data before and after filtering (61Hz highpass)')
    ax1 = subplot(2,1,1);
    plot(data),
    hold on, scatter(locs, max(datafilt)*1.1*ones(size(locs)), 'markeredgecolor', 'r')
    hold on, scatter(spikeTimes, max(datafilt)*ones(size(spikeTimes)), 'markeredgecolor', 'g')
    ax2 = subplot(2,1,2);
    plot(datafilt),
    hold on, scatter(locs, max(datafilt)*1.1*ones(size(locs)), 'markeredgecolor', 'r')
    hold on, scatter(spikeTimes, max(datafilt)*ones(size(spikeTimes)), 'markeredgecolor', 'g')
    linkaxes([ax1 ax2]);
end
end

function [thresh, falsePosRate, detectionRate, low_spk] = getThresh(pks, pnorm, doClip)
if nargin<2
    pnorm = 0.5; %sensitivity/selectivity tradeoff parameter; some number between 0 and 1 exclusive; larger values favor detections, smaller values favor rejecting false positives
end
spread = [min(pks)  max(pks)]; spread = spread+diff(spread).*[-0.05 0.05];
low_spk = false;
pts = linspace(spread(1), spread(2), 2001);
[f,xi] = ksdensity(pks, pts);
center = find(xi>median(pks),1,'first');
%[~,center] = max(f);

fmodel = [f(1:center) fliplr(f(1:(center-1)))];
if length(fmodel)<length(f)
    fmodel((length(fmodel)+1):length(f)) = min(fmodel);
else
    fmodel = fmodel(1:length(f));
end
%adjust the model so it doesn't exceed the data:
csf = cumsum(f)./sum(f); csmodel = cumsum(fmodel)./max(sum(f), sum(fmodel));
lastpt = find(csf(1:end-1)>csmodel(1:end-1)+eps & csf(2:end)<csmodel(2:end), 1, 'last');
if isempty(lastpt)
    lastpt = center;
end
fmodel(1:lastpt) = f(1:lastpt);
fmodel(lastpt:end) = min(fmodel(lastpt:end), f(lastpt:end));

csf = cumsum(f); csmodel = cumsum(fmodel);
%csf = max(csf,csmodel); %the cumulative distribution of the model should be less than or equal to the measurements

csf2 = csf(end)-csf;
csmodel2 = csmodel(end)-csmodel;

obj = csf2.^pnorm - csmodel2.^pnorm;
[~,maxind] = max(obj);
thresh = xi(maxind);
if sum(pks>thresh)<30
    low_spk = true;
    disp(['Very few spikes were detected at the desired sensitivity/specificity tradeoff. Adjusting threshold to take 30 largest spikes']);
    thresh = prctile(pks, 100*(1- 30/length(pks)));
elseif nargin>2 && ~isempty(doClip) && sum(pks>thresh)>doClip
    disp(['Selecting top ' int2str(doClip) ' spikes for template']);
    thresh = prctile(pks, 100*(1- doClip/length(pks)));
end
[~, ix] = min(abs(xi-thresh));
falsePosRate = csmodel2(ix)./csf2(ix);
detectionRate = (csf2(ix)-csmodel2(ix))./max(csf2-csmodel2);
end