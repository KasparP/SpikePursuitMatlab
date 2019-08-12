function output = spikePursuit(dr, fns, ROIs, cell_ids)
%The only required argument is ROIs, a [Y-by-X-by-#ROIs] binary array containing the initial ROIs. 
%The easiest way to call this function is spikePursuit([],[], ROIs)

%todo:
%intialize with a learned temporal filter (~line 70)
%allow for multiple templates in denoiseSpikes (k means clustering?)

%If the spike-triggered average is clearly noise (e.g. the points adjacent
%to the center point don't vary significantly from zero), set a flag/abort early

if isempty(dr)
    [fns, dr] = uigetfile('*.mat', 'Select your data chunks', 'multiselect', 'on');
    if ~iscell(fns)
        fns = {fns};
    end
end

%options
opts.doCrossVal = false; %cross-validate to optimize regression regularization parameters?
opts.doGlobalSubtract = false; %optional; did not improve performance in simulations unless there was very serious global noise 
        
opts.contextSize = 50; % 65; %number of pixels surrounding the ROI to use as context
opts.censorSize = 12; %number of pixels surrounding the ROI to censor from the background PCA; roughly the spatial scale of scattered/dendritic neural signals, in pixels.
opts.nPC_bg = 8; %number of principle components used for background subtraction
opts.tau_lp = 3; %time window for lowpass filter (seconds); signals slower than this will be ignored
opts.tau_pred = 1; % time window in seconds for high pass filtering to make predictor for regression
opts.sigmas = [1 1.5 2]; %spatial smoothing radius imposed on spatial filter;
opts.nIter = 5; %number of iterations alternating between estimating temporal and spatial filters.
opts.localAlign = false; 
opts.globalAlign = true;
opts.highPassRegression = false; %regress on a high-passed version of the data. Slightly improves detection of spikes, but makes subthreshold unreliable.


for fn_ix = 1:length(fns)

    disp(['Loading data batch: ' fns{fn_ix}])
    struct = load([dr filesep fns{fn_ix}], 'data', 'sampleRate');
    dataAll = struct.data;
    sampleRate = struct.sampleRate;
    disp(['sampleRate: ' int2str(sampleRate)])
    
    opts.windowLength =  sampleRate*0.02; %window length for spike templates
    
    %Compute global PCs with ROIs masked out 
    if opts.doGlobalSubtract
        mask = ~imdilate(any(ROIs,3), strel('disk', opts.censorSize));
        data = reshape(dataAll, [], size(dataAll, 3));
        data = double(data(mask(:),:));
        disp('Performing highpass filtering');
        tic
        data = highpassVideo(data', 1/opts.tau_lp, sampleRate)'; %takes ~2-3 minutes
        toc
        disp('Performing PCA...')
        tic
        [~,~,Vg_hp] = svds(data, opts.nPC_bg); %takes ~2-3 minutes
        toc
        Vg_pred = highpassVideo(Vg_hp, 1/opts.tau_pred, sampleRate); %filter Vg
    end
    
    opts.windowLength = sampleRate*0.02;
    if nargin < 4
         cell_ids = 1:size(ROIs,3);
    end
    
    for cellN = 1:size(ROIs,3)
        disp(['Processing cell:' int2str(cellN)]);
        tic
        bw = ROIs(:,:,cellN);
        
        %extract relevant region and align
        bwexp = imdilate(bw, ones(opts.contextSize));
        Xinds = find(any(bwexp,1),1,'first'):find(any(bwexp,1),1,'last');
        Yinds= find(any(bwexp,2),1,'first'):find(any(bwexp,2),1,'last');
        bw = bw(Yinds,Xinds);
        notbw = ~imdilate(bw, strel('disk', opts.censorSize));
        data = dataAll(Yinds, Xinds, :);   
        bw = logical(bw);
        notbw = logical(notbw);

        ref = median(double(data(:,:,1:500)),3);
%         figure('name', 'ROI selection'),
%         subplot(1,3,1); imagesc(ref); axis image; xlabel('mean Intensity')
%         subplot(1,3,2), imagesc(bw); axis image; xlabel('intial ROI');
%         subplot(1,3,3), imagesc(notbw); axis image; xlabel('background');
        
        if opts.localAlign
            refFFT = fft2(ref);
            drawnow;
            pobj = parpool;
            parfor f = 1:size(data,3)
                if ~mod(f,5000)
                    disp(['registering frame: ' int2str(f)])
                end
                %align data
                [~, G] = dftregistration(refFFT,fft2(data(:,:,f)),4);
                data(:,:,f) = real(ifft2(G));
            end
            delete(pobj);
            disp('Done alignment')
        end
        
        output.meanIM = mean(data,3);
        data = reshape(data, [], size(data,3));
        
        data = double(data);
        data = double(data-mean(data,2)); %mean subtract
        data = data-mean(data,2); %do it again because of numeric issues
        
        %remove low frequency components
        data_hp = highpassVideo(data', 1/opts.tau_lp, sampleRate)';
        data_lp = data-data_hp;
       
        if opts.highPassRegression
            data_pred =  highpassVideo(data', 1/opts.tau_pred, sampleRate)';
        else
            data_pred = data_hp;
        end

        close all;
        
        t = nanmean(double(data_hp(bw(:),:)),1); %initial trace is just average of ROI
        t = t-mean(t);
        output.t = t;
        
        %remove any variance in trace that can be predicted from the background PCs
        [Ub,Sb,Vb] = svds(double(data_hp(notbw(:),:)), opts.nPC_bg);
        b = regress(t', Vb);
        t = (t'-Vb*b); %initial trace
        
        %Initial spike estimate       
        [Xspikes, spikeTimes, guessData, output.rawROI.falsePosRate, output.rawROI.detectionRate, output.rawROI.templates, low_spk] = denoiseSpikes(-t', opts.windowLength, sampleRate,true, 100);
        Xspikes = -Xspikes;
        output.rawROI.X = t';
        output.rawROI.Xspikes = Xspikes; 
        output.rawROI.spikeTimes = spikeTimes;
        output.rawROI.spatialFilter = bw;
        output.rawROI.X = output.rawROI.X.*(mean(t(output.rawROI.spikeTimes))/mean(output.rawROI.X(output.rawROI.spikeTimes)));%correct shrinkage
        
        %prebuild the regression matrix
        pred = [ones(1,size(data_pred,2)); reshape(imgaussfilt(reshape(data_pred, [size(ref) size(data,2)]), 1.5), size(data))]'; %generate a predictor for ridge regression
 
        % To do: if not enough spikes, take spatial filter from previous block
                               
        % Cross-validation of regularized regression parameters
        lambdamax = norm(pred(2:end,:),'fro').^2;
        lambdas = lambdamax*logspace(-4,-2,3); %if you want multiple values of lambda
        %lambdas = lambdas(2); %fixing default value for speed
        I0 = eye(size(pred,2)); I0(1)=0;

        if opts.doCrossVal
            num_batches = 3;
            batchsize = floor(size(data,2)/num_batches);
            for batch = 1:num_batches
                disp(['crossvalidating lambda, batch ' int2str(batch) ' of ' int2str(num_batches)])
                select = false(size(guessData));
                select((batch-1)*batchsize + (1:batchsize)) = true;
                for s_ix = 1:length(opts.sigmas)
                    pred = [ones(1,size(data_pred,2)); reshape(imgaussfilt(reshape(data_pred, [size(ref) size(data_pred,2)]), opts.sigmas(s_ix)), size(data_pred))]';
                    for l_ix = 1:length(lambdas)
                        kk2= (pred(~select, :)'*pred(~select, :)+lambdas(l_ix)*I0)\pred(~select, :)';
                        weights = kk2*(guessData(~select))'; %regression
                        corrs2(l_ix, s_ix, batch) = corr(pred(select, :)*weights, guessData(select)');   %#ok<AGROW>
                    end
                end
            end
            [l_max, s_max] = find(nanmean(corrs2, 3) == nanmax(nanmax(nanmean(corrs2, 3))));
            opts.lambda = lambdas(l_max);
            opts.lambda_ix = l_max;
            opts.sigma = opts.sigmas(s_max);
            if isempty(s_max)
                disp('a cell had no spikes.... continuing')
                continue
            end
        else %fix the values:
            s_max = 2;
            l_max = 3;
            opts.lambda = lambdas(l_max);
            opts.sigma = opts.sigmas(s_max);
            opts.lambda_ix = l_max;
        end
        
        selectPred = true(1,size(data,2));
        if opts.highPassRegression
            selectPred([1:(sampleRate/2+1) (end-sampleRate/2):end]) = false; %discard data at edges to avoid any filtering artefacts; optional
        end
        
        pred = [ones(1,size(data_pred,2)); reshape(imgaussfilt(reshape(data_pred, [size(ref) size(data_pred,2)]), opts.sigmas(s_max)), size(data_pred))]';
        recon = [ones(1,size(data_hp,2)); reshape(imgaussfilt(reshape(data_hp, [size(ref) size(data_hp,2)]), opts.sigmas(s_max)), size(data_hp))]';
        kk = (pred(selectPred,:)'*pred(selectPred,:) +lambdas(l_max)*I0)\pred(selectPred,:)';
        
        for iter = 1:opts.nIter
            doPlot = false;
            if iter==opts.nIter
                doPlot = true;
            end
            
            disp('Identifying spatial filters') %identify spatial filters with regularized regression
            gD = guessData(selectPred); select = gD~=0;
            weights = kk(:,select)*gD(select)'; %regression
            X = double(recon*weights)';
            X = X-mean(X);
            
            spatialFilter = imgaussfilt(reshape(weights(2:end), size(ref)),opts.sigmas(s_max));

            if iter < opts.nIter
                b = regress(X', Vb);  %remove background contamination in intermediate iterations. This probably helps but not tested thoroughly.
                if doPlot
                    figure('name', 'Denoised trace vs background'), plot(X), hold on, plot(Vb*b)
                end
                X = X-(Vb*b)';
            else
                if opts.doGlobalSubtract
                    b = regress(X', Vg_hp);
                    X = X-(Vg_hp*b)';

                    b = regress(X', Vg_pred);
                    X = X-(Vg_pred*b)';
                    output.Vg = Vg_hp; %global background components
                    output.b = b; %weights of global background components that were subtracted to produce y
                end
            end
            X = X.*(mean(t(spikeTimes))/mean(X(spikeTimes)));%correct shrinkage

            %generate the new trace and the new denoised trace
            [Xspikes, spikeTimes, guessData, falsePosRate, detectionRate, templates, ~] = denoiseSpikes(-X, opts.windowLength, sampleRate,doPlot);
        end
        
        %ensure that the maximum of the spatial filter is within the ROI
        IMcorr = corr(-guessData', pred(:, 2:end));
        maxCorrInROI = max(IMcorr(bw(:)));
        if any(IMcorr(~bw(:))>maxCorrInROI)
            output.passedLocalityTest = false;
        else
            output.passedLocalityTest = true;
        end

        %compute SNR
        selectSpikes = false(length(Xspikes),1); selectSpikes(spikeTimes) = true;
        signal = mean(Xspikes(selectSpikes));
        noise = std(Xspikes(~selectSpikes));
        snr = signal/noise;
        output.snr = snr;
        
        %output
        output.y = X;
        output.yFilt = -Xspikes;
        output.ROI = [Xinds([1 end])'  Yinds([1 end])'];
        output.ROIbw = bw;
        output.spatialFilter = spatialFilter;
        output.falsePosRate = falsePosRate;
        output.detectionRate = detectionRate;
        output.templates = templates;
        output.spikeTimes = spikeTimes;
        output.opts = opts;
        output.F0 = nanmean(double(data_lp(bw(:),:))+output.meanIM(bw(:)),1);
        output.dFF = X(:)./output.F0(:);
        output.rawROI.dFF = output.rawROI.X(:)./output.F0(:);
        output.Vb = Vb; %local background components
        output.low_spk = low_spk;
        
        drawnow;
        save([dr filesep fns{fn_ix}(1:end-4) 'Cell' int2str(cell_ids(cellN))], 'output');
        toc
    end
end
end

function videoFilt = highpassVideo(video, freq, sampleRate)
        normFreq = freq/(sampleRate/2);
        [b,a] = butter(3,normFreq, 'high');
        videoFilt = filtfilt(b,a,video); 
end