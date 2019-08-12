drs = 'Z:\Whole cell recordings\Img_5_raw_data_tifs' ;
drSave = 'Z:\Whole cell recordings\Img_5_processed';
n = 1;
total_blocks = 0;
sampleRate = 500;
first_frame = 1;
last_frame = 112899;
%discard = 5;

%     dr = uigetdir;
%     [fnSave, drSave] = uiputfile([dr filesep '..' filesep 'dataset.mat']);
 for i = 1:n
    
    dr = drs(i, :);
%     block_id = block_ids(i);
    D = dir([dr filesep '*.tif'])';

    fsize = [D.bytes];
    select = abs(fsize-median(fsize))<0.1*median(fsize);
    D = D(select);
    fNames = sort_nat({D.name});
    % fNames = fNames(discard*sampleRate:end);
    fNames = fNames(first_frame:last_frame);
    clear D
    
    A = tiffread2([dr filesep fNames{1}]);
    preview = A.data;
    boundX = 1:size(preview,2);
    boundY = 1:size(preview,1);

    blocklength = min(length(fNames), ceil(length(fNames)/ ceil(length(fNames)/50000)));
    blockstarts = 1:blocklength:length(fNames); 
    
%     blocklength = 50000;
%     blockstarts = 1;
    
%     if do_all_blocks
%         last_block = length(blockstarts);
%     end
    for blockN = total_blocks + 1:length(blockstarts)
        total_blocks = total_blocks + 1;
        disp(['Loading datasetblock' int2str(total_blocks)]);
      
%         try
%             load([drSave filesep 'datasetblock' int2str(total_blocks) '.mat'])
%         catch
        
            data = zeros([length(boundY) length(boundX) blocklength], 'single');
            for fnum = blockstarts(blockN):min(length(fNames), blockstarts(blockN)+blocklength-1)
                if ~mod(fnum,10000)
                    disp(['reading file:' int2str(fnum)])
                end
                A = tiffread2([dr filesep fNames{fnum}]);
                data(:,:,fnum-blockstarts(blockN)+1) = single(A.data(boundY,boundX));
            end
        
%         end
        
        disp('performing global alignment... ');
        refG = median(double(data(:,:,1:500)),3);
        alignShifts = nan(2,size(data,3));
        refFFT = fft2(refG);
        pobj = parpool;
        parfor f = 1:size(data,3)
            if ~mod(f,5000)
                disp(['registering frame: ' int2str(f)])
            end
            %align data
            [alignOut, G] = dftregistration(refFFT,fft2(data(:,:,f)),4);
            data(:,:,f) = real(ifft2(G));
            alignShifts(:,f) = alignOut(3:4);
        end
        delete(pobj);
        disp('Done alignment')
    
        %save data
        figure, plot(reshape(mean(mean(data, 1), 2), [1, size(data, 3)]))
        title(['Block ' int2str(blockN) ' average intensity over time'])
        savefast([drSave filesep 'datasetblock' int2str(blockN) '.mat'], 'data', 'sampleRate');
        
    end

 end