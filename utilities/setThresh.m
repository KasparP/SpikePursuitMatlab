function thresh = setThresh(Y, pDes, maxSpikes)

if nargin<2
    pDes = 0.999;
end
if nargin<3
    maxSpikes = 500;
end

minThresh = prctile(Y, 100*(1-(maxSpikes/length(Y))));

%We want a threshold that corresponds to ((true positives/total detections)>pDes)
S.mu = [min(Y) ; max(Y)];
Sigma = std(Y);
S.Sigma = reshape([Sigma;Sigma], 1, 1, 2); %Sigma; %[Sigma 0; 0 Sigma];
S.ComponentProportion = [ 0.95 ; 0.05];
GMM = fitgmdist(Y,2, 'Start', S);  %'SharedCovariance', true

% iter = 0; 
% while iter<10
%     warning('');
%     GMM = fitgmdist(Y,2, S);
%     [warnMsg, warnId] = lastwarn;
%     if strfind(warnId, 'gmdistribution')
%         iter = iter+1;
%         if iter>9
%             error('couldn''t fit GMM')
%         end
%     end
% end
x = linspace(min(Y), max(Y), 1000);

nGauss1 = GMM.ComponentProportion(1)*(1-normcdf(x, GMM.mu(1),GMM.Sigma(1)));
nGauss2 = GMM.ComponentProportion(2)*(1-normcdf(x,GMM.mu(2),GMM.Sigma(1)));
if GMM.mu(1)<GMM.mu(2)
    tmp = nGauss1;
    nGauss1 = nGauss2; nGauss2 = tmp; clear tmp;
end
p = nGauss1./(nGauss1+nGauss2);
figure, plot(x,p)
thresh = x(find(p>pDes,1,'first'));

nTruePos = nGauss1(find(p>pDes,1,'first'))*length(Y);
if isempty(thresh)
    error('We couldn''t find any spikes')
elseif nTruePos<30
    warning(['You likely only have less than ' num2str(nTruePos) ' true spikes above threshold']);
end
  
% pdf1 = GMM.ComponentProportion(1)*(normpdf(x, GMM.mu(1),GMM.Sigma(1)));
% pdf2= GMM.ComponentProportion(2)*(normpdf(x, GMM.mu(2),GMM.Sigma(1)));
%figure, plot(x, pdf1); hold on, plot(x, pdf2)

thresh = max(thresh, minThresh);
end