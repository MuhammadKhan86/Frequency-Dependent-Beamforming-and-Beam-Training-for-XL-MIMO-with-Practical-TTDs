function Figure_3

close all; clear all; clc;

% parpool('HPCServer',50);
%%% parameters
params.centralFreq      = 100e9;                                            centralFreq    = params.centralFreq;
params.sampFreq         = 10e9;                                             sampFreq       = params.sampFreq;
params.nTx              = 512;                                              nTx            = params.nTx;
params.nTTD             = nTx;                                              nTTD           = params.nTTD;
params.waveLength       = physconst('LightSpeed')/centralFreq ;             waveLength     = params.waveLength;
params.antElemSpacing   = 1/2 * waveLength;                                 antElemSpacing = params.antElemSpacing;
params.maxBeamAngle     =  60;                                              maxBeamAngle   = params.maxBeamAngle;
params.minBeamAngle     = -60;                                              minBeamAngle   = params.minBeamAngle;

params.antArraySize     = nTx*antElemSpacing;

%%% TTD params
ttdTimeDelayRange = 15e-9;
ttdTimeDelayRes   = 3e-12;

beamDirection    = 30;
chanDirectionVec = -65:0.01:65;



%%% Subcarrier Frequencies
nCarr           = ceil((1/(sqrt(2)*1.782/nTx)* (sind(maxBeamAngle) - sind(minBeamAngle))) + 1);


fCarrAll        = (centralFreq + sampFreq / (nCarr-1) .* ((0:nCarr-1)-(nCarr - 1) / 2));
centralFreqRep  = repelem(centralFreq,1,nCarr);

timeDelayIdeal  = 0.5/(centralFreq*sampFreq) * (fCarrAll(1)*sind(minBeamAngle) - fCarrAll(end)*sind(maxBeamAngle));
fixedAngle      = asind(0.5/centralFreq * ((centralFreq - sampFreq/2)*sind(minBeamAngle) + (centralFreq + sampFreq/2)*sind(maxBeamAngle)));

optCodebookAngles = asind((2*centralFreqRep.*(centralFreqRep-fCarrAll)*timeDelayIdeal + centralFreqRep*sind(fixedAngle)) ./ (fCarrAll));
for locIdx = 1:length(chanDirectionVec)
  [~, optBeamIdx(locIdx)] = min(abs(optCodebookAngles-chanDirectionVec(locIdx)));
end


%%%% Frequency Dependent Beam Generation

%%% Parallel
timeDelayVecIdealPDPP         = timeDelayIdeal * (0:nTx-1);
timeDelayVecRangeConsPDPP     = min(abs(timeDelayVecIdealPDPP), abs(ttdTimeDelayRange)) .* sign(timeDelayVecIdealPDPP);
timeDelayVecResConsPDPP       = round(timeDelayVecIdealPDPP./ttdTimeDelayRes)*ttdTimeDelayRes;
timeDelayVecRangeResConsPDPP  = round(timeDelayVecRangeConsPDPP./ttdTimeDelayRes)*ttdTimeDelayRes;

%%% Serial
timeDelayVecIdealCDPP         = zeros(1, nTx);
timeDelayVecRangeConsCDPP     = zeros(1, nTx);
timeDelayVecResConsCDPP       = zeros(1, nTx);
timeDelayVecRangeResConsCDPP  = zeros(1, nTx);
for txIdx = 2:nTx
  timeDelayVecIdealCDPP(txIdx)        = timeDelayVecIdealCDPP(txIdx-1)        + timeDelayIdeal;
  timeDelayVecRangeConsCDPP(txIdx)    = timeDelayVecRangeConsCDPP(txIdx-1)    + min(abs(timeDelayIdeal), abs(ttdTimeDelayRange)) .* sign(timeDelayIdeal);
  timeDelayVecResConsCDPP(txIdx)      = timeDelayVecResConsCDPP(txIdx-1)      + round(timeDelayIdeal./ttdTimeDelayRes)*ttdTimeDelayRes;
  timeDelayVecRangeResConsCDPP(txIdx) = timeDelayVecRangeResConsCDPP(txIdx-1) + round(min(abs(timeDelayIdeal), abs(ttdTimeDelayRange)) .* sign(timeDelayIdeal)./ttdTimeDelayRes)*ttdTimeDelayRes;
end


%%% MSDPP
%%% Ideal
nBinaryBits               = ceil(log2(nTx));
idealDelays               = 2.^(0:1:nBinaryBits-1).'*timeDelayIdeal;
nIdealStages              = size(idealDelays,1);
idealScalingPerStage      = 2.^(0:nIdealStages-1);

for antIdx = 0:nTx-1
  ttdNumber = antIdx;
  for stageIdx = nIdealStages:-1:1
    idealCoeffs(antIdx+1, stageIdx) = floor(ttdNumber / idealScalingPerStage(stageIdx));
    ttdNumber = ttdNumber - idealCoeffs(antIdx+1, stageIdx) * idealScalingPerStage(stageIdx);
  end
end

timeDelayVecIdealMSDPP = (idealCoeffs * idealDelays).';

%%% Range Constraint
rangeConsDelays          = idealDelays(abs(idealDelays) <= abs(ttdTimeDelayRange));
nRangeConsStages         = size(rangeConsDelays,1);
rangeConsScalingPerStage = 2.^(0:nRangeConsStages-1);

for antIdx = 0:nTx-1
  ttdNumber = antIdx;
  for stageIdx = nRangeConsStages:-1:1
    rangeConsCoeffs(antIdx+1, stageIdx) = floor(ttdNumber / rangeConsScalingPerStage(stageIdx));
    ttdNumber = ttdNumber - rangeConsCoeffs(antIdx+1, stageIdx) * rangeConsScalingPerStage(stageIdx);
  end
end

timeDelayVecRangeConsMSDPP = (rangeConsCoeffs * rangeConsDelays).';

maxCascading                = max(sum(rangeConsCoeffs,2));

%%% Resolution Constraint
timeDelayVecResConsMSDPP = (idealCoeffs * (round(idealDelays./ttdTimeDelayRes)*ttdTimeDelayRes)).';

%%% Range and Resolution Constraint
timeDelayVecRangeResConsMSDPP = (rangeConsCoeffs * (round(rangeConsDelays./ttdTimeDelayRes)*ttdTimeDelayRes)).';

%%% APS
txPrecoderAPS               = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(beamDirection)*(0:nTx-1).');
%%% Parallel
txPrecoderIdealPDPP         = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(fixedAngle)*(0:nTx-1).') .* ...
  exp(-1j*2*pi* (centralFreqRep-fCarrAll) .* timeDelayVecIdealPDPP.');
txPrecoderRangeConsPDPP     = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(fixedAngle)*(0:nTx-1).') .* ...
  exp(-1j*2*pi* (centralFreqRep-fCarrAll) .* timeDelayVecRangeConsPDPP.');
txPrecoderResConsPDPP       = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(fixedAngle)*(0:nTx-1).') .* ...
  exp(-1j*2*pi* (centralFreqRep-fCarrAll) .* timeDelayVecResConsPDPP.');
txPrecoderRangeResConsPDPP  = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(fixedAngle)*(0:nTx-1).') .* ...
  exp(-1j*2*pi* (centralFreqRep-fCarrAll) .* timeDelayVecRangeResConsPDPP.');
%%% Serial
txPrecoderIdealCDPP         = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(fixedAngle)*(0:nTx-1).') .* ...
  exp(-1j*2*pi* (centralFreqRep-fCarrAll) .* timeDelayVecIdealCDPP.');
txPrecoderRangeConsCDPP     = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(fixedAngle)*(0:nTx-1).') .* ...
  exp(-1j*2*pi* (centralFreqRep-fCarrAll) .* timeDelayVecRangeConsCDPP.');
txPrecoderResConsCDPP       = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(fixedAngle)*(0:nTx-1).') .* ...
  exp(-1j*2*pi* (centralFreqRep-fCarrAll) .* timeDelayVecResConsCDPP.');
txPrecoderRangeResConsCDPP  = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(fixedAngle)*(0:nTx-1).') .* ...
  exp(-1j*2*pi* (centralFreqRep-fCarrAll) .* timeDelayVecRangeResConsCDPP.');
%%% MSDPP
txPrecoderIdealMSDPP        = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(fixedAngle)*(0:nTx-1).') .* ...
  exp(-1j*2*pi* (centralFreqRep-fCarrAll) .* timeDelayVecIdealMSDPP.');
txPrecoderRangeConsMSDPP    = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(fixedAngle)*(0:nTx-1).') .* ...
  exp(-1j*2*pi* (centralFreqRep-fCarrAll) .* timeDelayVecRangeConsMSDPP.');
txPrecoderResConsMSDPP      = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(fixedAngle)*(0:nTx-1).') .* ...
  exp(-1j*2*pi* (centralFreqRep-fCarrAll) .* timeDelayVecResConsMSDPP.');
txPrecoderRangeResMSDPP     = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(fixedAngle)*(0:nTx-1).') .* ...
  exp(-1j*2*pi* (centralFreqRep-fCarrAll) .* timeDelayVecRangeResConsMSDPP.');


%%% Channel and Array Gain
parfor chanIdx = 1:size(chanDirectionVec,2)
  chanDirection = chanDirectionVec(chanIdx);

  chanArrayRespFD = 1/sqrt(nTx) * exp(-1j*2*pi * fCarrAll/centralFreq *           antElemSpacing/waveLength * sind(chanDirection) .* (0:nTx-1).');
  chanArrayRespFI = 1/sqrt(nTx) * exp(-1j*2*pi * centralFreqRep./centralFreqRep * antElemSpacing/waveLength * sind(chanDirection) .* (0:nTx-1).');

  %%% Ideal
  arrayGainIdeal(chanIdx,:)               = diag(chanArrayRespFI' * repmat(txPrecoderAPS,1,nCarr));
  %%% APS
  arrayGainAPS(chanIdx,:)                 = diag(chanArrayRespFD' * repmat(txPrecoderAPS,1,nCarr));
  %%% Parallel
  arrayGainIdealPDPP(chanIdx,:)           = diag(chanArrayRespFD' * txPrecoderIdealPDPP);
  arrayGainRangeConsPDPP(chanIdx,:)       = diag(chanArrayRespFD' * txPrecoderRangeConsPDPP);
  arrayGainResConsPDPP(chanIdx,:)         = diag(chanArrayRespFD' * txPrecoderResConsPDPP);
  arrayGainRangeResConsPDPP(chanIdx,:)    = diag(chanArrayRespFD' * txPrecoderRangeResConsPDPP);
  %%% Serial
  arrayGainIdealCDPP(chanIdx,:)           = diag(chanArrayRespFD' * txPrecoderIdealCDPP);
  arrayGainRangeConsCDPP(chanIdx,:)       = diag(chanArrayRespFD' * txPrecoderRangeConsCDPP);
  arrayGainResConsCDPP(chanIdx,:)         = diag(chanArrayRespFD' * txPrecoderResConsCDPP);
  arrayGainRangeResConsCDPP(chanIdx,:)    = diag(chanArrayRespFD' * txPrecoderRangeResConsCDPP);
  %%% MSDPP
  arrayGainIdealMSDPP(chanIdx,:)          = diag(chanArrayRespFD' * txPrecoderIdealMSDPP);
  arrayGainRangeConsMSDPP(chanIdx,:)      = diag(chanArrayRespFD' * txPrecoderRangeConsMSDPP);
  arrayGainResConsMSDPP(chanIdx,:)        = diag(chanArrayRespFD' * txPrecoderResConsMSDPP);
  arrayGainRangeResConsMSDPP(chanIdx,:)   = diag(chanArrayRespFD' * txPrecoderRangeResMSDPP);
end

fontSize = 7;
yLimit = [-15 0.5];
xLimit = [0 0.9];
xTicks = [0:0.1:0.9];
yTicks = [min(yLimit):3:max(yLimit)];


fontSize = 7;
IEEE_FIG(fontSize, 4, 19)

T= tiledlayout(1,3,'TileSpacing','tight','Padding','compact');
ylabel(T, 'Beamforming gain (g(\theta, f_{k}))',  'fontsize', fontSize, 'FontName', 'times');


t0 = tiledlayout(T,1,1);
t0.Layout.Tile = 1;
t0.Layout.TileSpan = [1 1];
ax0 = axes(t0); 

b = plot(ax0, sind(chanDirectionVec), 10*log10(abs(arrayGainRangeResConsPDPP(:,[1, 16: 16:end]))),  LineStyle="-", LineWidth=0.5);  hold on;
a = plot(ax0, sind(chanDirectionVec), 10*log10(abs(arrayGainIdealPDPP(:,[1, 16:16:end]))), Color=[0 0 0], LineStyle=":", LineWidth=0.5);  hold on;

ax0.YLim  = yLimit;
ax0.XLim  = xLimit;
ax0.YTick = yTicks;
ax0.XTick = xTicks;
grid on;
xlabel(ax0, 'Beam angle (sin \theta)','FontSize',fontSize);


t1 = tiledlayout(T,1,1);
t1.Layout.Tile = 2;
t1.Layout.TileSpan = [1 1];
ax1 = axes(t1);

b = plot(ax1, sind(chanDirectionVec), 10*log10(abs(arrayGainRangeResConsCDPP(:,[1, 16:16:end]))),  LineStyle="-", LineWidth=0.5);  hold on;
a = plot(ax1, sind(chanDirectionVec), 10*log10(abs(arrayGainIdealPDPP(:,[1, 16:16:end]))), Color=[0 0 0], LineStyle=":", LineWidth=0.5);  hold on;

ax1.YLim  = yLimit;
ax1.XLim  = xLimit;
ax1.YTick = yTicks;
% ax1.YTickLabel = [];
ax1.XTick = xTicks;
grid on;
xlabel(ax1, 'Beam angle (sin \theta)','FontSize',fontSize);

t2 = tiledlayout(T,1,1);
t2.Layout.Tile = 3;
t2.Layout.TileSpan = [1 1];
ax2 = axes(t2);

b = plot(ax2, sind(chanDirectionVec), 10*log10(abs(arrayGainRangeResConsMSDPP(:,[1, 16:16:end]))),  LineStyle="-", LineWidth=0.5);  hold on;
a = plot(ax2, sind(chanDirectionVec), 10*log10(abs(arrayGainIdealPDPP(:,[1, 16:16:end]))), Color=[0 0 0], LineStyle=":", LineWidth=0.5);  hold on;

ax2.YLim  = yLimit;
ax2.XLim  = xLimit;
ax2.YTick = yTicks;
ax2.XTick = xTicks;
grid on;
xlabel(ax2, 'Beam angle (sin \theta)','FontSize',fontSize);




aa =1
