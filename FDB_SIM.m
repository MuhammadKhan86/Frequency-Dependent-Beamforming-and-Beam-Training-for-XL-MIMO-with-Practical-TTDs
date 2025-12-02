function FDB_SIM

close all; clear all; clc;

%%% parameters
params.centralFreq      = 100e9;                                              centralFreq    = params.centralFreq;
params.sampFreq         = 10e9;                                            sampFreq       = params.sampFreq;
params.nTx              = 128;                                               nTx            = params.nTx;
params.nTTD             = nTx;                                              nTTD           = params.nTTD;
params.waveLength       = physconst('LightSpeed')/centralFreq ;             waveLength     = params.waveLength;
params.antElemSpacing   = 1/2 * waveLength;                                 antElemSpacing = params.antElemSpacing;
params.maxBeamAngle     =  60;                                              maxBeamAngle   = params.maxBeamAngle;
params.minBeamAngle     = -60;                                              minBeamAngle   = params.minBeamAngle;

params.antArraySize     = nTx*antElemSpacing;

%%% TTD params
ttdTimeDelayRange = 15e-9;
ttdTimeDelayRes   = 3e-12;

%%% for random nmbers
% rng(1)
% nSamples = 10000;
% chanDirectionVec = (minBeamAngle + (maxBeamAngle-minBeamAngle).*rand(nSamples,1)).';

% nCarr            = 33;
beamDirection    = 30;
chanDirectionVec = -60:0.01:60;
% load('optCodebook.MAT')

% for locIdx = 1:length(chanDirectionVec)
%   [~, optBeamIdx(locIdx)] = min(abs(optCodebookAngles-chanDirectionVec(locIdx)));
% end

%%% Subcarrier Frequencies
nCarr           = ceil((1/(sqrt(3)*2*pi/nTx)*2*pi*(1 - (sampFreq.^2/(4*centralFreq.^2)) ) * (sind(maxBeamAngle) - sind(minBeamAngle))) + 1);
nCarr           = ceil((1/(sqrt(3)*2*pi/nTx)*2*pi * (sind(maxBeamAngle) - sind(minBeamAngle))) + 1);

% nCarr           = ceil(120/(2*asind(0.891*0.5/nTx)))
% nCarr           = ceil(120/(3*asind(2/pi * sqrt(6*(1-0.5) / (nTx.^2-1) ))))
fCarrAll        = (centralFreq + sampFreq / (nCarr-1) .* ((0:nCarr-1)-(nCarr - 1) / 2));
centralFreqRep  = repelem(centralFreq,1,nCarr);


%%%% Frequency Dependent Beam Generation

timeDelayIdeal  = 0.5/(centralFreq*sampFreq) * (fCarrAll(1)*sind(minBeamAngle) - fCarrAll(end)*sind(maxBeamAngle));
fixedAngle      = asind(0.5/centralFreq * ((centralFreq - sampFreq/2)*sind(minBeamAngle) + (centralFreq + sampFreq/2)*sind(maxBeamAngle)));

optCodebookAngles = asind((2*centralFreqRep.*(centralFreqRep-fCarrAll)*timeDelayIdeal + centralFreqRep*sind(fixedAngle)) ./ (fCarrAll));
for locIdx = 1:length(chanDirectionVec)
  [~, optBeamIdx(locIdx)] = min(abs(optCodebookAngles-chanDirectionVec(locIdx)));
end

%%% PDPP
timeDelayVecIdealPDPP         = timeDelayIdeal * (0:nTx-1);
timeDelayVecRangeConsPDPP     = min(abs(timeDelayVecIdealPDPP), abs(ttdTimeDelayRange)) .* sign(timeDelayVecIdealPDPP);
timeDelayVecResConsPDPP       = round(timeDelayVecIdealPDPP./ttdTimeDelayRes)*ttdTimeDelayRes;
timeDelayVecRangeResConsPDPP  = round(timeDelayVecRangeConsPDPP./ttdTimeDelayRes)*ttdTimeDelayRes;

%%% CDPP
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


%%% GDPP
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
%%% PDPP
txPrecoderIdealPDPP         = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(fixedAngle)*(0:nTx-1).') .* ...
  exp(-1j*2*pi* (centralFreqRep-fCarrAll) .* timeDelayVecIdealPDPP.');
txPrecoderRangeConsPDPP     = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(fixedAngle)*(0:nTx-1).') .* ...
  exp(-1j*2*pi* (centralFreqRep-fCarrAll) .* timeDelayVecRangeConsPDPP.');
txPrecoderResConsPDPP       = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(fixedAngle)*(0:nTx-1).') .* ...
  exp(-1j*2*pi* (centralFreqRep-fCarrAll) .* timeDelayVecResConsPDPP.');
txPrecoderRangeResConsPDPP  = 1/sqrt(nTx) * exp(-1j*2*pi*antElemSpacing/waveLength * sind(fixedAngle)*(0:nTx-1).') .* ...
  exp(-1j*2*pi* (centralFreqRep-fCarrAll) .* timeDelayVecRangeResConsPDPP.');
%%% CDPP
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
  %%% PDPP
  arrayGainIdealPDPP(chanIdx,:)           = diag(chanArrayRespFD' * txPrecoderIdealPDPP);
  arrayGainRangeConsPDPP(chanIdx,:)       = diag(chanArrayRespFD' * txPrecoderRangeConsPDPP);
  arrayGainResConsPDPP(chanIdx,:)         = diag(chanArrayRespFD' * txPrecoderResConsPDPP);
  arrayGainRangeResConsPDPP(chanIdx,:)    = diag(chanArrayRespFD' * txPrecoderRangeResConsPDPP);
  %%% CDPP
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

IEEE_FIG(fontSize, 5, 6.3)
yLimit = [-21 0.5];
xLimit = [0 0.9];
xTicks = [0:0.1:0.9];
yTicks = [min(yLimit):3:max(yLimit)];
% markerIndices = [1 length(ttDTimeDelayRangeVec)/10:length(ttDTimeDelayRangeVec)/10:length(ttDTimeDelayRangeVec)];
figure(1);
ax0 =axes; 
a = plot(sind(chanDirectionVec), 20*log10(abs(arrayGainIdealPDPP(:,[1,32:32:end]))), Color=[0 0 0], LineStyle=":", LineWidth=1);  hold on;
b = plot(sind(chanDirectionVec), 20*log10(abs(arrayGainRangeResConsPDPP(:,[1,32:32:end]))),  LineStyle="-", LineWidth=1);  hold on;

ax0.YLim  = yLimit;
ax0.XLim  = xLimit;
ax0.YTick = yTicks;
ax0.XTick = xTicks;
grid on;

figure(2);
ax1 =axes; 
a = plot(sind(chanDirectionVec), 20*log10(abs(arrayGainIdealPDPP(:,[1,32:32:end]))), Color=[0 0 0], LineStyle=":", LineWidth=1);  hold on;
b = plot(sind(chanDirectionVec), 20*log10(abs(arrayGainRangeResConsCDPP(:,[1,32:32:end]))),  LineStyle="-", LineWidth=1);  hold on;

ax1.YLim  = yLimit;
ax1.XLim  = xLimit;
ax1.YTick = yTicks;
ax1.XTick = xTicks;
grid on;

figure(3);
ax2 =axes; 
a = plot(sind(chanDirectionVec), 20*log10(abs(arrayGainIdealMSDPP(:,[1,32:32:end]))), Color=[0 0 0], LineStyle=":", LineWidth=1);  hold on;
b = plot(sind(chanDirectionVec), 20*log10(abs(arrayGainRangeResConsCDPP(:,[1,32:32:end]))),  LineStyle="-", LineWidth=1);  hold on;

ax2.YLim  = yLimit;
ax2.XLim  = xLimit;
ax2.YTick = yTicks;
ax2.XTick = xTicks;
grid on;

% a.Color     = RGB('Black');


figure(1); plot(chanDirectionVec, 20*log10(abs(arrayGainAPS)));

figure(2); plot(chanDirectionVec, 20*log10(abs(arrayGainIdealPDPP)),  "LineStyle", "-");  hold on;
figure(2); plot(chanDirectionVec, 20*log10(abs(arrayGainIdealCDPP)),  "LineStyle", "--"); hold on;
figure(2); plot(chanDirectionVec, 20*log10(abs(arrayGainIdealMSDPP)), "LineStyle", ":");  hold on;

figure(3); plot(chanDirectionVec, 20*log10(abs(arrayGainRangeConsPDPP)),  "LineStyle", "-");  hold on;
figure(3); plot(chanDirectionVec, 20*log10(abs(arrayGainRangeConsCDPP)),  "LineStyle", "--"); hold on;
figure(3); plot(chanDirectionVec, 20*log10(abs(arrayGainRangeConsMSDPP)), "LineStyle", ":");  hold on;

figure(4); plot(chanDirectionVec, 20*log10(abs(arrayGainResConsPDPP)),  "LineStyle", "-");  hold on;
figure(4); plot(chanDirectionVec, 20*log10(abs(arrayGainResConsCDPP)),  "LineStyle", "--"); hold on;
figure(4); plot(chanDirectionVec, 20*log10(abs(arrayGainResConsMSDPP)), "LineStyle", ":");  hold on;

figure(5); plot(chanDirectionVec, 20*log10(abs(arrayGainIdealPDPP(:,[1,32:32:end]))),  "LineStyle", "--", Color=[0 0 0]);  hold on;
figure(5); plot(chanDirectionVec, 20*log10(abs(arrayGainRangeResConsPDPP(:,[1,32:32:end]))),  "LineStyle", "-");  hold on;
figure(5); plot(chanDirectionVec, 20*log10(abs(arrayGainRangeResConsCDPP(:,[1,32:32:end]))),  "LineStyle", "-"); hold on;
figure(5); plot(chanDirectionVec, 20*log10(abs(arrayGainRangeResConsMSDPP(:,[1,32:32:end]))), "LineStyle", ":");  hold on;


figure(10); 
a = polarpattern(chanDirectionVec, 20*log10(abs(arrayGainIdealPDPP(:,[1,32:32:end]))));hold on

b = polarpattern(chanDirectionVec, 20*log10(abs(arrayGainIdealPDPP(:,[1,32:32:end]))));
aa=1
