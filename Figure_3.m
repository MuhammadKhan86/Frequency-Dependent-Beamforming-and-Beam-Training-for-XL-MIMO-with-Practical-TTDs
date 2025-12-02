function FDB_Beam_Pattern

close all; clear all; clc;

% parpool('HPCServer',50);
%%% parameters
params.centralFreq      = 100e9;                                              centralFreq    = params.centralFreq;
params.sampFreq         = 10e9;                                            sampFreq       = params.sampFreq;
params.nTx              = 512;                                               nTx            = params.nTx;
params.nTTD             = nTx;                                              nTTD           = params.nTTD;
params.waveLength       = physconst('LightSpeed')/centralFreq ;             waveLength     = params.waveLength;
params.antElemSpacing   = 1/2 * waveLength;                                 antElemSpacing = params.antElemSpacing;
params.maxBeamAngle     =  60;                                              maxBeamAngle   = params.maxBeamAngle;
params.minBeamAngle     = -60;                                              minBeamAngle   = params.minBeamAngle;

params.antArraySize     = nTx*antElemSpacing;

%%% TTD params
ttdTimeDelayRange = 15e-9;
ttdTimeDelayRes   = 3e-12;

% ttdTimeDelayRange = 1.28e-9;
% ttdTimeDelayRes   = 0.5e-12;

beamDirection    = 30;
chanDirectionVec = 0:0.01:65;



%%% Subcarrier Frequencies
% nCarr           = ceil((1/(sqrt(3)*2*pi/nTx)*2*pi*(1 - (sampFreq.^2/(4*centralFreq.^2)) ) * (sind(maxBeamAngle) - sind(minBeamAngle))) + 1);
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

  chanArrayRespFD = GenerateChannel(nTx, nCarr, fCarrAll, centralFreq, antElemSpacing, waveLength, chanDirection);
  chanArrayRespFI = GenerateChannel(nTx, nCarr, centralFreqRep, centralFreq, antElemSpacing, waveLength, chanDirection);


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

fontSize = 7;
yLimit = [-15 0];
xLimit = [0.7 0.9];
xTicks = [0:0.1:0.9];
yTicks = [min(yLimit):3:max(yLimit)];

if 1
  IEEE_FIG(fontSize, 6, 9)

  C = colororder([...
    RGB('MatlabBlue'); RGB('MatlabRed'); RGB('MatlabYellow'); RGB('MatlabPurple'); RGB('MatlabLightGreen');...
    RGB('MatlabLightBlue'); RGB('MatlabDarkRed'); RGB('MatlabDarkGreen'); RGB('HotPink'); RGB('OrangeRed');...
    RGB('Gold'); RGB('Tan'); RGB('Cyan'); RGB('Violet'); RGB('DarkSlateBlue'); RGB('BlueViolet'); RGB('Crimson');...
    RGB('IndianRed'); RGB('Salmon')]);

  figure(1);
  ax0 =axes;
  % ax0.colororder(C(1:end,:))

  a = plot(ax0, sind(chanDirectionVec(4001:end)), 10*log10(abs(arrayGainIdealPDPP(4001:end,[321:16:end]))), Color=[0 0 0], LineStyle=":", LineWidth=0.5);  hold on;
  ax0.ColorOrder = C;
  ax0.ColorOrderIndex = 1;
  b = plot(ax0, sind(chanDirectionVec(4001:end)), 10*log10(abs(arrayGainRangeResConsPDPP(4001:end,[321: 16:end]))),  LineStyle="-", LineWidth=0.5);  hold on;
  ax0.YLim  = yLimit;
  ax0.XLim  = xLimit;
  ax0.YTick = yTicks;
  ax0.XTick = xTicks;
  grid on;
  xlabel(ax0, 'Beam angle (sin \theta)','FontSize',fontSize);
  ylabel(ax0, 'Beamforming gain (g(\theta, f_{k}))',  'fontsize', fontSize, 'FontName', 'times')

  ax=gcf;
  % exportgraphics(ax,'Results\BaemPattern_Parallel.pdf','Resolution',1600)
  % exportgraphics(ax,'Results\BaemPattern_Parallel.png','Resolution',1600)
  % exportgraphics(ax,'Results\BaemPattern_Parallel.eps','Resolution',1600)
  % axis off
  % set(gca, 'Position', [0 0 1 0.999])
  % exportgraphics(gca,'Results\Beam_Pattern_Parallel.pdf','Resolution',2400, 'ContentType', 'vector', 'BackgroundColor', 'none')
  % exportgraphics(gca,'Results\Beam_Pattern_Parallel.eps','Resolution',2400, 'ContentType', 'vector', 'BackgroundColor', 'none')

  figure(2);
  ax1 =axes;
  a = plot(sind(chanDirectionVec(4001:end)), 10*log10(abs(arrayGainIdealPDPP(4001:end,[321:16:end]))), Color=[0 0 0], LineStyle=":", LineWidth=0.5);  hold on;
  ax1.ColorOrder = C;
  ax1.ColorOrderIndex = 1;
  b = plot(sind(chanDirectionVec(4001:end)), 10*log10(abs(arrayGainRangeResConsCDPP(4001:end,[321:16:end]))),  LineStyle="-", LineWidth=0.5);  hold on;

  ax1.YLim  = yLimit;
  ax1.XLim  = xLimit;
  ax1.YTick = yTicks;
  ax1.XTick = xTicks;
  grid on;
  xlabel(ax1, 'Beam angle (sin \theta)','FontSize',fontSize);
  ylabel(ax1, 'Beamforming gain (g(\theta, f_{k}))',  'fontsize', fontSize, 'FontName', 'times')

  ax=gcf;
  % exportgraphics(ax,'Results\BaemPattern_Series.pdf','Resolution',1600)
  % exportgraphics(ax,'Results\BaemPattern_Series.png','Resolution',1600)
  % exportgraphics(ax,'Results\BaemPattern_Series.eps','Resolution',1600)
  axis off
  % set(gca, 'Position', [0 0 1 0.999])
  % exportgraphics(gca,'Results\Beam_Pattern_Serial.pdf','Resolution',2400, 'ContentType', 'vector', 'BackgroundColor', 'none')
  % exportgraphics(gca,'Results\Beam_Pattern_Serial.eps','Resolution',2400, 'ContentType', 'vector', 'BackgroundColor', 'none')


  figure(3);
  ax2 =axes;
  a = plot(sind(chanDirectionVec(4001:end)), 10*log10(abs(arrayGainIdealPDPP(4001:end,[321:16:end]))), Color=[0 0 0], LineStyle=":", LineWidth=0.5);  hold on;
  ax2.ColorOrder = C;
  ax2.ColorOrderIndex = 1;
  b = plot(sind(chanDirectionVec(4001:end)), 10*log10(abs(arrayGainRangeResConsMSDPP(4001:end,[321:16:end]))),  LineStyle="-", LineWidth=0.5);  hold on;

  ax2.YLim  = yLimit;
  ax2.XLim  = xLimit;
  ax2.YTick = yTicks;
  ax2.XTick = xTicks;
  grid on;
  xlabel(ax2, 'Beam angle (sin \theta)','FontSize',fontSize);
  ylabel(ax2, 'Beamforming gain (g(\theta, f_{k}))',  'fontsize', fontSize, 'FontName', 'times')
end
end

function chanArrayResp = GenerateChannel(nTx, nCarr, fCarrAll, centralFreq, antElemSpacing, waveLength,chanDirection)
nPaths = 7;
rng(1);
%%% 1. Define Multipath Parameters
% Angles of arrival for each path (Degrees)
pathAngles = [chanDirection randi([-60,60],1,nPaths-1)];
% Complex Gains (Magnitude * Phase)
pathGains  = [1.0, randi([0,10],1,nPaths-1)*1/100];
pathGains  = [1.0, randi([0,30],1,1)*1/1000*exp(1j*pi/3), ...
  randi([0,30],1,1)*1/1000*exp(1j*pi/4), randi([0,30],1,1)*1/1000*exp(-1j*pi/5),...
  randi([0,30],1,1)*1/1000*exp(1j*pi/2), randi([0,30],1,1)*1/1000*exp(-1j*pi/2),...
  randi([0,30],1,1)*1/1000*exp(-1j*pi/9)];
% Path Delays (Seconds) - creates frequency selectivity
pathDelays = [0, randi([0,100],1,nPaths-1)*1e-9];

% pathDelays = [0, randi([0,0],1,nPaths-1)*1e-9];

% Ensure dimensions match
nPaths = length(pathAngles);

%%% 2. Construct Multipath Channel Response
% Initialize channel matrix: [nTx x nCarr]
chanArrayResp = zeros(nTx, nCarr);

for p = 1:nPaths
  % 1. Compute Spatial Steering Vector (Frequency Dependent Squint)
  % Size: [nTx x nCarr]
  steerVec = 1/sqrt(nTx) * exp(-1j*2*pi * (fCarrAll/centralFreq) .* ...
    (antElemSpacing/waveLength) * sind(pathAngles(p)) .* (0:nTx-1).');

  % 2. Compute Delay Term (Frequency Selective Fading)
  % Size: [1 x nCarr]
  delayTerm = exp(-1j*2*pi * fCarrAll * pathDelays(p));

  % 3. Accumulate weighted path into total channel
  % We use implicit expansion or repmat to multiply [nTx x nCarr] by [1 x nCarr]
  chanArrayResp = chanArrayResp + ...
    (pathGains(p) * steerVec .* delayTerm);
end
end



