function Figure_4_5

close all; clear all; clc;

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

%%% for random nmbers
% rng(1)
nSamples = 1000;
chanDirectionVec = (minBeamAngle + (maxBeamAngle-minBeamAngle).*rand(nSamples,1)).';

% nCarr            = 33;
beamDirection    = 30;
% chanDirectionVec = r;
% load('optCodebook.MAT')

%%% Subcarrier Frequencies
nCarr           = ceil((1/(sqrt(2)*1.782/nTx)* (sind(maxBeamAngle) - sind(minBeamAngle))) + 1);
fCarrAll        = (centralFreq + sampFreq / (nCarr-1) .* ((0:nCarr-1)-(nCarr - 1) / 2));
centralFreqRep  = repelem(centralFreq,1,nCarr);


%%%% Frequency Dependent Beam Generation

timeDelayIdeal  = 0.5/(centralFreq*sampFreq) * (fCarrAll(1)*sind(minBeamAngle) - fCarrAll(end)*sind(maxBeamAngle));
fixedAngle      = asind(0.5/centralFreq * ((centralFreq - sampFreq/2)*sind(minBeamAngle) + (centralFreq + sampFreq/2)*sind(maxBeamAngle)));

optCodebookAngles = asind((2*centralFreqRep.*(centralFreqRep-fCarrAll)*timeDelayIdeal + centralFreqRep*sind(fixedAngle)) ./ (fCarrAll));
for locIdx = 1:length(chanDirectionVec)
  [~, optBeamIdx(locIdx)] = min(abs(optCodebookAngles-chanDirectionVec(locIdx)));
end

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
  %%% PerfectCSI
  perfectGainIdealCSI(chanIdx,:)          = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * chanArrayRespFD(:, optBeamIdx(chanIdx))); 
  %%% Paralel
  arrayGainIdealPDPP(chanIdx,:)           = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * txPrecoderIdealPDPP(:, optBeamIdx(chanIdx)));
  arrayGainRangeConsPDPP(chanIdx,:)       = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * txPrecoderRangeConsPDPP(:, optBeamIdx(chanIdx)));
  arrayGainResConsPDPP(chanIdx,:)         = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * txPrecoderResConsPDPP(:, optBeamIdx(chanIdx)));
  arrayGainRangeResConsPDPP(chanIdx,:)    = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * txPrecoderRangeResConsPDPP(:, optBeamIdx(chanIdx)));
  %%% Serial
  arrayGainIdealCDPP(chanIdx,:)           = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * txPrecoderIdealCDPP(:, optBeamIdx(chanIdx)));
  arrayGainRangeConsCDPP(chanIdx,:)       = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * txPrecoderRangeConsCDPP(:, optBeamIdx(chanIdx)));
  arrayGainResConsCDPP(chanIdx,:)         = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * txPrecoderResConsCDPP(:, optBeamIdx(chanIdx)));
  arrayGainRangeResConsCDPP(chanIdx,:)    = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * txPrecoderRangeResConsCDPP(:, optBeamIdx(chanIdx)));
  %%% MSDPP
  arrayGainIdealMSDPP(chanIdx,:)          = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * txPrecoderIdealMSDPP(:, optBeamIdx(chanIdx)));
  arrayGainRangeConsMSDPP(chanIdx,:)      = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * txPrecoderRangeConsMSDPP(:, optBeamIdx(chanIdx)));
  arrayGainResConsMSDPP(chanIdx,:)        = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * txPrecoderResConsMSDPP(:, optBeamIdx(chanIdx)));
  arrayGainRangeResConsMSDPP(chanIdx,:)   = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * txPrecoderRangeResMSDPP(:, optBeamIdx(chanIdx)));
end


arrayGainPeerfectCSI  = 10*log10(abs(perfectGainIdealCSI));
arrayGainIdealdB      = 10*log10(abs(arrayGainIdealPDPP));
arrayGainPDPP         = 10*log10(abs(arrayGainRangeResConsPDPP));
arrayGainCDPP         = 10*log10(abs(arrayGainRangeResConsCDPP));
arrayGainMSDPP        = 10*log10(abs(arrayGainRangeResConsMSDPP));

arryGainDiffIdeal   = abs(arrayGainPeerfectCSI - arrayGainIdealdB);
arryGainDiffPDPP    = abs(arrayGainPeerfectCSI - arrayGainPDPP);
arryGainDiffCDPP    = abs(arrayGainPeerfectCSI - arrayGainCDPP);
arryGainDiffMSDPP   = abs(arrayGainPeerfectCSI - arrayGainMSDPP);

[y0, x0] = ecdf(arryGainDiffPDPP);
markerpoints0 = 0:1:max(x0);
y00 = interp1(x0(2:end),y0(2:end),markerpoints0);


[y1, x1] = ecdf(arryGainDiffCDPP);
markerpoints1 = 0:1:max(x1);
y11 = interp1(x1(2:end),y1(2:end),markerpoints1);

[y2, x2] = ecdf(arryGainDiffMSDPP);
markerpoints2 = 0:1:max(x2);
y22 = interp1(x2(2:end),y2(2:end),markerpoints2);

[y3, x3] = ecdf(arryGainDiffIdeal);
markerpoints3 = 0:1:max(x3);
y33 = interp1(x3(2:end),y3(2:end),markerpoints3);


fontSize = 7;
IEEE_FIG(fontSize, 4, 10)
figure(4);
ax = axes;

a  = plot(ax, nan, nan); hold on;
a.LineWidth = 1;
a.LineStyle = "-";
a.Color     = RGB('MatlabBlue');
a.Marker    = "o";
a.MarkerSize=4;
a.MarkerFaceColor = RGB('MatlabBlue');
a.MarkerEdgeColor = RGB('MatlabBlue');

b  = plot(ax, nan, nan); hold on;
b.LineWidth = 1;
b.LineStyle = "-";
b.Color     = RGB('MatlabRed');
b.Marker    = "square";
b.MarkerSize=4;
b.MarkerFaceColor = RGB('MatlabRed');
b.MarkerEdgeColor = RGB('MatlabRed');

c  = plot(ax, nan, nan); hold on;
c.LineWidth = 1;
c.LineStyle = "-";
c.Color     = RGB('Green');
c.Marker    = "hexagram";
c.MarkerSize=4;
c.MarkerFaceColor = RGB('Green');
c.MarkerEdgeColor = RGB('Green');

grid on;
ax.XLim = [0 9];
ax.YLim = [0 1.02];
ax.XTick = [0:1:9];
ax.YTick = [0:0.1:1];
xlabel(ax, 'Beamforming gain difference [dB]', 'FontSize',fontSize, 'FontName', 'times');
ylabel(ax, 'CDF',   'FontSize', fontSize, 'FontName', 'times')

aa  = plot(x3,y3);hold on;
aa.LineWidth = 1;
aa.LineStyle = "-";
aa.Color     = RGB('MatlabYellow');


a0  = plot(x0,y0);hold on;
a0.LineWidth = 1;
a0.LineStyle = "-";
a0.Color     = RGB('MatlabBlue');


a00 = plot(markerpoints0,y00); hold on;
a00.LineStyle = "none";
a00.Marker    = "o";
a00.MarkerSize=4;
a00.MarkerFaceColor = RGB('MatlabBlue');
a00.MarkerEdgeColor = RGB('MatlabBlue');



a1 = plot(x1,y2);hold on;
a1.LineWidth = 1;
a1.LineStyle = "-";
a1.Color     = RGB('MatlabRed');

a11 = plot(markerpoints1,y11); hold on;
a11.LineStyle = "none";
a11.Marker    = "square";
a11.MarkerSize=4;
a11.MarkerFaceColor = RGB('MatlabRed');
a11.MarkerEdgeColor = RGB('MatlabRed');



a2 = plot(x2,y2);hold on;
a2.LineWidth = 1;
a2.LineStyle = "-";
a2.Color     = RGB('Green');

a11 = plot(markerpoints2,y22); hold on;
a11.LineStyle = "none";
a11.Marker    = "hexagram";
a11.MarkerSize=4;
a11.MarkerFaceColor = RGB('Green');
a11.MarkerEdgeColor = RGB('Green');

legend('Parallel [8]','Serial [7]', 'MSDPP',  'fontsize', fontSize)
n0 = legend;
set(n0,'NumColumns',1,'FontSize',fontSize);
n0.ItemTokenSize = 4*[5 , 5, 5 ,5 , 5, 5, 5];
n0.Position = [0.68    0.63    0.2044    0.2517];



snrVecDb = [-6:0.5:21];
snrVecLin = 10.^(snrVecDb/10);

for snrIdx = 1:length(snrVecDb)
  snrLin   = snrVecLin(snrIdx);
  noisePwr = 1 / snrLin;
  se_PerfectCSI(snrIdx,:)  = log2(1 +  mean((abs(perfectGainIdealCSI)))/ noisePwr);
  se_ideal(snrIdx,:)  = log2(1 +  mean((abs(arrayGainIdealPDPP)))/ noisePwr);
  se_PDPP(snrIdx,:)   = log2(1 +  mean((abs(arrayGainRangeResConsPDPP)))  / noisePwr);
  se_CDPP(snrIdx,:)   = log2(1 +  mean((abs(arrayGainRangeResConsCDPP)))   / noisePwr);
  se_MSDPP(snrIdx,:)  = log2(1 +  mean((abs(arrayGainRangeResConsMSDPP)))   / noisePwr);
end
maxSe = max(se_PerfectCSI);


figure(5);
ax = axes;

a  = plot(ax, nan, nan); hold on;
a.LineWidth = 1;
a.LineStyle = ":";
a.Color     = RGB('Black');
a.Marker    = "*";
a.MarkerSize=4;
a.MarkerFaceColor = RGB('Black');
a.MarkerEdgeColor = RGB('Black');


a1 = plot(ax, snrVecDb, se_PDPP./maxSe); hold on;
a1.LineWidth = 1;
a1.LineStyle = "-";
a1.Color     = RGB('MatlabBlue');
a1.Marker    = "o";
a1.MarkerSize=4;
a1.MarkerFaceColor = RGB('MatlabBlue');
a1.MarkerEdgeColor = RGB('MatlabBlue');
a1.MarkerIndices = [1:6:length(snrVecDb)];


a2 = plot(ax, snrVecDb, se_CDPP./maxSe); hold on;
a2.LineWidth = 1;
a2.LineStyle = "-";
a2.Color     = RGB('MatlabRed');
a2.Marker    = "square";
a2.MarkerSize=4;
a2.MarkerFaceColor = RGB('MatlabRed');
a2.MarkerEdgeColor = RGB('MatlabRed');
a2.MarkerIndices = [1:6:length(snrVecDb)];


a3 = plot(ax, snrVecDb, se_MSDPP./maxSe); hold on;
a3.LineWidth = 1;
a3.LineStyle = "-";
a3.Color     = RGB('Green');
a3.Marker    = "hexagram";
a3.MarkerSize=4;
a3.MarkerFaceColor = RGB('Green');
a3.MarkerEdgeColor = RGB('Green');
a3.MarkerIndices = [1:6:length(snrVecDb)];

a0 = plot(ax, snrVecDb, se_PerfectCSI./maxSe); hold on; grid on;
a0.LineWidth = 1;
a0.LineStyle = ":";
a0.Color     = RGB('Black');
a0.Marker    = "*";
a0.MarkerSize=4;
a0.MarkerFaceColor = RGB('Black');
a0.MarkerEdgeColor = RGB('Black');
a0.MarkerIndices = [1:6:length(snrVecDb)];



ax.XLim   = [min(snrVecDb) max(snrVecDb)];
ax.YLim   = [0 1.02];
ax.XTick  = [min(snrVecDb):3:max(snrVecDb)];
ax.YTick  = [0:0.1:1];
legend('Perfect CSI', 'Parallel [8]','Serial [7]', 'MSDPP', '',  'fontsize', fontSize)
n0 = legend;
set(n0,'NumColumns',1,'FontSize',fontSize);
n0.ItemTokenSize = 4*[5 , 5, 5 ,5 , 5, 5, 5];
n0.Position = [0.17    0.5553    0.2044    0.3278];


xlabel(ax, 'SNR [dB]', 'FontSize',fontSize, 'FontName', 'times');
ylabel(ax, 'Normalized SE [bps/Hz]',   'FontSize', fontSize, 'FontName', 'times')



