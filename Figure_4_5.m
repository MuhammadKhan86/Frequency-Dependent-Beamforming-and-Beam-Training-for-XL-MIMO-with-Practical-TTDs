function FDB_SIM_ArryGain
close all; clear all; clc;
if 1
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
ttdTimeDelayRange = 15.0e-9;
ttdTimeDelayRes   = 3e-12;

%%% --- NEW PHYSICS PARAMETERS ---
% 1. Thermal Noise Calculation (PSD based)
k_boltz     = physconst('Boltzmann');
T_sys       = 290; % Kelvin
NF_dB       = 6;   % Noise Figure
NoisePower_Watts = k_boltz * T_sys * sampFreq * 10^(NF_dB/10); % P_noise = kTB * F

% 2. Path Loss (Fixed Distance)
targetDistance = 50; % meters
PathLoss_Lin   = (waveLength / (4 * pi * targetDistance))^2; % FSPL

% 3. Circuit Power for EE Calculation
P_circuit_Watts = 35; % Assumed static power consumption (e.g., LNA, ADC, TTDs)
%%% -----------------------------

%%% for random nmbers
% rng(1)
nSamples = 1000;
chanDirectionVec = (minBeamAngle + (maxBeamAngle-minBeamAngle).*rand(nSamples,1)).';
chanDirectionVec = -60:0.05:60;


beamDirection    = 30;


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
  
  chanArrayRespFD = GenerateChannel(nTx, nCarr, fCarrAll, centralFreq, antElemSpacing, waveLength, chanDirection);
  chanArrayRespFI = GenerateChannel(nTx, nCarr, centralFreqRep, centralFreq, antElemSpacing, waveLength, chanDirection);
  %%% Ideal
  arrayGainIdeal(chanIdx,:)               = diag(chanArrayRespFI' * repmat(txPrecoderAPS,1,nCarr));
  %%% APS
  arrayGainAPS(chanIdx,:)                 = diag(chanArrayRespFD' * repmat(txPrecoderAPS,1,nCarr));
  %%% PerfectCSI
  perfectGainIdealCSI(chanIdx,:)          = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * chanArrayRespFD(:, optBeamIdx(chanIdx))); 
  %%% PDPP
  arrayGainIdealPDPP(chanIdx,:)           = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * txPrecoderIdealPDPP(:, optBeamIdx(chanIdx)));
  arrayGainRangeConsPDPP(chanIdx,:)       = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * txPrecoderRangeConsPDPP(:, optBeamIdx(chanIdx)));
  arrayGainResConsPDPP(chanIdx,:)         = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * txPrecoderResConsPDPP(:, optBeamIdx(chanIdx)));
  arrayGainRangeResConsPDPP(chanIdx,:)    = (chanArrayRespFD(:, optBeamIdx(chanIdx))' * txPrecoderRangeResConsPDPP(:, optBeamIdx(chanIdx)));
  %%% CDPP
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
IEEE_FIG(fontSize, 8.9, 4)
figure(1);
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
c.Color     = RGB('MatlabDarkGreen');
c.Marker    = "hexagram";
c.MarkerSize=4;
c.MarkerFaceColor = RGB('MatlabDarkGreen');
c.MarkerEdgeColor = RGB('MatlabDarkGreen');
grid on;
ax.XLim = [0 9];
ax.YLim = [0 1];
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
a1 = plot(x1,y1);hold on;
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
a2.Color     = RGB('MatlabDarkGreen');
a11 = plot(markerpoints2,y22); hold on;
a11.LineStyle = "none";
a11.Marker    = "hexagram";
a11.MarkerSize=4;
a11.MarkerFaceColor = RGB('MatlabDarkGreen');
a11.MarkerEdgeColor = RGB('MatlabDarkGreen');



%%% --- SE vs TX POWER PLOT ---
% Define X-axis: Transmit Power in dBm (e.g. -10 to 50 dBm)
txPowerVecdBm = 0:0.1:40;
txPowerVecLin = 10.^(txPowerVecdBm/10);
% Pre-calculate Mean Channel Gain (Normalized) from the simulation results
% We use the mean of the absolute squared values from the loop
meanGainPerfect = mean(abs(perfectGainIdealCSI), 'all');
meanGainIdeal   = mean(abs(arrayGainIdealPDPP), 'all');
meanGainPDPP    = mean(abs(arrayGainRangeResConsPDPP), 'all');
meanGainCDPP    = mean(abs(arrayGainRangeResConsCDPP), 'all');
meanGainMSDPP   = mean(abs(arrayGainRangeResConsMSDPP), 'all');
% Calculate SE for each Tx Power step
for idx = 1:length(txPowerVecLin)
  P_tx = txPowerVecLin(idx);
  
  % SNR = (TxPower * PathLoss * ChannelGain) / NoisePower
  SNR_Perfect = (P_tx * PathLoss_Lin * meanGainPerfect) / NoisePower_Watts;
  SNR_Ideal   = (P_tx * PathLoss_Lin * meanGainIdeal)   / NoisePower_Watts;
  SNR_PDPP    = (P_tx * PathLoss_Lin * meanGainPDPP)    / NoisePower_Watts;
  SNR_CDPP    = (P_tx * PathLoss_Lin * meanGainCDPP)    / NoisePower_Watts;
  SNR_MSDPP   = (P_tx * PathLoss_Lin * meanGainMSDPP)   / NoisePower_Watts;
  se_PerfectCSI(idx) = log2(1 + SNR_Perfect);
  se_ideal(idx)      = log2(1 + SNR_Ideal);
  se_PDPP(idx)       = log2(1 + SNR_PDPP);
  se_CDPP(idx)       = log2(1 + SNR_CDPP);
  se_MSDPP(idx)      = log2(1 + SNR_MSDPP);
  
  %%% ENERGY EFFICIENCY (EE vs Tx Power)
  % EE = (SE ) / TotalPower
  % TotalPower = TxPower + CircuitPower
  P_total = P_tx + P_circuit_Watts;
  ee_tx_Perfect(idx) = (se_PerfectCSI(idx) * sampFreq) / P_total;
  ee_tx_PDPP(idx)    = (se_PDPP(idx) * sampFreq) / P_total;
  ee_tx_CDPP(idx)    = (se_CDPP(idx) * sampFreq) / P_total;
  ee_tx_MSDPP(idx)   = (se_MSDPP(idx) * sampFreq) / P_total;


    ee_tx_Perfect(idx) = (se_PerfectCSI(idx) ) / P_total;
  ee_tx_PDPP(idx)    = (se_PDPP(idx) ) / P_total;
  ee_tx_CDPP(idx)    = (se_CDPP(idx) ) / P_total;
  ee_tx_MSDPP(idx)   = (se_MSDPP(idx) ) / P_total;
end
maxSe = max(se_PerfectCSI);

figure(2);
ax = axes;
a  = plot(ax, nan, nan); hold on;
a.LineWidth = 1;
a.LineStyle = "-";
a.Color     = RGB('Black');
a.Marker    = "none";
a.MarkerSize=4;
a.MarkerFaceColor = RGB('Black');
a.MarkerEdgeColor = RGB('Black');
a1 = plot(ax, txPowerVecdBm, se_PDPP./maxSe); hold on;
a1.LineWidth = 1;
a1.LineStyle = "-";
a1.Color     = RGB('MatlabBlue');
a1.Marker    = "none";
a1.MarkerSize=4;
a1.MarkerFaceColor = RGB('MatlabBlue');
a1.MarkerEdgeColor = RGB('MatlabBlue');
% a1.MarkerIndices = [1:6:length(snrVecDb)];
a2 = plot(ax, txPowerVecdBm, se_CDPP./maxSe); hold on;
a2.LineWidth = 1;
a2.LineStyle = "-";
a2.Color     = RGB('MatlabRed');
a2.Marker    = "none";
a2.MarkerSize=4;
a2.MarkerFaceColor = RGB('MatlabRed');
a2.MarkerEdgeColor = RGB('MatlabRed');
% a2.MarkerIndices = [1:6:length(snrVecDb)];
a3 = plot(ax, txPowerVecdBm, se_MSDPP./maxSe); hold on;
a3.LineWidth = 1;
a3.LineStyle = "-";
a3.Color     = RGB('MatlabDarkGreen');
a3.Marker    = "none";
a3.MarkerSize=4;
a3.MarkerFaceColor = RGB('MatlabDarkGreen');
a3.MarkerEdgeColor = RGB('MatlabDarkGreen');
% a3.MarkerIndices = [1:6:length(snrVecDb)];
a0 = plot(ax, txPowerVecdBm, se_PerfectCSI./maxSe); hold on; grid on;
a0.LineWidth = 1;
a0.LineStyle = "-";
a0.Color     = RGB('Black');
a0.Marker    = "none";
a0.MarkerSize=4;
a0.MarkerFaceColor = RGB('Black');
a0.MarkerEdgeColor = RGB('Black');
% a0.MarkerIndices = [1:6:length(snrVecDb)];
ax.XLim   = [min(txPowerVecdBm) max(txPowerVecdBm)];
ax.YLim   = [0 1.02];
% ax.XTick  = [min(snrVecDb):3:max(snrVecDb)];
ax.YTick  = [0:0.1:1];
legend('Perfect CSI', 'Parallel [8]','Serial [7]', 'MSDPP', '',  'fontsize', fontSize)
n0 = legend;
set(n0,'NumColumns',1,'FontSize',fontSize);
n0.ItemTokenSize = 4*[5 , 5, 5 ,5 , 5, 5, 5];
n0.Position = [0.17    0.5553    0.2044    0.3278];
xlabel(ax, 'Transmit Power [dBm]', 'FontSize',fontSize, 'FontName', 'times');
ylabel(ax, 'Normalized SE [bps/Hz]',   'FontSize', fontSize, 'FontName', 'times')

%%% --- FIGURE 3: EE vs TX POWER ---
figure(3);
ax4 = axes;
maxEe = max(ee_tx_Perfect);


plot(ax4, txPowerVecdBm, ee_tx_Perfect/maxEe, '-', 'Color', RGB('black')', 'LineWidth', 1, 'Marker', 'none', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on;
plot(ax4, txPowerVecdBm, ee_tx_PDPP/maxEe,    '-', 'Color', RGB('MatlabBlue'), 'LineWidth', 1, 'MarkerSize', 4, 'MarkerFaceColor', RGB('MatlabBlue'), 'MarkerEdgeColor', RGB('MatlabBlue'));
plot(ax4, txPowerVecdBm, ee_tx_CDPP/maxEe,    '-', 'Color', RGB('MatlabRed'), 'LineWidth', 1, 'MarkerSize', 4, 'MarkerFaceColor', RGB('MatlabRed'), 'MarkerEdgeColor', RGB('MatlabRed'));
plot(ax4, txPowerVecdBm, ee_tx_MSDPP/maxEe,   '-', 'Color', RGB('MatlabDarkGreen'), 'LineWidth', 1, 'MarkerSize', 4, 'MarkerFaceColor', RGB('MatlabDarkGreen'), 'MarkerEdgeColor', RGB('MatlabDarkGreen'));
grid on;
xlabel(ax4, 'Transmit Power [dBm]', 'FontSize', fontSize, 'FontName', 'times');
ylabel(ax4, 'Energy Efficiency [bits//Joule]', 'FontSize', fontSize, 'FontName', 'times');
legend('Perfect CSI', 'Parallel [8]', 'Serial [7]', 'MSDPP', 'Location', 'northeast', 'FontSize', fontSize);

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
function c = RGB(colorName)
    switch colorName
        case 'MatlabBlue', c = [0 0.4470 0.7410];
        case 'MatlabRed', c = [0.8500 0.3250 0.0980];
        case 'MatlabYellow', c = [0.9290 0.6940 0.1250];
        case 'MatlabPurple', c = [0.4940 0.1840 0.5560];
        case 'MatlabGreen', c = [0.4660 0.6740 0.1880];
        case 'MatlabDarkGreen', c = [0.4660 0.6740 0.1880]*0.7;
        case 'MatlabLightBlue', c = [0.3010 0.7450 0.9330];
        case 'Black', c = [0 0 0];
        otherwise, c = [0 0 0];
    end
end
function IEEE_FIG(fontSize, width, height)
    set(0,'DefaultAxesFontName','Times New Roman');
    set(0,'DefaultTextFontName','Times New Roman');
    set(0,'DefaultAxesFontSize',fontSize);
    set(0,'DefaultTextFontSize',fontSize);
    f = figure(1);
    f.Units = 'centimeters';
    f.Position = [10 10 width height];
    f = figure(2);
    f.Units = 'centimeters';
    f.Position = [20 10 width height];
    f = figure(3);
    f.Units = 'centimeters';
    f.Position = [30 10 width height];
end