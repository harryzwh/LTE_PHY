%%
clc
clear
close all
%#ok<*UNRCH>

%%
Num_Figure=0;
DEBUG=0;

%%
UL_RMC='A5-5';
DL_RMC='R.7';
Duplex_Mode='FDD';
Num_SubFrame=20;
Num_Delay_Sample=25;
Chan_Config.DelayProfile = 'EPA';
Chan_Config.NRxAnts = 1;
Chan_Config.DopplerFreq = 5;
Chan_Config.MIMOCorrelation = 'Low';
Chan_Config.Seed = 1;
Chan_Config.InitPhase = 'Random';
Chan_Config.ModelType = 'GMEDS';
Chan_Config.NTerms = 16;
Chan_Config.NormalizeTxAnts = 'On';
Chan_Config.NormalizePathGains = 'On';
Chan_Config.InitTime=0;

%%
Spectrum_Scope=dsp.SpectrumAnalyzer();
Constellation_Diagram=comm.ConstellationDiagram();
Constellation_Diagram.ShowReferenceConstellation=false;

%%
UL_Config=lteRMCUL(UL_RMC,Duplex_Mode,Num_SubFrame);
UL_Num_Tx_Bit=sum(UL_Config.PUSCH.TrBlkSizes);
UL_Tx_Data=randi([0,1],UL_Num_Tx_Bit,1);
[UL_Tx_Waveform,UL_Tx_Grid,UL_Config]=lteRMCULTool(UL_Config,UL_Tx_Data);
if DEBUG
    UL_Tx_Spectrum_Scope=Spectrum_Scope.clone();
    UL_Tx_Spectrum_Scope.SampleRate=UL_Config.SamplingRate;
    UL_Tx_Spectrum_Scope.step(UL_Tx_Waveform);
    UL_Tx_Spectrum_Scope.Title='Uplink Tx Spectrum';
    UL_Tx_Constellation_Diagram=Constellation_Diagram.clone();
    UL_Tx_Constellation_Diagram.Title='Uplink Tx Constellation_Diagramm';
    UL_Tx_Constellation_Diagram.step(UL_Tx_Grid(:))
end

%%
DL_Config=lteRMCDL(DL_RMC,Duplex_Mode,Num_SubFrame);
DL_Num_Rx_Bit=sum(DL_Config.PDSCH.TrBlkSizes);
DL_Tx_Data=randi([0,1],DL_Num_Rx_Bit,1);
[DL_Tx_Waveform,DL_Grid,DL_Config]=lteRMCDLTool(DL_Config,DL_Tx_Data);
if DEBUG
    DL_Tx_Spectrum_Scope=Spectrum_Scope.clone();
    DL_Tx_Spectrum_Scope.SampleRate=DL_Config.SamplingRate;
    DL_Tx_Spectrum_Scope.step(DL_Tx_Waveform);
    DL_Tx_Spectrum_Scope.Title='Downlink Tx Spectrum';
    DL_Tx_Constellation_Diagram=Constellation_Diagram.clone();
    DL_Tx_Constellation_Diagram.Title='Downlink Tx Constellation_Diagramm';
    DL_Tx_Constellation_Diagram.step(DL_Grid(:))
end

%%
UL_Tx_Waveform_Dec=resample(UL_Tx_Waveform,2,3);
if DEBUG
    UL_Tx_Spectrum_Scope_Dec=Spectrum_Scope.clone();
    UL_Tx_Spectrum_Scope_Dec.SampleRate=UL_Config.SamplingRate/3*2;
    UL_Tx_Spectrum_Scope_Dec.step(UL_Tx_Waveform_Dec);
    UL_Tx_Spectrum_Scope_Dec.Title='Uplink Tx Spectrum After resampling @ 10.24MHz';
end

%%
Chan_Config.SamplingRate=UL_Config.SamplingRate;
UL_Rx_Waveform=lteFadingChannel(Chan_Config,[UL_Tx_Waveform;zeros(Num_Delay_Sample,size(UL_Tx_Waveform,2))]);
if DEBUG
    UL_Rx_Spectrum_Scope=Spectrum_Scope.clone();
    UL_Rx_Spectrum_Scope.SampleRate=UL_Config.SamplingRate;
    UL_Rx_Spectrum_Scope.step(UL_Rx_Waveform);
    UL_Rx_Spectrum_Scope.Title='Uplink Rx Spectrum';
end

%%
DL_Rx_Waveform=lteFadingChannel(Chan_Config,[DL_Tx_Waveform;zeros(Num_Delay_Sample,size(DL_Tx_Waveform,2))]);
if DEBUG
    DL_Rx_Spectrum_Scope=Spectrum_Scope.clone();
    DL_Rx_Spectrum_Scope.SampleRate=DL_Config.SamplingRate;
    DL_Rx_Spectrum_Scope.step(DL_Rx_Waveform);
    DL_Rx_Spectrum_Scope.Title='Downlink Rx Spectrum';
end

%%
UL_Frame_Offset=lteULFrameOffset(UL_Config,UL_Config.PUSCH,UL_Rx_Waveform);
UL_Rx_Waveform(1:UL_Frame_Offset)=[];

%%
DL_Frame_Offset=lteDLFrameOffset(DL_Config,DL_Rx_Waveform);
DL_Rx_Waveform(1:DL_Frame_Offset)=[];

%%
UL_Rx_Grid=lteSCFDMADemodulate(UL_Config,UL_Rx_Waveform);
if DEBUG
    UL_Rx_Constellation_Diagram=Constellation_Diagram.clone();
    UL_Rx_Constellation_Diagram.Title='Uplink Tx Constellation_Diagramm Before Equalization';
    UL_Rx_Constellation_Diagram.step(UL_Rx_Grid(:))
end

%%
DL_Rx_Grid=lteOFDMDemodulate(DL_Config,DL_Rx_Waveform);
if DEBUG
    DL_Rx_Constellation_Diagram=Constellation_Diagram.clone();
    DL_Rx_Constellation_Diagram.Title='Downlink Tx Constellation_Diagramm Before Equalization';
    DL_Rx_Constellation_Diagram.step(DL_Rx_Grid(:))
end

%%
% Estimator_Config=struct('FreqWindow',7,'TimeWindow',1,'InterpType','cubic','Window','Right');
% UL_H_Est=complex(zeros(size(UL_Rx_Grid)));
% pucch1.ResourceIdx  = 0;
% pucch1.DeltaShift   = 1;
% pucch1.CyclicShifts = 0;
% for SubFrame_Idx=1:Num_SubFrame
% [UL_H_Est(:,(SubFrame_Idx-1)*14+(1:14)), UL_Noise_Est]=lteULChannelEstimate(UL_Config,UL_Config.PUSCH,Estimator_Config,UL_Rx_Grid(:,(SubFrame_Idx-1)*14+(1:14)));
% end
% H = lteULPerfectChannelEstimate(UL_Config, Chan_Config);
% mesh(abs(H).^2)
% mesh(abs(UL_H_Est).^2)

%%
DL_Estimator_Config.FreqWindow = 1;
DL_Estimator_Config.TimeWindow = 1;
DL_Estimator_Config.InterpType = 'cubic';
DL_Estimator_Config.PilotAverage = 'UserDefined';
DL_Estimator_Config.InterpWinSize = 3;
DL_Estimator_Config.InterpWindow = 'Causal';
[DL_H_Est,DL_Noise_Est]=lteDLChannelEstimate(DL_Config,DL_Estimator_Config,DL_Rx_Grid);
if DEBUG
    Num_Figure=Num_Figure+1;
    figure(Num_Figure)
    mesh(abs(DL_H_Est).^2)
    title('Downlink Channel Estimation')
end

%%
DL_Eq_Grid=lteEqualizeZF(DL_Rx_Grid,DL_H_Est);
if DEBUG
    DL_Eq_Constellation_Diagram=Constellation_Diagram.clone();
    DL_Eq_Constellation_Diagram.Title='Downlink Tx Constellation_Diagramm After Equalization';
    DL_Eq_Constellation_Diagram.step(DL_Eq_Grid(:))
end

%%
DL_EVM=lteEVM(DL_Grid,DL_Eq_Grid);
if DEBUG
    [PDSCH_EVM,PDSCH_EVM_Plot]=hPDSCHEVM(DL_Config,struct('PilotAverage','TestEVM'),DL_Rx_Waveform(1:length(DL_Tx_Waveform)));
end