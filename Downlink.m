%%
clc
clear
close all
%#ok<*UNRCH>

%%
Num_Figure=0;
DEBUG=0;
DECIMATION=0;

%%
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
DL_Config=lteRMCDL(DL_RMC,Duplex_Mode,Num_SubFrame);
DL_Num_Rx_Bit=sum(DL_Config.PDSCH.TrBlkSizes);
DL_Tx_Data=randi([0,1],DL_Num_Rx_Bit,1);
[DL_Tx_Waveform,DL_TX_Grid,DL_Config]=lteRMCDLTool(DL_Config,DL_Tx_Data);
if DEBUG
    DL_Tx_Spectrum_Scope=Spectrum_Scope.clone();
    DL_Tx_Spectrum_Scope.SampleRate=DL_Config.SamplingRate;
    DL_Tx_Spectrum_Scope.step(DL_Tx_Waveform);
    DL_Tx_Spectrum_Scope.Title='Downlink Tx Spectrum';
    DL_Tx_Constellation_Diagram=Constellation_Diagram.clone();
    DL_Tx_Constellation_Diagram.Title='Downlink Tx Constellation_Diagramm';
    DL_Tx_Constellation_Diagram.step(DL_TX_Grid(:))
end

%%
if DECIMATION
    DL_Tx_Waveform=resample(DL_Tx_Waveform,2,3);
    if DEBUG
        DL_Tx_Spectrum_Scope_Dec=Spectrum_Scope.clone();
        DL_Tx_Spectrum_Scope_Dec.SampleRate=DL_Config.SamplingRate/3*2;
        DL_Tx_Spectrum_Scope_Dec.step(DL_Tx_Waveform);
        DL_Tx_Spectrum_Scope_Dec.Title='Downlink Tx Spectrum After resampling @ 10.24MHz';
    end
end

%%
if DECIMATION
    Chan_Config.SamplingRate=DL_Config.SamplingRate/3*2;
else
    Chan_Config.SamplingRate=DL_Config.SamplingRate;
end
DL_Rx_Waveform=lteFadingChannel(Chan_Config,[DL_Tx_Waveform;zeros(Num_Delay_Sample,size(DL_Tx_Waveform,2))]);
if DEBUG
    DL_Rx_Spectrum_Scope=Spectrum_Scope.clone();
    DL_Rx_Spectrum_Scope.SampleRate=Chan_Config.SamplingRate;
    DL_Rx_Spectrum_Scope.step(DL_Rx_Waveform);
    DL_Rx_Spectrum_Scope.Title='Downlink Rx Spectrum';
end

%%
if DECIMATION
    DL_Rx_Waveform=resample(DL_Rx_Waveform,3,2);
    if DEBUG
        DL_Tx_Spectrum_Scope_Interp=Spectrum_Scope.clone();
        DL_Tx_Spectrum_Scope_Interp.SampleRate=DL_Config.SamplingRate;
        DL_Tx_Spectrum_Scope_Interp.step(DL_Rx_Waveform);
        DL_Tx_Spectrum_Scope_Interp.Title='Downlink Rx Spectrum After resampling @ 15.36MHz';
    end
end

%%
DL_Frame_Offset=lteDLFrameOffset(DL_Config,DL_Rx_Waveform);
DL_Rx_Waveform(1:DL_Frame_Offset)=[];

%%
DL_Rx_Grid=lteOFDMDemodulate(DL_Config,DL_Rx_Waveform);
if DEBUG
    DL_Rx_Constellation_Diagram=Constellation_Diagram.clone();
    DL_Rx_Constellation_Diagram.Title='Downlink Rx Constellation_Diagramm Before Equalization';
    DL_Rx_Constellation_Diagram.step(DL_Rx_Grid(:))
end

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
    DL_Eq_Constellation_Diagram.Title='Downlink Rx Constellation_Diagramm After Equalization';
    DL_Eq_Constellation_Diagram.step(DL_Eq_Grid(:))
end

%%
DL_EVM=lteEVM(DL_TX_Grid,DL_Eq_Grid);
if DEBUG
    [PDSCH_EVM,PDSCH_EVM_Plot]=hPDSCHEVM(DL_Config,struct('PilotAverage','TestEVM'),DL_Rx_Waveform);
end
disp(DL_EVM.RMS*100)