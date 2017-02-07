%%
clc
clear
close all
%#ok<*UNRCH>

%%
Num_Figure=0;
DEBUG=1;
DECIMATION=0;

%%
DL_RMC='R.7';
Duplex_Mode='FDD';
Num_SubFrame=10;
Num_Delay_Sample=25;
Chan_Config.DelayProfile = 'EPA';
Chan_Config.NRxAnts = 1;
Chan_Config.DopplerFreq = 5;
Chan_Config.MIMOCorrelation = 'Low';
Chan_Config.Seed = 10;
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
DL_Num_Rx_Bit=sum(DL_Config.PDSCH.TrBlkSizes(mod(0:Num_SubFrame-1,10)+1));
DL_Tx_Data=randi([0,1],DL_Num_Rx_Bit,1);
Num_Frame=floor(Num_SubFrame/10);
Num_Last_SubFrame=mod(Num_SubFrame,10);
DL_Config.TotSubframes=Num_Last_SubFrame;
[DL_Tx_Waveform_Temp,DL_Tx_Grid_Temp,DL_Config]=lteRMCDLTool(DL_Config,DL_Tx_Data(Num_Frame*sum(DL_Config.PDSCH.TrBlkSizes)+1:end));
% [DL_Tx_Waveform,DL_Tx_Grid,DL_Config]=lteRMCDLTool(DL_Config,DL_Tx_Data);
DL_Config.TotSubframes=10;
Num_Sample_Per_Frame=1e-2*DL_Config.SamplingRate;
DL_Tx_Waveform=complex(zeros(Num_SubFrame*Num_Sample_Per_Frame/10,1));
DL_Tx_Grid=complex(zeros(DL_Config.NDLRB*12,Num_SubFrame*14));
for Frame_Idx=0:Num_Frame-1
    [DL_Tx_Waveform(Frame_Idx*Num_Sample_Per_Frame+(1:Num_Sample_Per_Frame)),DL_Tx_Grid(:,Frame_Idx*140+(1:140)),DL_Config]...
        =lteRMCDLTool(DL_Config,DL_Tx_Data(Frame_Idx*sum(DL_Config.PDSCH.TrBlkSizes)+(1:sum(DL_Config.PDSCH.TrBlkSizes))));
end
Frame_Idx=Frame_Idx+1;
DL_Tx_Waveform(Frame_Idx*Num_Sample_Per_Frame+1:end)=DL_Tx_Waveform_Temp;
DL_Tx_Grid(:,Frame_Idx*140+1:end)=DL_Tx_Grid_Temp;
DL_Config.TotSubframes=Num_SubFrame;
if DEBUG
    DL_Tx_Spectrum_Scope=Spectrum_Scope.clone();
    DL_Tx_Spectrum_Scope.SampleRate=DL_Config.SamplingRate;
    DL_Tx_Spectrum_Scope.step(DL_Tx_Waveform);
    DL_Tx_Spectrum_Scope.Title='Downlink Tx Spectrum';
    DL_Tx_Constellation_Diagram=Constellation_Diagram.clone();
    DL_Tx_Constellation_Diagram.Title='Downlink Tx Constellation_Diagramm';
    DL_Tx_Constellation_Diagram.step(DL_Tx_Grid(:))
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
% DL_Rx_Waveform=DL_Tx_Waveform;
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
    subplot(1,3,1)
    mesh(abs(DL_H_Est).^2)
    title('Downlink Channel Estimation')
    DL_H_Pre=lteDLPerfectChannelEstimate(DL_Config,Chan_Config);
    subplot(1,3,2)
    mesh(abs(DL_H_Pre).^2)
    title('Downlink Perfecr Channel')
    subplot(1,3,3)
    mesh(10*log10((abs(abs(DL_H_Pre).^2-abs(DL_H_Est).^2)./abs(DL_H_Pre).^2)*100))
    title('Downlink Channel Estimation Error')
end

%%
DL_Eq_Grid=lteEqualizeZF(DL_Rx_Grid,DL_H_Est);
if DEBUG
    DL_Eq_Constellation_Diagram=Constellation_Diagram.clone();
    DL_Eq_Constellation_Diagram.Title='Downlink Rx Constellation_Diagramm After Equalization';
    DL_Eq_Constellation_Diagram.step(DL_Eq_Grid(:))
end

%%
DL_EVM=lteEVM(DL_Tx_Grid(:),DL_Eq_Grid(:));
if DEBUG
    [PDSCH_EVM,PDSCH_EVM_Plot]=hPDSCHEVM(DL_Config,struct('PilotAverage','TestEVM'),DL_Rx_Waveform);
end
% disp(DL_EVM.RMS*100)

%%
DL_Config.PDSCH.NTxAnts=size(DL_Eq_Grid,3);
DL_Config.PDSCH.W=eye(size(DL_Eq_Grid,3));
DL_Rx_Data=zeros(length(DL_Tx_Data),1);
DL_Rx_Data_Idx=0;
for SubFrame_Idx=0:Num_SubFrame-1
    DL_Config.NSubframe=mod(SubFrame_Idx,10);
    Decode_Block_Size=DL_Config.PDSCH.TrBlkSizes(:,DL_Config.NSubframe+1);
    DL_Demod_Bit=ltePDSCHDecode(DL_Config, DL_Config.PDSCH,DL_Rx_Grid(:,SubFrame_Idx*14+(1:14)),DL_H_Est(:,SubFrame_Idx*14+(1:14)),DL_Noise_Est);
    [DL_Decoded_Bit,CRC]=lteDLSCHDecode(DL_Config, DL_Config.PDSCH, Decode_Block_Size, DL_Demod_Bit);
    DL_Rx_Data(DL_Rx_Data_Idx+(1:Decode_Block_Size))=cell2mat(DL_Decoded_Bit);
    isequal(cell2mat(DL_Decoded_Bit),DL_Tx_Data(DL_Rx_Data_Idx+(1:Decode_Block_Size)))
    DL_Rx_Data_Idx=DL_Rx_Data_Idx+Decode_Block_Size;
end

