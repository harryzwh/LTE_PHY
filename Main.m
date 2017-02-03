%%
clc
clear
close all

%%
Num_Figure=0;
DEBUG=true;

%%
UL_RMC='A5-5';
DL_RMC='R.7';
Duplex_Mode='FDD';
Num_SubFrame=20;

%%
Spectrum_Scope=dsp.SpectrumAnalyzer();
Constellation_Diagram=comm.ConstellationDiagram();
Constellation_Diagram.ShowReferenceConstellation=false;

%%
UL_Config=lteRMCUL(UL_RMC,Duplex_Mode,Num_SubFrame);
UL_Num_Tx_Bit=sum(UL_Config.PUSCH.TrBlkSizes);
UL_Tx_Data=randi([0,1],UL_Num_Tx_Bit,1);
[UL_Tx_Waveform,UL_Grid,UL_Config]=lteRMCULTool(UL_Config,UL_Tx_Data);
if DEBUG
    UL_Spectrum_Scope=Spectrum_Scope.clone();
    UL_Spectrum_Scope.SampleRate=UL_Config.SamplingRate;
    UL_Spectrum_Scope.step(UL_Tx_Waveform);
    UL_Spectrum_Scope.Title='Uplink Spectrum';
    UL_Constellation_Diagram=Constellation_Diagram.clone();
    UL_Constellation_Diagram.Title='Uplink Constellation_Diagramm';
    UL_Constellation_Diagram.step(UL_Grid(:))
end

%%
DL_Config=lteRMCDL(DL_RMC,Duplex_Mode,Num_SubFrame);
DL_Num_Rx_Bit=sum(DL_Config.PDSCH.TrBlkSizes);
DL_Rx_Data=randi([0,1],DL_Num_Rx_Bit,1);
[DL_Tx_Waveform,DL_Grid,DL_Config]=lteRMCDLTool(DL_Config,DL_Rx_Data);
if DEBUG
    DL_Spectrum_Scope=Spectrum_Scope.clone();
    DL_Spectrum_Scope.SampleRate=DL_Config.SamplingRate;
    DL_Spectrum_Scope.step(DL_Tx_Waveform);
    DL_Spectrum_Scope.Title='Downlink Spectrum';
    DL_Constellation_Diagram=Constellation_Diagram.clone();
    DL_Constellation_Diagram.Title='Downlink Constellation_Diagramm';
    DL_Constellation_Diagram.step(DL_Grid(:))
end
