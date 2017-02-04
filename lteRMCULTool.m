%lteRMCULTool Uplink RMC (FRC) waveform generation
%   [WAVEFORM,GRID,RMCCFGOUT] = lteRMCULTool(...) returns WAVEFORM, the
%   generated reference measurement channel waveform and GRID representing
%   the populated resource grid for all the physical channels specified in
%   TS 36.104 Annex A. RMCCFGOUT is a structure containing information
%   about the SC-FDMA modulated waveform as described in <a href="matlab:help('lteSCFDMAInfo')">lteSCFDMAInfo</a>
%   in addition to the RMC specific configuration parameters as described 
%   in <a href="matlab:help('lteRMCUL')">lteRMCUL</a>. The RMC waveform can be configured via a graphical user
%   interface (GUI) or by passing the required input parameters in a
%   function call.
%   
%   [WAVEFORM,GRID,RMCCFGOUT] = lteRMCULTool launches a GUI for the
%   parameterization and generation of the RMC waveforms. The main function
%   output variables are specified in the GUI but they can also be assigned
%   to the WAVEFORM, GRID and RMCCFGOUT function output variables (time
%   domain waveform, resource grid and RMC configuration structure
%   respectively).
%   
%   WAVEFORM is a T-by-P matrix where T is the number of time domain
%   samples and P is the number of antennas. GRID is a 3-dimensional array
%   of resource elements for a number of subframes across all configured
%   antenna ports, as described in the <a
%   href="matlab:web([docroot '/lte/gs/data-structures.html'])">Data Structures</a> documentation. 
%   RMCCFGOUT is a structure containing information about the SC-FDMA
%   modulated waveform as well as RMC configuration parameters.
%   
%   [WAVEFORM,GRID,RMCCFGOUT] = lteRMCULTool(RC,TRDATA,DUPLEXMODE,TOTSUBFRAMES) 
%   returns the WAVEFORM, GRID and RMCCFGOUT for the default reference
%   measurement channel defined by RC (for details, see <a href="matlab:doc('lteRMCULTool')">lteRMCULTool</a>
%   using the information bits TRDATA. DUPLEXMODE and TOTSUBFRAMES are
%   optional input parameters and define the duplex mode of the generated
%   waveform and total number of subframes that make up the GRID.
%      
%   The input RC is a string representing reference measurement channel
%   number (TS 36.104) and must be one of following:
%   ('A1-1','A1-2','A1-3','A1-4','A1-5','A2-1','A2-2','A2-3','A3-1','A3-2',
%   'A3-3','A3-4','A3-5','A3-6','A3-7','A4-1','A4-2','A4-3','A4-4','A4-5',
%   'A4-6','A4-7','A4-8','A5-1','A5-2','A5-3','A5-4','A5-5','A5-6','A5-7',
%   'A7-1','A7-2','A7-3','A7-4','A7-5','A7-6','A8-1','A8-2','A8-3', 'A8-4',
%   'A8-5','A8-6','A11-1','A3-2-9RB', 'A4-3-9RB').
%   Note that the RC 'A11-1' supports TTI bundling as described in the <a href="matlab:doc('lteRMCULTool')">doc</a>
%   for both FDD and TDD modes. RC ('A3-2-9RB', 'A4-3-9RB') are custom RMC
%   configured for non-standard bandwidths but with the same code rate as
%   the standardized versions.
%   
%   TRDATA is an array or cell array containing one or two vectors of bit
%   values where each vector is the information bits stream to be coded
%   across the duration of the generation i.e. representing multiple
%   concatenated transport blocks. Internally these vectors are looped if
%   the number of bits required across all subframes of the generation
%   exceeds the length of the vectors provided. This allows for the user to
%   enter a short pattern e.g. [1; 0; 0; 1] that will be repeated as the
%   input to the transport coding. In each subframe of generation, the
%   number of data bits taken from this stream is given by the elements of
%   the TrBlkSizes matrix, a field of the PUSCH substructure of the RMC
%   configuration structure RMCCFGOUT.
%   
%   DUPLEXMODE is an optional input representing the frame structure type
%   ('FDD'(default),'TDD').
%   
%   TOTSUBFRAMES is an optional representing the total number of subframes
%   to be generated (default 10).
%   
%   [WAVEFORM,GRID,RMCCFGOUT] = lteRMCULTool(RMCCFG,TRDATA) generates the
%   WAVEFORM, GRID and RMCCFGOUT in the same way as above except it takes
%   the user defined reference channel structure RMCCFG as input parameter.
%   The reference configuration structure with default parameters can
%   easily be created with the function <a href="matlab:
%   help('lteRMCUL')">lteRMCUL</a> which is designed to 
%   generate the various RMC configuration structures as per TS 36.104
%   Annex A. This configuration structure then can be modified as per
%   requirement and can be used in the generation of WAVEFORM.
%   
%   [WAVEFORM,GRID,RMCCFGOUT] = lteRMCULTool(RMCCFG,TRDATA,CQI,RI,ACK)
%   generates the WAVEFORM, GRID and RMCCFGOUT in the same way as above but
%   with support for control information transmission on PUSCH specified in
%   vectors CQI, RI and ACK (together these three fields form UCI). The
%   vectors CQI, RI and ACK can be empty vectors if these particular
%   control information bits are not present in this transmission. The UCI
%   is encoded for PUSCH transmission using the processing defined in TS
%   36.212 Section 5.2.4, consisting of UCI coding and channel
%   interleaving. The vectors CQI, RI and ACK are not treated as data
%   streams thus each subframe will contain the same CQI, RI and ACK
%   information bits.
%   
%   Example:
%   % Generates a time domain signal txWaveform and a 2-dimensional array  
%   % of the Resource Elements txGrid for the modified RC, A1-1 as 
%   % specified in  TS 36.104. This transmission will use 16QAM modulation
%   % of QPSK.
%   
%   frc = lteRMCUL('A1-1');
%   frc.PUSCH.Modulation = '16QAM';
%   [txWaveform,txGrid,rmcCfgOut] = lteRMCULTool(frc,[1;0;0;1]);
%
%   See also lteRMCUL, lteDLConformanceTestTool, lteRMCDLTool,
%   lteTestModelTool.

%   Copyright 2010-2016 The MathWorks, Inc.

function [timeDomainSig,frame_moded,rmcconfig] = lteRMCULTool(varargin)

    if(isempty(varargin))
        if(nargout >0 )
            % GUI call if no parameter is provided
            [timeDomainSig,frame,rmcconfig] = saLteRMCULToolGUI(nargout>0);
        else
            saLteRMCULToolGUI(nargout>0);
        end
    else
        if(isstruct(varargin{1}) && (nargin >= 2))
            rmc = varargin{1};
            transportStream = varargin{2};
            
            % Get CQI bits if provided
            if(nargin>2)
                cqiBits = varargin{3};
            else
                cqiBits = [];
            end
            
            % Get RI bits if provided
            if(nargin>3)
                riBits = varargin{4};
            else
                riBits = [];
            end
            
            % Get HARQ-ACK bits if provided
            if(nargin>4)
                ackBits = varargin{5};
            else
                ackBits = [];
            end
            
            if(~isfield(rmc,'RC'))
                error('lte:error','Input structure must contain the RC number');
            end
        elseif(ischar(varargin{1}) && (nargin >= 2))
            rmcNo = varargin{1};
            transportStream = varargin{2};

            % Get the number of subframes if supplied
            if(nargin == 4)
                totSubframes = varargin{4};
            else
                totSubframes = 10;
            end
            
            if(nargin > 2)
                duplexMode = varargin{3};
            else
                duplexMode = 'FDD';
            end
            cqiBits = [];
            riBits = [];
            ackBits = [];
            
            rmc.RC = rmcNo;
            rmc.TotSubframes = totSubframes;
            rmc.DuplexMode = duplexMode;
        else
            error('lte:error','Input argument(s) is not correct. Please see help.');
        end

        if(iscell(transportStream) || isempty(transportStream))
            [rows,cols] = size(transportStream);
            Ncodewords = max(rows,cols);            % Number of codewords
        elseif(isvector(transportStream)|| isempty(transportStream))
            Ncodewords = 1;
        end
        
        if(Ncodewords == 0 || isempty(transportStream))
            Ncodewords = 1;
        end
        
        % Get full RMC configuration for given inputs
        rmc = lteRMCUL(rmc);
 
        ttiPerBundle = 1;     % Set default value (No TTI bundling)
        if strcmpi(rmc.RC,'A11-1')
            ttiPerBundle = 4; % Enable TTI bundling for A11-1 RMC
        end
        
        PUSCH = rmc.PUSCH;
        
        % Input configuration structure is passed at output without any
        % parameter update so that the waveform can be regenerated with the
        % function output configuration
        rmcconfig = rmc;
        
        % Set up data generator source(s).
        if(iscell(transportStream))
            idx = cellfun(@isempty,transportStream);
            if(isempty(transportStream))
                dataSource1=saVectorDataSource([]);
                dataSource2=saVectorDataSource([]);
                PUSCH.TrBlkSizes(:) = 0;
                PUSCH.CodedTrBlkSizes(:) = 0;
            else
                dataSource1=saVectorDataSource(transportStream{1});
                if(numel(transportStream)>1)
                    dataSource2=saVectorDataSource(transportStream{2});
                elseif(idx)
                    PUSCH.TrBlkSizes(:) = 0;
                    PUSCH.CodedTrBlkSizes(:) = 0;
                end
            end
        else
            if(isempty(transportStream))
                PUSCH.TrBlkSizes(:) = 0;
                PUSCH.CodedTrBlkSizes(:) = 0;
            end
            dataSource1=saVectorDataSource(transportStream);
        end
        
        frame = [];
        reGrid = [];
        cBlkSizeIdx = 0;
        frame_moded=[];
        
        % Empty subframe construction for given configuration filled with all
        % zeros.
        subframe = lteULResourceGrid(rmc);
        subframe_moded=subframe;
        
        totalSubframes = rmc.TotSubframes;
        if(~isnumeric(totalSubframes) || isempty(totalSubframes) || totalSubframes<0)
            error('lte:error','The function call resulted in an error: The total number of subframes must be a positive integer');
        end
        
        % Check if DL cyclic prefix is specified
        if(~isfield(rmc,'CyclicPrefix'))
            rmc.CyclicPrefix='Normal';
        end
        
        % Set default value to SerialCat if its not defined.
        if(~isfield(rmc,'SerialCat'))
            rmc.SerialCat = true;
        end
        serialCat = rmc.SerialCat;
        
        if(~isfield(PUSCH,'RVSeq') || (isfield(PUSCH,'RVSeq') && isempty(PUSCH.RVSeq)))
            if(Ncodewords == 2)
                PUSCH.RVSeq = [0;0];
            elseif(Ncodewords == 1)
                PUSCH.RVSeq = 0;
            end
        end
        
        % HARQ setup initialization
        noHarqProcesses = PUSCH.NHARQProcesses;
        harqTableIdx = 1;
        % If TTI bundling is enabled ('A11-1'), each HARQ process will
        % transmit 'ttiPerBundle' times with different RV values before
        % moving on to the next HARQ process. Define the HARQ table with
        % all HARQ processes repeating 'ttiPerBundle' times
        harqTable = reshape(repmat(1:noHarqProcesses,ttiPerBundle,1),1,noHarqProcesses*ttiPerBundle);
        
        harqprocess.data = struct('blk1',[],'blk2',[]);
        harqprocess.RVIdx = ones(Ncodewords,1);
        harqProcesses(1:max(harqTable)) = harqprocess;
        newData = [1 1];
        
        [nRows,nColms] = size(PUSCH.TrBlkSizes);
        tempBlkSize = PUSCH.TrBlkSizes;
        if(Ncodewords==2)
            if(nColms==2)
                tempBlkSize = tempBlkSize.';
            end
        else
            if(nRows>2)
                tempBlkSize = tempBlkSize.';
            end
        end
        
        % Get the absolute subframe number from NSubframe and NFrame
        NSubframe = rmc.NFrame*10+rmc.NSubframe;
        for subframeIdx=NSubframe:(NSubframe+totalSubframes)-1

            % Update subframe number and clear subframe
            subframe(:)=0;
            subframe_moded(:)=0;
            
            rmc.NSubframe = mod(subframeIdx,10);
            rmc.NFrame = floor(subframeIdx/10);

            % If this subframe is a cell-specific SRS subframe, configure
            % the PUCCH for shortened transmission.
            if(isfield(rmc,'SRS') && isstruct(rmc.SRS))
                srsInfo = lteSRSInfo(rmc,rmc.SRS);
                rmc.Shortened = srsInfo.IsSRSSubframe;
            end
            
            info=lteDuplexingInfo(rmc);
            transportBlock = {[]}; % Initialize cw1
            if (info.NSymbolsUL)
                % Get the transport block size for the current subframe
                if(Ncodewords==2)
                    transportBlkSize = tempBlkSize(:,mod(rmc.NSubframe,size(tempBlkSize,2)) + 1);
                    if(size(transportBlkSize,1)==2)
                        transportBlkSize = transportBlkSize.';
                    elseif(size(transportBlkSize,1)==1 && size(transportBlkSize,2)==1)
                        transportBlkSize = [transportBlkSize transportBlkSize]; %#ok<AGROW>
                    end
                    
                    if(isempty(transportBlkSize(1)) || isempty(transportBlkSize(2)))
                        effectiveCodewords = 1;
                        transportBlkSize(1) = transportBlkSize(1)*layersPerCW(PUSCH.NLayers,effectiveCodewords);
                        transportBlkSize(2) = transportBlkSize(2)*layersPerCW2(PUSCH.NLayers,effectiveCodewords);
                    end
                    transportBlock{2} = [];% Initialize cw2
                else
                    transportBlkSize = tempBlkSize(1,mod(rmc.NSubframe,size(tempBlkSize,2)) + 1);
                end
                
                % Generate PUSCH, PUSCH can only be transmitted in UL subframes
                if(strcmpi(info.SubframeType,'Uplink')) && any(transportBlkSize~=0)
                    harqIdx =  harqTable(harqTableIdx);
                    tempRv(1,1) = PUSCH.RVSeq(1,harqProcesses(harqIdx).RVIdx(1));
                    if(Ncodewords==2)
                        tempRv(2,1) = PUSCH.RVSeq(end,harqProcesses(harqIdx).RVIdx(2));
                    end
                    PUSCH.RV = tempRv(:,1);
                    newData(:) = (harqProcesses(harqIdx).RVIdx == 1);

                    % UL-SCH transport block size of configured RMC as in TS 36.104.
                    if(Ncodewords==1 && (newData(1) || isempty(harqProcesses(harqIdx).data.blk1)))
                        harqProcesses(harqIdx).data.blk1 = dataSource1.getData(transportBlkSize);
                    elseif(Ncodewords == 2)
                        transportBlkSize(idx) = 0;
                        if(isempty(PUSCH.RV))
                            PUSCH.RV = 0;
                        end
                        if(newData(1) || isempty(harqProcesses(harqIdx).data.blk1))
                            harqProcesses(harqIdx).data.blk1 = dataSource1.getData(transportBlkSize(1));
                        end
                        if(newData(2) || isempty(harqProcesses(harqIdx).data.blk2))
                            harqProcesses(harqIdx).data.blk2 = dataSource2.getData(transportBlkSize(2));
                        end
                    end
                    if(transportBlkSize(1)~=0)
                        transportBlock = {harqProcesses(harqIdx).data.blk1};
                    end
                    if(Ncodewords == 2)
                        if(transportBlkSize(2)~=0)
                            transportBlock{2} = harqProcesses(harqIdx).data.blk2;
                        end
                    end
                    % If transport block size wasn't 0, update rv and
                    % harq table index
                    harqProcesses(harqIdx).RVIdx(transportBlkSize~=0) = mod(harqProcesses(harqIdx).RVIdx(transportBlkSize~=0) ,size(PUSCH.RVSeq,2))+1;
                    harqTableIdx = mod(harqTableIdx,length(harqTable)) + 1;
                    
                    % Generating UL-SCH coded bits by performing complete channel coding including
                    % CRC calculation, code block segmentation and CRC attachment, turbo
                    % coding, rate matching and code block concatenation.
                    codedTrBlock = lteULSCH(rmc,PUSCH,transportBlock,cqiBits,riBits,ackBits);

                    cBlkSizeIdx = cBlkSizeIdx+1;
                    if(Ncodewords==2)
                        PUSCH.TrBlkSizes(1,cBlkSizeIdx) = transportBlkSize(1);
                        PUSCH.TrBlkSizes(2,cBlkSizeIdx) = transportBlkSize(2);
                        PUSCH.CodedTrBlkSizes(1,cBlkSizeIdx) = length(codedTrBlock{1});
                        PUSCH.CodedTrBlkSizes(2,cBlkSizeIdx) = length(codedTrBlock{2});
                    else
                        PUSCH.TrBlkSizes(1,cBlkSizeIdx) = transportBlkSize;
                        PUSCH.CodedTrBlkSizes(1,cBlkSizeIdx) = length(codedTrBlock{1});
                    end

                    % Complex-valued modulated symbol generation for PUSCH. This involves
                    % scrambling, modulation and precoding processes.
                    if(Ncodewords==2 && isempty(codedTrBlock{1}) && ~isempty(codedTrBlock{2}))
                        % Swapping codeword and corresponding modulation scheme
                        [puschSymbols, puschSymbols_moded] = ltePUSCH(rmc,setfield(PUSCH,'Modulation',{PUSCH.Modulation{2} PUSCH.Modulation{1}}),{codedTrBlock{2} []}); %#ok<SFLD>
                    else
                        [puschSymbols, puschSymbols_moded] = ltePUSCH(rmc,PUSCH,codedTrBlock);
                    end
                    
                    % PUSCH symbols mapping on to the resource grid 
                    puschIndices = ltePUSCHIndices(rmc,PUSCH);
                    subframe(puschIndices) = puschSymbols;  
                    subframe_moded(puschIndices) = puschSymbols_moded; 
                    
                    % PUSCH DRS symbol creation and mapping on to resource grid
                    puschDrsSeq = ltePUSCHDRS(rmc,PUSCH);
                    puschDrsSeqIndices = ltePUSCHDRSIndices(rmc,PUSCH);                
                    subframe(puschDrsSeqIndices) = puschDrsSeq;        
                end
                
                if(isfield(rmc,'SRS') && isstruct(rmc.SRS) && srsInfo.IsSRSSubframe)
                    % Transmit SRS (if active under UE-specific SRS and cell-specific configuration)
                    [srsIndices,srsIndicesInfo] = lteSRSIndices(rmc,rmc.SRS);
                    srsSymbols = lteSRS(rmc,rmc.SRS);
                    if (rmc.SRS.NTxAnts==1 && rmc.NTxAnts>1)
                        subframe(offsetIndices(rmc,srsIndices,srsIndicesInfo.Port)) = srsSymbols;
                    else
                        subframe(srsIndices) = srsSymbols;
                    end
                end
            end
            % concatenate subframes to form a complete frame
            frame = cat(2,frame,subframe);
            frame_moded=cat(2,frame_moded,subframe_moded);
            if(~serialCat)
                reGrid(:,:,:,mod(subframeIdx-NSubframe,NSubframe+totalSubframes)+1) = subframe; %#ok<AGROW>
            end
        end
        
        if(~isfield(rmc,'Windowing'))
            rmc.Windowing = 0;
        end
        
        % Time Domain mapping by performing SC-FDMA modulation for uplink
        % symbols.
        [timeDomainSig,infoScfdma] = lteSCFDMAModulate(rmc,frame);
        if(~serialCat)
            frame = reGrid;
        end
        rmcconfig.SamplingRate = infoScfdma.SamplingRate;
        rmcconfig.Nfft = infoScfdma.Nfft;
        rmcconfig.Windowing = infoScfdma.Windowing;
        rmcconfig.PUSCH.HARQProcessSequence = getHARQTable(rmcconfig,ttiPerBundle);
    end
end

% This function calculates the number of layers per codeword.
function nLyrs = layersPerCW(layers,nCWs)
    nLyrs = floor(layers/nCWs);
end

function nLyrs = layersPerCW2(layers,nCWs)
    nLyrs = ceil(layers/nCWs);
end

% obtains indices 'out' corresponding to the indices 'in' offset to address
% antenna plane 'p'. 
function out = offsetIndices(ue,in,p)
    griddims=lteULResourceGridSize(ue);
    K=griddims(1);          % number of subcarriers
    L=griddims(2);          % number of OFDM symbols in a subframe    
    out=in+(K*L*p);
end

% Function to calculate the HARQ table assuming all uplink subframes are
% carrying data and have the same transport block size
function harqTable = getHARQTable(rmc,ttiPerBundle)
    noHarqProcesses = rmc.PUSCH.NHARQProcesses;
    % Define the default values for cyclic prefix
    if ~isfield(rmc,'CyclicPrefixUL')
        rmc.CyclicPrefixUL = 'Normal';
    end
    if ~isfield(rmc,'CyclicPrefix')
        rmc.CyclicPrefix = 'Normal';
    end
    info = arrayfun(@(x)lteDuplexingInfo(setfield(rmc,'NSubframe',x)),0:9); %#ok<SFLD>
    activesfs = arrayfun(@(x)strcmpi(x.SubframeType,'Uplink'),info);
    harqTable = ones(10,lcm(noHarqProcesses*ttiPerBundle,sum(activesfs))/sum(activesfs))*-1;
    harqTable(activesfs==0,:) = 0; % Non-transmititng subframes
    harqTable(harqTable==-1) = repmat(repmat(1:noHarqProcesses,ttiPerBundle,1),1,lcm(noHarqProcesses*ttiPerBundle,sum(activesfs))/(noHarqProcesses*ttiPerBundle));
	harqTable = harqTable(:).';
end