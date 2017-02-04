%ltePUSCH Physical uplink shared channel
%   SYM=ltePUSCH(UE,CHS,CWS) returns a vector containing the Physical
%   Uplink Shared Channel (PUSCH) complex symbols for UE-specific settings
%   UE, PUSCH channel specific configuration CHS and the codeword or
%   codewords contained in CWS. The size of the matrix SYM is N-by-P with N
%   being the number of modulation symbols for one antenna port, and P
%   being the number of transmission antennas.
%
%   UE must be a structure including the fields:
%   NCellID   - Physical layer cell identity
%   NSubframe - Subframe number
%   RNTI      - Radio Network Temporary Identifier (16-bit)
%   NTxAnts   - Optional. Number of transmission antennas (1(default),2,4).
%
%   CHS is the PUSCH channel-specific structure including the fields:
%   Modulation - A string specifying the modulation format for one codeword
%                or a cell array of one or two strings specifying the
%                modulation formats for one or two codewords:
%                'QPSK', '16QAM', '64QAM'
%   PRBSet     - A 1- or 2-column matrix, containing the Physical
%                Resource Block (PRB) indices corresponding to the slot
%                wise resource allocations for this PUSCH.
%   NLayers    - Optional. Number of transmission layers (1(default),2,3,4)
%   Only required for NTxAnts=2 or NTxAnts=4:
%      PMI     - Optional. Scalar Precoder Matrix Indication to be used
%                during precoding. (0(default)...23), depending on NTxAnts 
%                and NLayers. (See <a 
%                href="matlab:help('lteULPMIInfo')">lteULPMIInfo</a>)
%   
%   CWS is a vector of bit values for one codeword to be modulated, or a
%   cell array containing one or two vectors of bit values corresponding to
%   the one or two codewords to be modulated.
%   
%   For PRBSet, if a column vector is provided, the resource allocation is
%   the same in both slots of the subframe; the 2-column matrix can be used
%   to specify differing PRBs for each slot in a subframe. Note that the
%   PRB indices are 0-based.
%   
%   Example:
%   % Generate PUSCH symbols for TS 36.104 Uplink FRC A3-3 3MHz.
%
%   ue.NCellID = 1;
%   ue.NSubframe = 0;
%   ue.RNTI = 1;
%   pusch.Modulation = 'QPSK';
%   pusch.PRBSet = [0:14].';
%   pusch.RV = 0;
%   frc = lteRMCUL('A3-3');
%   trBlk  = randi([0,1],frc.PUSCH.TrBlkSizes(1),1);
%   cw = lteULSCH(ue,pusch,trBlk );
%   puschSym = ltePUSCH(ue,pusch,cw);
%   
%   See also ltePUSCHDecode, ltePUSCHIndices, ltePUSCHDRS,
%   ltePUSCHDRSIndices, ltePUSCHPrecode, lteULScramble, lteULPrecode,
%   lteULSCH.

%   Copyright 2010-2016 The MathWorks, Inc.

%   Add output of modulated symbol (layers)     2017/2/4    harryzwh

function [sym, layers] = ltePUSCH(ue,chs,cw)
    
    % If input codeword is empty then the output will be an empty matrix
    if(isempty(cw))
        sym = [];
    else
        % Scramble the codeword
        scrambled = lteULScramble(ue,cw);

        % Perform modulation of scrambled bits
        if(iscell(chs.Modulation))
            ncw = size(cw,2);
            modulated = cell(1,ncw);
            for cwIdx=1:ncw
                modulated{cwIdx} = lteSymbolModulate(scrambled{cwIdx},chs.Modulation{cwIdx});
            end
        else
            modulated = lteSymbolModulate(scrambled,chs.Modulation);
        end

        % Perform layer mapping of modulated symbols
        if (~isfield(ue,'NTxAnts'))
            ue.NTxAnts=1;
            lte.internal.defaultValueWarning('NTxAnts','1');
        end    
        if (~isfield(chs,'NLayers'))
            chs.NLayers=1;
            lte.internal.defaultValueWarning('NLayers','1');
        end        
        if (chs.NLayers<1 || chs.NLayers>4)
            error('lte:error','The function call (ltePUSCH) resulted in an error: For the parameter field NLayers, the value (%d) must be within the range [1,4]',chs.NLayers);
        end
        layers = lteLayerMap(chs,modulated);  

        % Perform SC-FDMA precoding of the layers
        scfdmaPrecoded = lteULPrecode(layers,size(chs.PRBSet,1));

        % Perform spatial precoding of SC-FDMA precoded symbols
        sym = ltePUSCHPrecode(ue,chs,scfdmaPrecoded);
    end
            
end
