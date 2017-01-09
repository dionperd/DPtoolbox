function [Labels, INDs, fc, nm] = DPconMat2conVecCFG(Nf,Nch,f,chanLabels,freqMode,chanMode)



if strcmpi(freqMode,'am2ph')
    
    Nloops = (Nf*(Nf-1)/2)*(Nch^2);
    INDs = nan(Nloops,4);
    fc = nan(Nloops,2);
    nm = nan(Nloops,2);
    Labels = cell(Nloops,4);
    iL=0;
    for iF1=2:Nf; %for all frequencies
        
        for iF2=1:(iF1-1); %for all lower frequencies
            
            for iC1=1:Nch; %for all channel pairs
                
                for iC2=1:Nch; %Self couplings are allowed, because it will always be amplitude (1st argument) to phase (2nd argument)
                    
                    iL=iL+1;
                    
                    INDs(iL,:) = [iF1 iF2 iC1 iC2]; %the indices of center frequencies and channels
                    
                    fc(iL,:) = [f(iF1) f(iF2)]; %the central frequencies
                    
                    Labels(iL,:) = {num2str(f(iF1)),num2str(f(iF2)),chanLabels{iC1},chanLabels{iC2}}; %frequency and channel labels
                end
            end
        end
    end

    
elseif strcmpi(freqMode,'CFC')
    
    Nn = Nf*Nch; %number of different nodes    
    Ch = repmat([1:Nch].',[Nf,1]);
    Fr = repmat([1:Nf],[Nch,1]);
    Fr = Fr(:);
    
    if strcmpi(chanMode,'directed')
        
        %All to all frequencies' and all to all channels' connections
        %expect for self-connections
        Nloops = Nn*(Nn-1);
        INDs = nan(Nloops,4);
        fc = nan(Nloops,2);
        nm = nan(Nloops,2);
        Labels = cell(Nloops,4);
        
        iL=0;
        for iN1=1:Nn; %for all nodes
            
            iF1 = Fr(iN1);
            iC1 = Ch(iN1);
            
            for iN2=1:Nn; %for all nodes
                
                if (iN1~=iN2) %No self couplings are allowed
                    
                    iL=iL+1;
                    
                    
                    iF2 = Fr(iN2);
                    iC2 = Ch(iN2);
                    
                    INDs(iL,:) = [iF1 iF2 iC1 iC2]; %the indices of center frequencies and channels
                    
                    fc(iL,:) = [f(iF1) f(iF2)]; %the central frequencies
                    
                    fint =  floor(10^6*fc(iL,:)); %to account for non integer frequencies with a precision of 10^(-6)
                    nm(iL,:) = fint([2,1])/gcd(fint(1),fint(2)); %The n:m ratios
                    
                    Labels(iL,:) = {num2str(f(iF1)),num2str(f(iF2)),chanLabels{iC1},chanLabels{iC2}}; %frequency and channel labels
                    
                end
                
            end
        end
        
        if nargout>3
            nm=ones(Nloops,2);
        end
        
    else
        
        %All to all frequencies' and one way channels' connections
        Nloops = Nn*(Nn-1)/2;
        INDs = nan(Nloops,4);
        fc = nan(Nloops,2);
        nm = nan(Nloops,2);
        Labels = cell(Nloops,4);

        iL=0;
        for iN1=1:Nn; %for all nodes
            
            iF1 = Fr(iN1);
            iC1 = Ch(iN1);
            
            for iN2=iN1+1:Nn; %for all the nodes with higher index
                
                    
                    iL=iL+1;
                    
                    iF2 = Fr(iN2);
                    iC2 = Ch(iN2);
                    
                    INDs(iL,:) = [iF1 iF2 iC1 iC2]; %the indices of center frequencies and channels
                    
                    fc(iL,:) = [f(iF1) f(iF2)]; %the central frequencies
                    
                    fint =  floor(10^6*fc(iL,:)); %to account for non integer frequencies with a precision of 10^(-6)
                    nm(iL,:) = fint([2,1])/gcd(fint(1),fint(2)); %The n:m ratios
                    
                    Labels(iL,:) = {num2str(f(iF1)),num2str(f(iF2)),chanLabels{iC1},chanLabels{iC2}}; %frequency and channel labels
                                    
            end
        end
        
    end
    
    
else %('SFC', i.e., same frequency phase to phase coupling)
    
    if strcmpi(chanMode,'directed')
        
        %All to all channels' connections expect for
        %self-connections, for each frequency band
        Nloops = Nf*(Nch*(Nch-1));
        INDs = nan(Nloops,3);
        fc = nan(Nloops,1);
        Labels = cell(Nloops,3);
        iL=0;
        for iF=1:Nf; %for all frequencies
            
            for iC1=1:Nch; %for all channel pairs
                
                for iC2=1:Nch;
                    
                    if (iC1~=iC2)%No self couplings are allowed
                        
                        iL=iL+1;
                        
                        INDs(iL,:) = [iF iC1 iC2]; %the indices of center frequencies and channels
                        
                        fc(iL,:) = f(iF); %the central frequencies
                                                
                        Labels(iL,:) = {num2str(f(iF)),chanLabels{iC1},chanLabels{iC2}}; %frequency and channel labels
                        
                    end
                end
            end
        end
        if nargout>3
            nm=ones(Nloops,2);
        end
        
    else
        %One way channels' connections, for each frequency band
        Nloops = Nf*(Nch*(Nch-1)/2);
        INDs = nan(Nloops,3);
        fc = nan(Nloops,1);
        Labels = cell(Nloops,3);
        iL=0;
        for iF=1:Nf; %for all frequencies
            
            for iC1=1:Nch; %for all channel pairs
                
                for iC2=iC1+1:Nch; %for one way channel connections
                    
                    iL=iL+1;
                    
                    INDs(iL,:) = [iF iC1 iC2]; %the indices of center frequencies and channels
                    
                    fc(iL,:) = f(iF); %the central frequencies
                                        
                    Labels(iL,:) = {num2str(f(iF)),chanLabels{iC1},chanLabels{iC2}}; %frequency and channel labels
                    
                end
            end
        end
        if nargout>3
            nm=ones(Nloops,2);
        end
    end
end




end

