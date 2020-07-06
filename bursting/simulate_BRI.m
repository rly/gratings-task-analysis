% script to model bursting as measured with SPNA/BRI approach
% this script was written to examine why we see such low BRI values. We
% need to make sure that the script is working properly so I will simulate
% bursting spike trains and feed them into the BRI pipeline

% SPNA takes binary data
% create long data
clear

for runi = 1:50
    tmp_data = zeros(250,800); % trials, time points
    burst = [1 0 0 0 1 0 0 1];
    for triali = 1:size(tmp_data,1)
        start_p = randperm(100,1) + randperm(200,1);
        if start_p + length(burst) < size(tmp_data,2)
            tmp_data(triali,start_p:(start_p+length(burst)-1)) = burst;
            start_p = find(tmp_data(triali,:),1,'last') + randperm(100,1) + randperm(200,1);
        end
        if start_p + length(burst) < size(tmp_data,2)
            tmp_data(triali,start_p:(start_p+length(burst)-1)) = burst;
            start_p = find(tmp_data(triali,:),1,'last') + randperm(100,1) + randperm(200,1);
        end
        if start_p + length(burst) < size(tmp_data,2)
            tmp_data(triali,start_p:(start_p+length(burst)-1)) = burst;
            start_p = find(tmp_data(triali,:),1,'last') + randperm(100,1) + randperm(200,1);
        end
        if start_p + length(burst) < size(tmp_data,2)
            tmp_data(triali,start_p:(start_p+length(burst)-1)) = burst;
            start_p = find(tmp_data(triali,:),1,'last') + randperm(100,1) + randperm(200,1);
        end
        if start_p + length(burst) < size(tmp_data,2)
            tmp_data(triali,start_p:(start_p+length(burst)-1)) = burst;
            start_p = find(tmp_data(triali,:),1,'last') + randperm(100,1) + randperm(200,1);
        end
        if start_p + length(burst) < size(tmp_data,2)
            tmp_data(triali,start_p:(start_p+length(burst)-1)) = burst;
            start_p = find(tmp_data(triali,:),1,'last') + randperm(100,1) + randperm(200,1);
        end
        if start_p + length(burst) < size(tmp_data,2)
            tmp_data(triali,start_p:(start_p+length(burst)-1)) = burst;
            start_p = find(tmp_data(triali,:),1,'last') + randperm(100,1) + randperm(200,1);
        end
    end
    
    SC_ptrial(runi,:) = sum(tmp_data,2);
    FR(runi) = mean(sum(tmp_data,2)/0.8);

    spikeFs = 1000;
    maxLag = spikeFs / 4; % 20ms
    [SPNA, meanAutoCorr, meanCrossCorr, sdCrossCorr] = computeSPNA(tmp_data, maxLag);

    % calculate BRI
    BRILagStartInd = maxLag + 1 + spikeFs * 0.001;
    BRILagEndInd = maxLag + 1 + spikeFs * 0.004;
    BRI(runi) = mean(SPNA(BRILagStartInd:BRILagEndInd));
end

figure
plot(BRI,FR,'o')
ylabel('Firing rate (Hz)'); xlabel('BRI')
title(['Burst [' num2str(burst) '] trial length = ' num2str(size(tmp_data,2)) 'ms'])

%% Add Poisson spikes
clear
runs_burst = randperm(50,25); 
freq_options = [5 13 13 15 15];
for runi = 1:110
    freqid = randperm(numel(freq_options),1);
  
    tmp_data = zeros(250,300); % trials, time points
    
    if ismember(runi,runs_burst)
        burst = [1 0 0 0 1];
        for triali = 1:5 % size(tmp_data,1) % randperm(250,5) %
            start_p = randperm(150,1);% + randperm(200,1);
%             if mod(triali,2)
                if start_p + length(burst) < size(tmp_data,2)
                    tmp_data(triali,start_p:(start_p+length(burst)-1)) = burst;
                    start_p = find(tmp_data(triali,:),1,'last') + randperm(100,1) + randperm(200,1);
                end
%             else
%                 start_p = start_p + 4;
%                 if start_p + length(burst) < size(tmp_data,2)
%                     tmp_data(triali,start_p:(start_p+length(burst)-1)) = burst;
%                     start_p = find(tmp_data(triali,:),1,'last') + randperm(100,1) + randperm(200,1);
%                 end
%             end
                    
    %         if start_p + length(burst) < size(tmp_data,2)
    %             tmp_data(triali,start_p:(start_p+length(burst)-1)) = burst;
    %             start_p = find(tmp_data(triali,:),1,'last') + randperm(100,1) + randperm(200,1);
    %         end
    %         if start_p + length(burst) < size(tmp_data,2)
    %             tmp_data(triali,start_p:(start_p+length(burst)-1)) = burst;
    %             start_p = find(tmp_data(triali,:),1,'last') + randperm(100,1) + randperm(200,1);
    %         end
    %         if start_p + length(burst) < size(tmp_data,2)
    %             tmp_data(triali,start_p:(start_p+length(burst)-1)) = burst;
    %             start_p = find(tmp_data(triali,:),1,'last') + randperm(100,1) + randperm(200,1);
    %         end
    %         if start_p + length(burst) < size(tmp_data,2)
    %             tmp_data(triali,start_p:(start_p+length(burst)-1)) = burst;
    %             start_p = find(tmp_data(triali,:),1,'last') + randperm(100,1) + randperm(200,1);
    %         end
    %         if start_p + length(burst) < size(tmp_data,2)
    %             tmp_data(triali,start_p:(start_p+length(burst)-1)) = burst;
    %             start_p = find(tmp_data(triali,:),1,'last') + randperm(100,1) + randperm(200,1);
    %         end
    %         if start_p + length(burst) < size(tmp_data,2)
    %             tmp_data(triali,start_p:(start_p+length(burst)-1)) = burst;
    %             start_p = find(tmp_data(triali,:),1,'last') + randperm(100,1) + randperm(200,1);
    %         end
        end
    end
    
    fr = freq_options(freqid);%max(sum(tmp_data,2));%1; % Hz
    dt = 1/1000; % s
    nBins = 300; % 10 ms spike train
    for triali = 1:250
        tmp_data(triali,:) = rand(1, nBins) < fr*dt;
    end
  
    tmp_data(126:end,1:end-2) = tmp_data(1:125,3:end);
    
%     figure; imagesc(tmp_data)
    SC_ptrial(runi,:) = sum(tmp_data,2);
    FR(runi) = mean(sum(tmp_data,2)/(nBins/1000));

    spikeFs = 1000;
    maxLag = spikeFs / 4; % 20ms
    [SPNA, meanAutoCorr, meanCrossCorr, sdCrossCorr] = computeSPNA(tmp_data, maxLag);

    % calculate BRI
    BRILagStartInd = maxLag + 1 + spikeFs * 0.001;
    BRILagEndInd = maxLag + 1 + spikeFs * 0.004;
    BRI(runi) = mean(SPNA(BRILagStartInd:BRILagEndInd));
end

figure
plot(BRI,FR,'o')
ylabel('Firing rate (Hz)'); xlabel('BRI')
title(['Burst [' num2str(burst) '] trial length = ' num2str(size(tmp_data,2)) 'ms'])

figure 
histogram(BRI,-0.6:0.1:0.6)



% eof