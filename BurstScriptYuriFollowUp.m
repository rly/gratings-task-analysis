clear
close all
outputDir = '/Volumes/HDD1/Bursting/LTS_HTB';
v = 14;

minFiringRateOverall = 5; % Hz

homeDir = '/Volumes/HDD1/raw_data/';
cd(homeDir)
rawDataLocs = {'C110524',	'C110614',	'C110708',	 'C110811',	'L101018',...
    'L101102',	'L101214',	'L110419',	'L110519',   'C110527',	'C110615',...
    'C110713',	'C110812',	'L101019',	'L101103',	 'L101215',	'L110420',...
    'L110523',  'C110531',	'C110616',	'C110720',	 'L101007',	'L101020',...
    'L101105',	'L101216',	'L110421',	'L110524',   'C110601',	'C110617',...
    'C110722',	'L101008',	'L101021',	'L101119',	 'L101221',	'L110422',...
    'L110531',  'C110603',	'C110623',	'C110727',	 'L101011',	'L101022',...
    'L101123',	'L101222',	'L110426',	'L110711',   'C110608',	'C110624',...
    'C110728',	'L101012',	'L101025',	'L101124',	 'L110411',	'L110429',...
    'L110811',  'C110610',	'C110630',	'C110804',	 'L101013',	'L101027',...
    'L101208',	'L110412',	'L110502',	'L110812',   'C110613',	'C110701',...
    'C110809',	'L101014',	'L101029',	'L101209',	 'L110413',	'L110503'};
rawDataLocs = sort(rawDataLocs);
for sessioni = 1:numel(rawDataLocs)
    clearvars -except sessioni homeDir rawDataLocs outputDir minFiringRateOverall v
    
    cd(rawDataLocs{sessioni})
    sessionInfo = dir('*_allSpikes.mat');
    if ~isempty(sessionInfo)
        load(sessionInfo.name,'sig00*')
        if exist('sig001a','var')
            units = whos('sig001*');
            nUnits = length(units);
            for uniti = 1:nUnits
                spikeTimes = eval(units(uniti).name);
                totalTimeOverall = spikeTimes(end) - spikeTimes(1);
                unitName = units(uniti).name;
                firingRateOverall = numel(spikeTimes) / totalTimeOverall;
                fprintf('SPNA %s (%d/%d = %d%%)... \n', unitName, uniti, ...
                    nUnits, round(uniti/nUnits*100));

                %     saveFileName = sprintf('%s/%s-%s-evokedSpiking-v%d.mat', ...
                %             processedDataDir, unitName, blockName, v);
                if firingRateOverall >= minFiringRateOverall
                    fprintf('\tOverall firing rate = %0.2f Hz > minimum firing rate = %0.2f Hz in these blocks.\n', ...
                        firingRateOverall, minFiringRateOverall);
                    %         fprintf('\tComputing evoked spiking and writing file %s...\n', saveFileName);

                    allDiffSpikeTimes = diff(spikeTimes); % get the ISIs
                    allDiffSpikeTimesPre = [NaN; allDiffSpikeTimes]; % miss align ISIs by adding NaN to create pre and postISI
                    allDiffSpikeTimesPost = [allDiffSpikeTimes; NaN];
                    allDiffSpikeTimesPre = log10(allDiffSpikeTimesPre*1000); % logISI in ms
                    allDiffSpikeTimesPost = log10(allDiffSpikeTimesPost*1000);

                    FigH = figure;
                    % return plot
                    subplot(4,3,[1 4 7 10])
                    plot(allDiffSpikeTimesPre,allDiffSpikeTimesPost,'.','MarkerSize', 1)
                    hold on
                    postTimesUnderThreshShort = find(allDiffSpikeTimesPost < log10(4));
                    x2 = allDiffSpikeTimesPre(postTimesUnderThreshShort);
                    y2 = allDiffSpikeTimesPost(postTimesUnderThreshShort);
                    bothThresh = find(allDiffSpikeTimesPre(postTimesUnderThreshShort) > log10(100));
                    plot(x2(bothThresh),y2(bothThresh),'k.','MarkerSize', 1)
                    clear x2 y2
                    postTimesUnderHTBThresh = find(allDiffSpikeTimesPost < log10(25));
                    x2 = allDiffSpikeTimesPre(postTimesUnderHTBThresh);
                    y2 = allDiffSpikeTimesPost(postTimesUnderHTBThresh);
                    bothThreshHTB = find(allDiffSpikeTimesPre(postTimesUnderHTBThresh) < log10(25));
                    plot(x2(bothThreshHTB),y2(bothThreshHTB),'r.','MarkerSize', 1)
                    title('Return plot: All data')
                    hold on
                    plot([2 2],[log10(1) log10(4)],'--k')
                    plot([2 4],[log10(4) log10(4)],'--k')
                    %             plot([3 3],[log10(1) log10(4)],'--k')
                    plot([log10(25) log10(25)],[0 log10(25)],'--k')
                    plot([0 log10(25)],[log10(25) log10(25)],'--k')
                    xlim([0 4])
                    ylim([0 4])
                    xticks(0:4)
                    xticklabels({'1', '10', '100', '1000'})
                    xlabel('PreISI in ms on log scale')
                    yticks(0:4)
                    yticklabels({'1', '10', '100', '1000'})
                    ylabel('PostISI in ms on log scale')
                    
                    subplot(4,3,[2 5 8])
                    title('LTS Burst')
                    spikeIndicesLTSBurstStart = postTimesUnderThreshShort(bothThresh);
                    for spikei = 1:length(spikeIndicesLTSBurstStart)
                        if spikeIndicesLTSBurstStart(spikei) > 2 && spikeIndicesLTSBurstStart(spikei) < size(spikeTimes,1) - 5
                        plot([(spikeTimes(spikeIndicesLTSBurstStart(spikei)-2:spikeIndicesLTSBurstStart(spikei)+5)) - spikeTimes(spikeIndicesLTSBurstStart(spikei)), ...
                            (spikeTimes(spikeIndicesLTSBurstStart(spikei)-2:spikeIndicesLTSBurstStart(spikei)+5)) - spikeTimes(spikeIndicesLTSBurstStart(spikei))],[spikei spikei+1],'k')
                        hold on
                        plot([spikeTimes(spikeIndicesLTSBurstStart(spikei):spikeIndicesLTSBurstStart(spikei)) - spikeTimes(spikeIndicesLTSBurstStart(spikei)), ...
                            spikeTimes(spikeIndicesLTSBurstStart(spikei):spikeIndicesLTSBurstStart(spikei)) - spikeTimes(spikeIndicesLTSBurstStart(spikei))],[spikei spikei+1],'r')
                        else
                            continue
                        end
                    end
                    xlim([-0.5 0.5])
                    
                    prepBurstToPSTH = zeros(length(spikeIndicesLTSBurstStart),8);
                    time = linspace(-0.5,0.5,1001);
                    binaryPrepBurstPSTH = zeros(length(spikeIndicesLTSBurstStart),size(time,2));
                    for spikei = 1:length(spikeIndicesLTSBurstStart)
                        if spikeIndicesLTSBurstStart(spikei) > 2 && spikeIndicesLTSBurstStart(spikei) < size(spikeTimes,1) - 5
                            prepBurstToPSTH(spikei,:) = ((spikeTimes(spikeIndicesLTSBurstStart(spikei)-2:spikeIndicesLTSBurstStart(spikei)+5)) - spikeTimes(spikeIndicesLTSBurstStart(spikei)))';
                            for spikeii = 1:8
                                [~,idx(spikeii)] = min(abs(time-prepBurstToPSTH(spikei,spikeii)));
                                binaryPrepBurstPSTH(spikei,idx) = 1;
                            end
                        end
                    end
                    subplot(4,3,11)
                    plot(time,sum(binaryPrepBurstPSTH,1))
                    smoothPSTH = smooth(sum(binaryPrepBurstPSTH,1),11,'lowess');
                    hold on; plot(time,smoothPSTH,'g','LineWidth',2)
                    ylim([-1 20]); xlim([-0.45 0.45])
                    
                    subplot(4,3,[3 6 9])
                    title('HT Burst')
                    spikeIndicesHTBBurstStart = postTimesUnderHTBThresh(bothThreshHTB);
                    newHTB = find(diff(spikeIndicesHTBBurstStart)>1);
                    for spikei = 1:length(newHTB)
                        if spikeIndicesHTBBurstStart(newHTB(spikei)) > 5 && spikeIndicesHTBBurstStart(newHTB(spikei)) < size(spikeTimes,1) - 2
                        plot([(spikeTimes(spikeIndicesHTBBurstStart(newHTB(spikei))-5:spikeIndicesHTBBurstStart(newHTB(spikei))+2)) - spikeTimes(spikeIndicesHTBBurstStart(newHTB(spikei))), ...
                            (spikeTimes(spikeIndicesHTBBurstStart(newHTB(spikei))-5:spikeIndicesHTBBurstStart(newHTB(spikei))+2)) - spikeTimes(spikeIndicesHTBBurstStart(newHTB(spikei)))],[spikei spikei+1],'k')
                        hold on
                        plot([spikeTimes(spikeIndicesHTBBurstStart(newHTB(spikei)):spikeIndicesHTBBurstStart(newHTB(spikei))) - spikeTimes(spikeIndicesHTBBurstStart(newHTB(spikei))), ...
                            spikeTimes(spikeIndicesHTBBurstStart(newHTB(spikei)):spikeIndicesHTBBurstStart(newHTB(spikei))) - spikeTimes(spikeIndicesHTBBurstStart(newHTB(spikei)))],[spikei spikei+1],'r')
                        else
                            continue
                        end
                    end
                    xlim([-0.2 0.2])
                    
                    prepHTBurstToPSTH = zeros(length(spikeIndicesHTBBurstStart),8);
                    time = linspace(-0.5,0.5,1001);
                    binaryPrepHTBurstPSTH = zeros(length(spikeIndicesHTBBurstStart),size(time,2));
                    for spikei = 1:length(spikeIndicesHTBBurstStart)
                        if spikeIndicesHTBBurstStart(spikei) > 5 && spikeIndicesHTBBurstStart(spikei) < size(spikeTimes,1) - 2
                            prepHTBurstToPSTH(spikei,:) = ((spikeTimes(spikeIndicesHTBBurstStart(spikei)-5:spikeIndicesHTBBurstStart(spikei)+2)) - spikeTimes(spikeIndicesHTBBurstStart(spikei)))';
                            for spikeii = 1:8
                                [~,idx(spikeii)] = min(abs(time-prepHTBurstToPSTH(spikei,spikeii)));
                                binaryPrepHTBurstPSTH(spikei,idx) = 1;
                            end
                        end
                    end
                    subplot(4,3,12)
                    plot(time,sum(binaryPrepHTBurstPSTH,1))
                    smoothHTBPSTH = smooth(sum(binaryPrepHTBurstPSTH,1),11,'lowess');
                    hold on; plot(time,smoothHTBPSTH,'g','LineWidth',2)
                    ylim([-1 1000]); xlim([-0.2 0.2])
                    
                    plotFileName = sprintf('%s/%s-sessionInd%d-Uniti-BurstingLTSHTBTvV-v%d.png', outputDir, rawDataLocs{sessioni}, uniti, v);
                    fprintf('\tSaving figure to file %s...\n', plotFileName);
                    saveas(FigH, plotFileName);
                    close all
                end
                cd(homeDir)
            end
        else
            cd(homeDir)
            continue
        end
    else
        cd(homeDir)
        continue
    end
end





