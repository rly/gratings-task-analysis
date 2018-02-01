function [Amatrix,Binvector,Unit_order,As_order]=assembly_assignment_matrix(As_across_bins,nneu,BinSizes,display,plotTitle,unitIDChars)
% Display the Assembly assignment matrix.
% USAGE: [Amatrix,Binvector,Unit_order]=assembly_assignment_matrix(As_across_bins,nneu,BinSizes,display)
%
% ARGUMENTS:
%   As_across_bins  - structure containing assembly information
%   nneu            - # of units
%   BinSizes          - list of investigated bin sizes
%   display         - Possible display options:
%           1) display='raw'  ->  no rearrangement made, assemblies and units displayed as 
%                                 from the detection algorithm
%           2) display='ordunit'  ->  Assembly order is left unchanged, units rearranged 
%                                     in order to separate assembly units from units not taking part to any assembly.
%           3) display='clustered'  -> Assembly and unit order changed in order to
%                                      display close to one other assemblies with similar elements
%
% RETURNS:
%       Amatrix     - ASSEMBLY ASSIGNMENT MATRIX: each column is an assembly 
%       each raw a unit. If unit i takes part in assembly j, Amatrix(i,j) displays
%       the firing delay (lag) of i with respect to the first assembly unit.
%       Binvector   - bin size. Binvector(j) is the bin size at which
%       Amatrix(:,j) has been detected.
%       Unit_order  - unit order in Amatrix
%
%
%
%
%  © 2016 Russo, Durstewitz.
%  for information please contact eleonora.russo@zi-mannheim.de; daniel.durstewitz@zi-mannheim.de.
%
%  last update 11/01/2016
%%
if nargin<4 , display='raw'; end;  
if nargin == 4
    if ~strcmp(display,'raw') && ~strcmp(display,'ordunit') && ~strcmp(display,'clustered')
        fprintf('''raw'', ''ordunit'', ''clustered'' are the only acceptable style specifications\n')
        return
    end
end




nAss_final=size(As_across_bins,2);
AAT=nan(nneu,nAss_final);
Binvector=nan(1,nAss_final);

for i=1:nAss_final
    AAT(As_across_bins{i}.elements,i)=As_across_bins{i}.lag;
    Binvector(i)=(As_across_bins{i}.bin);
end

switch display
            
    case 'raw'
        Unit_order=1:nneu;
        As_order=1:length(As_across_bins);
        
    case 'ordunit'
        aus=all(isnan(AAT),2);
        idx_nan=find(aus==1);
        idx_activ=find(aus==0);       
%         AAT=[AAT(~~sum(~isnan(AAT),2),:);AAT(~sum(~isnan(AAT),2),:)];   
%         Unit_order=[idx_activ;idx_nan];
%         As_order=1:length(As_across_bins);
        AAT=AAT(~~sum(~isnan(AAT),2),:); % remove nan cells
        Unit_order=idx_activ;
        As_order=1:length(As_across_bins);
        
    case 'clustered'
        aus=all(isnan(AAT),2);
        idx_nan=find(aus==1);
        idx_activ=find(aus==0);
        aus=~all(isnan(AAT),2);
        AAT_activ=AAT(aus,:);

        %%% order on the base of units co-occurrence 
        A01=~isnan(AAT_activ);
        M_assemb=zeros(size(AAT_activ,1));

        for n=1:size(A01,2)
            aus=find(A01(:,n)==1);
            M_assemb(aus,aus)=M_assemb(aus,aus)+1;
        end

        M_assemb(M_assemb==0)=0.0001;
        d_assemb=ones(size(AAT_activ,1))./M_assemb;
        d_assemb=d_assemb-diag(diag(d_assemb));
        Q=linkage(d_assemb,'average');
        [~,~,perm]=dendrogram(Q,0);
        perm1=idx_activ(perm);
        
        %%% order on the base of assemblies distance assemblies
        D = pdist(double(~isnan(AAT_activ))');
        Z = squareform(D);
        Q2=linkage(Z,'average');
        [~,~,perm2]=dendrogram(Q2,0);

        AAT=AAT_activ(perm,perm2);
        AAT=[AAT; nan(length(idx_nan),size(AAT,2))];
        Unit_order=[perm1;idx_nan];
        Binvector=Binvector(perm2);         % I apply the same order on the bin size vector
        As_order=perm2;

end

Amatrix=AAT;
binmat=log10(Binvector);
AAT(end+1,end+1)=0;
binmat(end+1,end+1)=0;




subplot(4,1,1:3)

h=pcolor(AAT);
set(h, 'EdgeColor', [0.8 0.8 0.8]);
caxis([0, max(AAT(:))])
colormap jet
ylabel('Unit #')
set(gca, 'XTick', []);
% RL commented temporarily
hcb=colorbar;
hcb.Label.String = 'Time lag \it l \rm (# bins)';
hcbPos = get(hcb, 'Position');
set(hcb,'Position',hcbPos + [hcbPos(3)*3 0 0 0]);
yy = 1:size(AAT,1);
set(gca,'YTick',yy(1)+0.5:yy(end),'Yticklabel',unitIDChars(Unit_order)) 
set(gcf, 'Color', [1,1,1]);
title(plotTitle);
        
        
        


ax2 = subplot(4,1,4);

pcolor(binmat)
yy = 1:size(AAT,2);
set(gca,'XTick',yy(1)+0.5:2:yy(end)+1+0.05,'Xticklabel',As_order(1:2:end)) 
xlabel('Assembly #')
set(gca, 'YTick', []);
currPos = get(gca, 'Position');
hC = colorbar;
graymap = colormap(ax2, 'gray'); % RL added
graymap = graymap(round(linspace(1, size(graymap, 1), length(BinSizes))),:);
colormap(ax2, graymap);
if length(BinSizes)>1
    caxis([log10(min(BinSizes)), log10(max(BinSizes))]);
else
    caxis([log10(BinSizes-0.001), log10(BinSizes+0.001)]);
end
set(hC,'Ytick',log10(BinSizes),'YTicklabel',BinSizes);
set(hC,'Location','southoutside')
set(ax2, 'Position', currPos + [0 0.1 0 -0.1]);
hcPos = get(hC, 'Position');
set(hC, 'Position', hcPos + [0 -0.05 0 0.03]);
hC.Label.String = 'Temporal precision \Delta (sec)';





end
