clc; clear variables; close all; warning off; disp('Initiated')
n0 = 6;
%% Gene to SNPs
G2S = [];
for i = 1:15
    S2Gi = readtable(['N' num2str(i) '.txt']);
    S2Gi = S2Gi(:,[6 1]);
    S2Gi.Properties.VariableNames = {'G','S'};
    G2S = [G2S;S2Gi];
end
G2S = unique(G2S);
G2S = table2cell(G2S);
Temp = strcmp('-',G2S(:,1));
G2S(Temp,:) = [];
disp('Gene to SNPs Done')

%% Pathway to SNPs
Pathways = dir('Coclustering/Pathways');
Pathways = Pathways(~ismember({Pathways.name},{'.','..'}));
Pathways = struct2cell(Pathways); Pathways = Pathways(1,:)';

P2S = cell(length(Pathways),1);
for i = 1:length(Pathways)
    P2Gi = readtable(['Coclustering\Pathways\' Pathways{i}]);
    P2Gi = P2Gi(:,1);
    P2Gi = table2cell(P2Gi);
    P2Si = cell(0);
    for j = 1:length(P2Gi)
        P2Si = [P2Si;G2S(strcmp(P2Gi(j),G2S(:,1)),2)];
    end
    P2S{i} = P2Si;
end
disp('Pathway to SNPs Done')

%% Storing SNPs

disp('Storing SNPs Done')

