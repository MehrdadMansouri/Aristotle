function G2S = Gene2SNP(DataDirectory)

Files = dir([DataDirectory '\SNP to Gene']);
Files = Files(~ismember({Files.name},{'.','..'}));
Files = struct2cell(Files); Files = Files(1,:)';

G2S = [];
for i = 1:length(Files)
    S2Gi = readtable([DataDirectory '\SNP to Gene\' Files{i}]);
    S2Gi = S2Gi(:,[6 1]);
    S2Gi.Properties.VariableNames = {'G','S'};
    G2S = [G2S;S2Gi];
end
G2S = unique(G2S);
G2S = table2cell(G2S);
Temp = strcmp('-',G2S(:,1));
G2S(Temp,:) = [];
