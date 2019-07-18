function P2S = Pathway2SNP(DataDirectory,G2S)

Files = dir([DataDirectory '\Pathway to Gene']);
Files = Files(~ismember({Files.name},{'.','..'}));
Files = struct2cell(Files); Files = Files(1,:)';

P2S = cell(length(Files),1);
for i = 1:length(Files)
    P2Gi = readtable([DataDirectory '\Pathway to Gene\' Files{i}]);
    P2Gi = P2Gi(:,1);
    P2Gi = table2cell(P2Gi);
    P2Si = cell(0);
    for j = 1:length(P2Gi)
        P2Si = [P2Si;G2S(strcmp(P2Gi(j),G2S(:,1)),2)];
    end
    P2S{i} = P2Si;
end
