clc; format short g; clear variables; close all; warning off;
Alpha = 4;
Dir = 'Aristotle (Supervised Clustering)\Data';
Files = dir([Dir '\Weights']);
Files = Files(~ismember({Files.name},{'.','..'}));
Files = struct2cell(Files); Files = Files(1,:)';
PG = []; SNP = []; Weight = []; W = []; Pathway = []; Number = [];
for i = 1:length(Files)
    disp(['iteration ' num2str(i)])
    Wi = csvread([Dir '\Weights\' Files{i}]);
    PGi = csvread([Dir '\Pathway Patients\' Files{i}(1:end-6) '.csv']);
    [~,NGi] = xlsread([Dir '\Pathway SNPs\' Files{i}(1:end-15) '_SNPs.xlsx']);
    F = isoutlier(Wi,'median','ThresholdFactor',Alpha); F = F(Wi>median(Wi));
    F = find(F); Number = [Number;F];
    PGi = PGi(:,F); PG = [PG PGi];
    NGi = NGi(F); SNP = [SNP;NGi'];
    Weight = [Weight;Wi(F)]; W = [W;Wi];
    QGi = i*ones(size(F)); Pathway = [Pathway;QGi];
end
Table = table(SNP,Weight,Pathway,Number);
[~,Temp,~] = unique(Table.SNP,'stable');
PG1 = PG(:,Temp);
csvwrite('Profile.txt',PG)
writetable(Table,'Feature_Selections.txt')
