clc; format short g; clear variables; close all; warning off;
L = 4; P0 = 0.05; n0 = 8; Dir = 'Aristotle (Supervised Clustering)\Data';

Files = dir([Dir '\Weights']);
Files = Files(~ismember({Files.name},{'.','..'}));
Files = struct2cell(Files); Files = Files(1,:)';
[M,~] = xlsread([Dir '\ClinicalData.xlsx']);
Temp=M(:,5); M(:,5)=(Temp-min(Temp))/(max(Temp)-min(Temp));
Temp=M(:,6); M(:,6)=(Temp-min(Temp))/(max(Temp)-min(Temp));
Z = M(:,4:10); dP = pdist2(Z,Z,'cityblock');
PC = csvread([Dir '\Patient Clusters\Consensus_kp.txt']); R1 = [1 5 8];
% R = M(:,3); D = [ones(434,1) M(:,1:2)];

Candidates = table; NumHyp = 0;
for i = 1:length(Files)
    disp(['iteration ' num2str(i)])
    w = csvread([Dir '\Weights\' Files{i}]);
    F = isoutlier(w,'median','ThresholdFactor',L); F = F.*[w>median(w)];
    PGS = csvread([Dir '\Pathway Patients\' Files{i}(1:end-6) '.csv']);
    [~,NameGS] = xlsread([Dir '\Pathway SNPs\' Files{i}(1:end-15) '_SNPs.xlsx']);
    PGF = PGS(:,find(F)'); NameGF = NameGS(find(F)');
    P = zeros(length(NameGF),1);
    for j = 1:length(NameGF)
        G = PGF(:,j); U = find(G); V = find(1-G); dUV = dP(U,V);
        [P1,P2] = ind2sub(size(dUV),find(munkres(dUV)));
        if isempty(P1)==0
            Temp3 = P1; Temp2 = U'; Temp1 = 1:length(U);
            for k = 1:numel(Temp2); Temp3(P1 == Temp1(k)) = Temp2(k); end
            P1 = Temp3;
            Temp3 = P2; Temp2 = V'; Temp1 = 1:length(V);
            for k = 1:numel(Temp2); Temp3(P2 == Temp1(k)) = Temp2(k); end
            P2 = Temp3;
            pk = ones(length(R1),1);
            for kk = R1
                kR = (PC==kk);
                NumHyp = NumHyp + (sum(kR)>n0);
                T = zeros(2,2);
                for k = 1:length(P1)
                    T = T + [kR(P1(k)); 1-kR(P1(k))] * [kR(P2(k)) 1-kR(P2(k))];
                end
                T1 = T(1,2); T2 = T(2,1);
                pk(kk) = 1 - chi2cdf((T1 - T2)^2 / (T1 + T2),1);
            end
            pk(pk==0) = 1;
            P(j) = min(pk);
        end
    end
    % SNP = NameGF(P<P0)'; P_Value = P(P<P0);
    % Pathway = str2double(Files{i}(8:end-15))*ones(nnz(P<P0),1);
    SNP = NameGF'; P_Value = P;
    Pathway = str2double(Files{i}(8:end-15))*ones(size(P));
    Candidates = [Candidates; table(SNP,P_Value,Pathway)];
end
disp('Analyzed')
Candidates = sortrows(Candidates,2);
[~,New,~] = unique(Candidates(:,1),'stable');
Predictions = table;
for i = 1:length(New)-1
    Same = Candidates(New(i):New(i+1)-1,:);
    PathsofSame = mat2str(table2array(Same(:,3)));
    PathsofSame = erase(PathsofSame,["[","]"]);
    PathsofSame = cellstr(string(PathsofSame));
    Predictions = [Predictions; Candidates(New(i),1:2) PathsofSame];
end
P_values = table2array(Predictions(:,2));
n = length(P_values);
loglog(1:n, (P0/NumHyp)*[1:n], 1:n, P_values, '.-');
xlim([1 n]); grid minor; xlabel('SNP'); ylabel('P Value');
title(['Aristotle 1.1 with L = ' num2str(L)])
legend('BH Baseline','Hypotheses','Location','northwest')

writetable(Predictions,'AristotleFinalCausesFiltered.xlsx')


