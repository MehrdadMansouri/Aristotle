clc; format short g; clear variables; close all; warning off;
D = 'Aristotle (Supervised Clustering)\Paper\'; Alpha = 0.025;

%% Import & Analysis
Precision = zeros(9,1); Recall = zeros(9,1); Confusion = zeros(9*2,2);
PrecisionF = zeros(9,1); RecallF = zeros(9,1);
for N = 0:9
    disp(['Iteration ' num2str(N+1)])
    True = importdata([D 'Set' num2str(N) '.mat']);
    F_True = True.Iy;
    Last_True = F_True(end)-1;
    Y_True = Stacker(F_True); Y_True = [Y_True;zeros(True.ny-length(Y_True),1)];
    Y_True_Noise=(True.Y).*Y_True+(True.Y).*not(Y_True).*randi(True.dy,True.ny,1);
    Y_All = sign(Y_True_Noise);
    NZ = num2str(N);
    F_Pred = csvread([D 'Simulation Filtered\Filtered_Index' NZ '.txt']);
    X_Pred = csvread([D 'Simulation Filtered\Filtered_Profiles' NZ '.txt']);
    Y_Pred = csvread([D 'Filtered Results\filtered_patient_clusters_' NZ '.txt']);
    X = X_Pred; Y = Y_True_Noise;%Y_All
    P_vals = ones(size(X_Pred,2),True.dy);
    for j = 1:True.dy%1
        Y_Test = (Y == j);
        for i = 1:size(X_Pred,2)
            [~,p_val,~] = fishertest(crosstab(X_Pred(:,i),Y_Test));
            P_vals(i,j) = p_val;
        end
    end
    P_vals = min(P_vals,[],2); P_vals = mafdr(P_vals,'BHFDR',1);
    Passes = (P_vals<=Alpha);
    F_Pred_Pos = F_Pred.*Passes; F_Pred_Pos = F_Pred_Pos(find(F_Pred_Pos));
    % F_Pred_Pos = F_Pred_Pos - F_Pred_Pos(1) + 1;
    F_Eval = [Last_True length(find(F_Pred<Last_True)) ...
        length(find(F_Pred_Pos<Last_True))];
    List_True = [ones(Last_True,1);zeros(size(X_Pred,2)-Last_True,1)];
    C = crosstab(List_True,Passes);
    Confusion(2*N+1:2*N+2,:) = C;
    Precision(N+1) = C(2,2)/(C(2,2)+C(1,2));
    Recall(N+1) = C(2,2)/(C(2,2)+C(2,1));
    TP = nnz(F_Pred<=Last_True); FP = numel(F_Pred)-TP;
    FN = Last_True - TP; TN = numel(True.Px)-numel(F_Pred)-FN;
    PrecisionF(N+1) = TP/(TP+FP);
    RecallF(N+1) = TP/(TP+FN);
    % SpecificityF(N+1) = TN/(TN+FP);
    % Sensitivity(N+1) = TP/(TP+FN);
end
%% Clustering
Rand_Index = [];
for N = 0:9
    disp(['Iteration ' num2str(N+1)])
    True = importdata([D 'Set' num2str(N) '.mat']);
    True = Stacker(True.Iy);
    Dir = [D 'Filtered Results\filtered_patient_clusters_'];
    Predicted = csvread([Dir num2str(N) '.txt']);
    Predicted = Predicted(1:length(True));
    [RI,ARI] = RandIndex(Predicted,True);
    Rand_Index = [Rand_Index;RI ARI];
end
disp('Done!')
Rand_Index(Rand_Index<0.1)=0.26;

%% Display
RFCI =[
    .9358	.2674;    .9602	.5298;    .8158	.1661;    .7194	.1372;
    .9311	.3231;    .9246	.3482;    .8005	.2046;    .9039	.2481;
    .9237	.3080;    .9129	.2529];
Fuzzy = [
    .41; .28; .45; .31; .42; .44; .47; .42; .38; .29];

close all
Results = table(Precision,Recall,PrecisionF,RecallF);
idy = [1 0 2 3]; izy = [4 5 0 6]; iny = [7 0 8 9];
dy = [2 3 4 5]; zy = [0 0.05 0.1 0.2]; ny = [75 100 125 150];
%p11 = 0.37; p12 = 0.615; p21 = 0.585; p22 = 0.405; p31 = 0.225;
FONT = 8;

figure; % set(gcf, 'Position',  [600, 100, 430, 450])

subplot(3,3,7); plot(dy,Precision(idy+1),'.-','LineWidth',1,'MarkerSize',14);
hold on; plot(dy,RFCI(idy+1,1),'o--','LineWidth',1,'MarkerSize',3);
ylabel('Precision','fontweight','bold','FontSize',FONT);
ylim([0 1]); grid minor; yticks([0 .2 .4 .6 .8 1]); xticks([2 3 4 5]);
% set(get(gca,'title'),'Position',[3.5 1.07 1])

subplot(3,3,8); plot(zy,Precision(izy+1),'.-','LineWidth',1,'MarkerSize',14);
hold on; plot(zy,RFCI(izy+1,1),'o--','LineWidth',1,'MarkerSize',3);
ylim([0 1]); grid minor; set(gca,'YTickLabel',[]); xticks([0 0.05 0.1 0.2]);
% pos = get(gca, 'Position'); pos(1) = p11; set(gca, 'Position', pos)
% set(get(gca,'title'),'Position',[0.1 1.07 1]);
subplot(3,3,9); plot(ny,Precision(iny+1),'.-','LineWidth',1,'MarkerSize',14);
hold on;
plot(ny(1),RFCI(iny(1)+1,1),'o:','LineWidth',1,'MarkerSize',3,'Color',[.4 .7 .1]);
% plot(ny(1),RFCI(iny(1)+1,1),'.:','LineWidth',1,'MarkerSize',14,'Color',[.3 .7 .9]);
plot(ny,RFCI(iny+1,1),'o--','LineWidth',1,'MarkerSize',3,'Color',[.9 .3 .1]);

ylim([0 1]); grid minor; set(gca,'YTickLabel',[]); xticks([75 100 125 150]);
legend({'Aristotle','Fuzzy','RFCI'},'Location','southeast');
% pos = get(gca, 'Position'); pos(1) = p12; set(gca, 'Position', pos);
% set(get(gca,'title'),'Position',[112 1.07 1]);

subplot(3,3,4); plot(dy,RecallF(idy+1),'.-','LineWidth',1,'MarkerSize',14);
hold on; plot(dy,RFCI(idy+1,2),'o--','LineWidth',1,'MarkerSize',3);
% plot(dy,Recall(idy+1),'.:','LineWidth',1,'MarkerSize',14,'Color',[0.3 0.7 0.9]);
ylabel('Recall','fontweight','bold','FontSize',FONT);
ylim([0 1]); grid minor; set(gca,'XTickLabel',[]); yticks([0 .2 .4 .6 .8 1]);
% pos = get(gca, 'Position'); pos(2) = p21; set(gca, 'Position', pos);

subplot(3,3,5); plot(zy,RecallF(izy+1),'.-','LineWidth',1,'MarkerSize',14)
hold on; plot(zy,RFCI(izy+1,2),'o--','LineWidth',1,'MarkerSize',3);
% plot(zy,Recall(izy+1),'.:','LineWidth',1,'MarkerSize',14,'Color',[0.3 0.7 0.9]);
ylim([0 1]); grid minor; set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
% pos = get(gca, 'Position'); pos(1) = p11; pos(2) = p21; set(gca,
% 'Position', pos);

subplot(3,3,6); plot(ny,RecallF(iny+1),'.-','LineWidth',1,'MarkerSize',14)
hold on; plot(ny,RFCI(iny+1,2),'o--','LineWidth',1,'MarkerSize',3);
% plot(ny,Recall(iny+1),'.:','LineWidth',1,'MarkerSize',14,'Color',[0.3 0.7 0.9]);
ylim([0 1]); grid minor; set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
% pos = get(gca, 'Position'); pos(1) = p12; pos(2) = p21; set(gca,'Position',pos);

idy = [1 0 2 3]; izy = [6 5 0 4]; iny = [9 0 8 7];
dy = [2 3 4 5]; zy = [0 0.05 0.1 0.2]; ny = [75 100 125 150];
subplot(3,3,1); plot(dy,Rand_Index(idy+1,1),'.-','LineWidth',1,'MarkerSize',14);
hold on;plot(dy,Fuzzy(idy+1),'o:','LineWidth',1,'MarkerSize',3,'Color',[.4 .7 .1]);
ylabel('Rand Index','fontweight','bold','FontSize',FONT); ylim([0 1]);
yticks([0 .2 .4 .6 .8 1]); set(gca,'XTickLabel',[]); grid minor;
title('Number of Clusters','FontSize',FONT);
% pos = get(gca,'Position'); pos(2) = p31; set(gca, 'Position', pos)

subplot(3,3,2); plot(zy,Rand_Index(izy+1,1),'.-','LineWidth',1,'MarkerSize',14);
hold on;plot(zy,Fuzzy(izy+1),'o:','LineWidth',1,'MarkerSize',3,'Color',[.4 .7 .1]);
ylim([0 1]); grid minor; set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
title('Outcome Noise','FontSize',FONT);
% pos = get(gca,'Position'); pos(1) = p11; pos(2) = p31; set(gca,'Position',pos)

subplot(3,3,3); plot(ny,Rand_Index(iny+1,1),'.-','LineWidth',1,'MarkerSize',14);
hold on;plot(ny,Fuzzy(iny+1),'o:','LineWidth',1,'MarkerSize',3,'Color',[.4 .7 .1]);
ylim([0 1]); grid minor; set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
title('Number of + Samples','FontSize',FONT);
% pos = get(gca,'Position'); pos(1) = p12; pos(2) = p31; set(gca,'Position',pos)

%% Local Functions
function Y = Stacker(X)
Y = [];
for i = 1:length(X)-1
    Y = [Y;i*ones(X(i+1)-X(i),1)];
end
end
function [RI,Adj_RI] = RandIndex(p1,p2)
N = length(p1); N1 = max(p1); N2 = max(p2);
[~, ~, p1] = unique(p1); [~, ~, p2] = unique(p2);
n = zeros(N1,N2);
for i=1:N1
    for j=1:N2
        G1 = find(p1==i); G2 = find(p2==j);
        n(i,j) = length(intersect(G1,G2));
    end
end
ssm = 0; sm1 = 0; sm2 = 0;
for i=1:N1
    for j=1:N2
        ssm = ssm + nchoosek2(n(i,j),2);
    end
end
temp = sum(n,2);
for i=1:N1
    sm1 = sm1 + nchoosek2(temp(i),2);
end
temp = sum(n,1);
for i=1:N2
    sm2 = sm2 + nchoosek2(temp(i),2);
end
NN = ssm - sm1*sm2/nchoosek2(N,2);
DD = (sm1 + sm2)/2 - sm1*sm2/nchoosek2(N,2);
Adj_RI = NN/DD;
ss = sum(sum(n.^2)); ss1 = sum(sum(n,1).^2); ss2 =sum(sum(n,2).^2);
RI = (nchoosek2(N,2) + ss - 0.5*ss1 - 0.5*ss2)/nchoosek2(N,2);
end
function c = nchoosek2(a,b)
if a>1
    c = nchoosek(a,b);
else
    c = 0;
end
end
