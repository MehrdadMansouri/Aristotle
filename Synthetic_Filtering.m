clc; format short g; clear variables; close all; warning off;
Threshold = 3;
for Sim_Num = 0:9
    disp(['Simulation Number ' num2str(Sim_Num)])
    Dir = ['Aristotle Tests\Filtered\' num2str(Sim_Num)];
    load(['Aristotle Tests\Set' num2str(Sim_Num) '.mat'])
    Files = dir(Dir);
    Files = Files(~ismember({Files.name},{'.','..'}));
    Files = struct2cell(Files); Files = Files(1,:)';
    Files = Files(startsWith(Files,'feature_weights'));
    Files = natsort(Files);
    HH = [];
    for i = 1:length(Files)
        w = csvread([Dir '\' Files{i}]);
        F = isoutlier(w,'median','ThresholdFactor',Threshold);
        F = F(w>median(w)); F = find(F);
        G = find(Pathway==i);
        H = G(F);
        HH = [HH;H];
    end
    HH = sort(HH);
    XX = X(:,HH);
    csvwrite(['Filtered_Profiles' num2str(Sim_Num) '.txt'],XX)
    csvwrite(['Filtered_Index' num2str(Sim_Num) '.txt'],HH)
end
disp('Done !')
