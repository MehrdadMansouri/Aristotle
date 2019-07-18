clc; close all;
q = 0.08;
P = P_values; m = length(P); j = 0; k = 0; SS = zeros(m,1); m1 = 100000;
for i = 1:m
    BH(i,1) = i*q/m; BH1(i,1) = i*q/m1;
    S(i)  = (1-P(i))/(m+1-i); % S(i)  = round(S(i),4,'significant');
    if i > 1 && j == 0
        SS(i) = S(i) - S(i-1);
        if SS(i) < 0
            j = i;
            m0 = 1 + 1/S(i);
        end
    end
end
for i = m:-1:1
    BH0(i,1) = i*q/m0;
    if P(i)<= i*q/m0 && k == 0
        k = i;
    end
end
[j m0 S(j) k P(k)]
figure; loglog([P BH BH1 BH0]);
xlim([1 m]); grid minor; xlabel('SNP'); ylabel('P Value');
title(['Aristotle Multiple Hypothesis Testing'])
legend('P Values','Candidates','All','Estimated')
