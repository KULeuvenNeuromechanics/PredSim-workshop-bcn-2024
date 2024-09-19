clear
close all
clc
allsegments_KE_L = load('CP4_T0_MEH_L.mat','allsegments');
allsegments_PF_L = load('CP4_T0_PF_L.mat','allsegments');
allsegments_KE_R = load('CP4_T0_MEH_R.mat','allsegments');
allsegments_PF_R = load('CP4_T0_PF_R.mat','allsegments');

f1 = fields(allsegments_KE_L.allsegments);
f2 = fields(allsegments_PF_L.allsegments);
f3 = fields(allsegments_PF_R.allsegments);
f4 = fields(allsegments_KE_R.allsegments);

equalKE = all(ismember(f1,f4));
ismember(f2,f3);
ismember(f3,f2);
ismember(f2,f1) == ismember(f3,f1);
% emg4 is not there in allsegments_KE_L. Add it
for f = 1:length(f3)
    if ~ismember(f3{f},f1)
        for i=1:length(allsegments_KE_L.allsegments)
            allsegments_KE_L.allsegments(i).(f3{f}) = [];
            if equalKE
                allsegments_KE_R.allsegments(i).(f3{f}) = [];
            end
        end
    end
end

% Always put KE first, then KF, then DF
allsegments = [allsegments_KE_L.allsegments allsegments_KE_R.allsegments allsegments_PF_L.allsegments allsegments_PF_R.allsegments];
save('IPSA_data.mat','allsegments')
