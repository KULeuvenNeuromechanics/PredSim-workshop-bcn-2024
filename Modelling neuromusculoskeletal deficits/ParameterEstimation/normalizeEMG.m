function DatStore = normalizeEMG(Misc,DatStore,Results)
% --------------------------------------------------------------------------
%normalizeEMG
%    This function normalizes EMG data and updates the saved results
% 
% INPUT:
%     Misc
%     Miscellaneous info used through the code
% 
%     DatStore
%     Structure of all data
%     
%     Results
%     Structure of all results
%     
% Original author: Dhruv Gupta
% Original date: Feb 07, 2023
%
% Last edit by: Dhruv Gupta
% Last edit date: Feb 07, 2023
% --------------------------------------------------------------------------
%% Normalize EMG
% load(fullfile(Misc.OutPath,'emgMax_all.mat'))
load(fullfile(Misc.DataPath,'EMG','emgMax_all.mat'))
% Misc.trials_sel = Misc.trials_gait;
if Misc.boolEMG
    maxMRS = nan(Misc.nAllMuscList,1);
    for t = Misc.trials_sel
        for m=1:Misc.NMuscles(t)
            idx_m = Misc.idx_allMuscleList{t}(m);
            if Misc.noSOMRS
                maxMRS(idx_m,1) = 1;
            else
                if Misc.normalizeToMRS
                    maxMRS(idx_m,1) = max(maxMRS(idx_m,1),max(Results.MActivation(t).genericMRS(m,:)));
                else
                    maxMRS(idx_m,1) = 1;
                end
            end
        end
    end
end
% save(fullfile(Misc.OutPath,'normMax_all.mat'),'maxEMG','maxMRS')

for t = Misc.trials_sel
    for m=1:DatStore(t).EMG.nEMG
        idx_m = Misc.idx_EMGsel{t}(m,1);
        DatStore(t).EMG.EMGsel(:,m) = DatStore(t).EMG.EMGsel(:,m)./(maxEMG(idx_m,1)/maxMRS(idx_m,1));
    end
end
% for t = Misc.trials_sel
%     DatStore(t).EMG.time_delay100 = DatStore(t).EMG.time + 0.1;
%     for m=1:DatStore(t).EMG.nEMG
%         DatStore(t).EMG.EMGsel_delay100(:,m) = interp1(DatStore(t).EMG.time,DatStore(t).EMG.EMGsel(:,m),DatStore(t).EMG.time_delay100);
%     end
% end
