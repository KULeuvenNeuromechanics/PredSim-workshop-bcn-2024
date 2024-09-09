clear
close all
clc

trialInfo = readtable('timepoints_IPSA.xlsx');

figure
for i=1:6
    clear labelsEMG
    labelsEMG = {'LMEH' 'LBIF' 'LREF'};
    clear EMG
    EMG = read_motionFile_v40(['../EMG/filtEMG_IPSA/' trialInfo.EMGfile{i}]);
    if contains(trialInfo.movement{i},'slow')
        for l=1:length(labelsEMG)
            subplot(length(labelsEMG),2,(l-1)*2+1)
            plot(EMG.data(:,1),EMG.data(:,find(ismember(EMG.labels,labelsEMG{l}))))
            hold on
            if l==1
                title({'Slow hams IPSA Left',labelsEMG{l}})
            else
                title(labelsEMG{l})
            end
        end
    elseif contains(trialInfo.movement{i},'fast')
        for l=1:length(labelsEMG)
            subplot(length(labelsEMG),2,(l-1)*2+2)
            plot(EMG.data(:,1),EMG.data(:,find(ismember(EMG.labels,labelsEMG{l}))))
            hold on
            if l==1
                title({'Slow hams IPSA Left',labelsEMG{l}})
            else
                title(labelsEMG{l})
            end
        end
    end
end
        
figure
for i=7:14
    clear labelsEMG
    labelsEMG = {'RMEH' 'RBIF' 'RREF'};
    clear EMG
    EMG = read_motionFile_v40(['../EMG/filtEMG_IPSA/' trialInfo.EMGfile{i}]);
    if contains(trialInfo.movement{i},'slow')
        for l=1:length(labelsEMG)
            subplot(length(labelsEMG),2,(l-1)*2+1)
            plot(EMG.data(:,1),EMG.data(:,find(ismember(EMG.labels,labelsEMG{l}))))
            hold on
            if l==1
                title({'Slow hams IPSA Right',labelsEMG{l}})
            else
                title(labelsEMG{l})
            end
        end
    elseif contains(trialInfo.movement{i},'fast')
        for l=1:length(labelsEMG)
            subplot(length(labelsEMG),2,(l-1)*2+2)
            plot(EMG.data(:,1),EMG.data(:,find(ismember(EMG.labels,labelsEMG{l}))))
            hold on
            if l==1
                title({'Slow hams IPSA Right',labelsEMG{l}})
            else
                title(labelsEMG{l})
            end
        end
    end
end
        
figure
for i=15:22
    clear labelsEMG
    labelsEMG = {'LTIA' 'LSOL' 'LGAS' 'LLATGAS'};
    clear EMG
    EMG = read_motionFile_v40(['../EMG/filtEMG_IPSA/' trialInfo.EMGfile{i}]);
    if contains(trialInfo.movement{i},'slow')
        for l=1:length(labelsEMG)
            subplot(length(labelsEMG),2,(l-1)*2+1)
            plot(EMG.data(:,1),EMG.data(:,find(ismember(EMG.labels,labelsEMG{l}))))
            hold on
            if l==1
                title({'Slow pf IPSA Left',labelsEMG{l}})
            else
                title(labelsEMG{l})
            end
        end
    elseif contains(trialInfo.movement{i},'fast')
        for l=1:length(labelsEMG)
            subplot(length(labelsEMG),2,(l-1)*2+2)
            plot(EMG.data(:,1),EMG.data(:,find(ismember(EMG.labels,labelsEMG{l}))))
            hold on
            if l==1
                title({'Slow pf IPSA Left',labelsEMG{l}})
            else
                title(labelsEMG{l})
            end
        end
    end
end
    
figure
for i=23:30
    clear labelsEMG
    labelsEMG = {'RTIA' 'RSOL' 'RGAS' 'RLATGAS'};
    clear EMG
    EMG = read_motionFile_v40(['../EMG/filtEMG_IPSA/' trialInfo.EMGfile{i}]);
    if contains(trialInfo.movement{i},'slow')
        for l=1:length(labelsEMG)
            subplot(length(labelsEMG),2,(l-1)*2+1)
            plot(EMG.data(:,1),EMG.data(:,find(ismember(EMG.labels,labelsEMG{l}))))
            hold on
            if l==1
                title({'Slow pf IPSA Right',labelsEMG{l}})
            else
                title(labelsEMG{l})
            end
        end
    elseif contains(trialInfo.movement{i},'fast')
        for l=1:length(labelsEMG)
            subplot(length(labelsEMG),2,(l-1)*2+2)
            plot(EMG.data(:,1),EMG.data(:,find(ismember(EMG.labels,labelsEMG{l}))))
            hold on
            if l==1
                title({'Slow pf IPSA Right',labelsEMG{l}})
            else
                title(labelsEMG{l})
            end
        end
    end
end
    