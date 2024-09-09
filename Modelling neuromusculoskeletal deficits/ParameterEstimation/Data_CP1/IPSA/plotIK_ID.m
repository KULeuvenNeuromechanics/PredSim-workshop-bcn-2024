clear
close all
clc

trialInfo = readtable('timepoints_IPSA.xlsx');

figure
for i=1:6
    clear IK
    clear ID
    IK = read_motionFile_v40(['./JointAngles/' trialInfo.IKfile{i}]);
    ID = read_motionFile_v40(['./ID/' trialInfo.IDfile{i}]);
    if contains(trialInfo.movement{i},'slow')
        subplot(2,2,1)
        plot(IK.data(:,1),IK.data(:,find(ismember(IK.labels,trialInfo.ipsaDOF{i}))))
        hold on
        title({'Slow hams IPSA','knee angle left'})
        subplot(2,2,3)
        plot(ID.data(:,1),ID.data(:,find(ismember(ID.labels,[trialInfo.ipsaDOF{i} '_moment']))))
        hold on
        title('knee moment left')
    elseif contains(trialInfo.movement{i},'fast')
        subplot(2,2,2)
        plot(IK.data(:,1),IK.data(:,find(ismember(IK.labels,trialInfo.ipsaDOF{i}))))
        hold on
        title({'Fast hams IPSA','knee angle left'})
        subplot(2,2,4)
        plot(ID.data(:,1),ID.data(:,find(ismember(ID.labels,[trialInfo.ipsaDOF{i} '_moment']))))
        hold on
        title('knee moment left')
    end
end
        
figure
for i=7:14
    clear IK
    clear ID
    IK = read_motionFile_v40(['./JointAngles/' trialInfo.IKfile{i}]);
    ID = read_motionFile_v40(['./ID/' trialInfo.IDfile{i}]);
    if contains(trialInfo.movement{i},'slow')
        subplot(2,2,1)
        plot(IK.data(:,1),IK.data(:,find(ismember(IK.labels,trialInfo.ipsaDOF{i}))))
        hold on
        title({'Slow hams IPSA','knee angle right'})
        subplot(2,2,3)
        plot(ID.data(:,1),ID.data(:,find(ismember(ID.labels,[trialInfo.ipsaDOF{i} '_moment']))))
        hold on
        title('knee moment right')
    elseif contains(trialInfo.movement{i},'fast')
        subplot(2,2,2)
        plot(IK.data(:,1),IK.data(:,find(ismember(IK.labels,trialInfo.ipsaDOF{i}))))
        hold on
        title({'Fast hams IPSA','knee angle right'})
        subplot(2,2,4)
        plot(ID.data(:,1),ID.data(:,find(ismember(ID.labels,[trialInfo.ipsaDOF{i} '_moment']))))
        hold on
        title('knee moment right')
    end
end

figure
for i=15:22
    clear IK
    clear ID
    IK = read_motionFile_v40(['./JointAngles/' trialInfo.IKfile{i}]);
    ID = read_motionFile_v40(['./ID/' trialInfo.IDfile{i}]);
    if contains(trialInfo.movement{i},'slow')
        subplot(2,2,1)
        plot(IK.data(:,1),IK.data(:,find(ismember(IK.labels,trialInfo.ipsaDOF{i}))))
        hold on
        title({'Slow pf IPSA','ankle angle left'})
        subplot(2,2,3)
        plot(ID.data(:,1),ID.data(:,find(ismember(ID.labels,[trialInfo.ipsaDOF{i} '_moment']))))
        hold on
        title('ankle moment left')
    elseif contains(trialInfo.movement{i},'fast')
        subplot(2,2,2)
        plot(IK.data(:,1),IK.data(:,find(ismember(IK.labels,trialInfo.ipsaDOF{i}))))
        hold on
        title({'Fast pf IPSA','ankle angle left'})
        subplot(2,2,4)
        plot(ID.data(:,1),ID.data(:,find(ismember(ID.labels,[trialInfo.ipsaDOF{i} '_moment']))))
        hold on
        title('ankle moment left')
    end
end
    
figure
for i=23:30
    clear IK
    clear ID
    IK = read_motionFile_v40(['./JointAngles/' trialInfo.IKfile{i}]);
    ID = read_motionFile_v40(['./ID/' trialInfo.IDfile{i}]);
    if contains(trialInfo.movement{i},'slow')
        subplot(2,2,1)
        plot(IK.data(:,1),IK.data(:,find(ismember(IK.labels,trialInfo.ipsaDOF{i}))))
        hold on
        title({'Slow pf IPSA','ankle angle right'})
        subplot(2,2,3)
        plot(ID.data(:,1),ID.data(:,find(ismember(ID.labels,[trialInfo.ipsaDOF{i} '_moment']))))
        hold on
        title('ankle moment right')
    elseif contains(trialInfo.movement{i},'fast')
        subplot(2,2,2)
        plot(IK.data(:,1),IK.data(:,find(ismember(IK.labels,trialInfo.ipsaDOF{i}))))
        hold on
        title({'Fast pf IPSA','ankle angle right'})
        subplot(2,2,4)
        plot(ID.data(:,1),ID.data(:,find(ismember(ID.labels,[trialInfo.ipsaDOF{i} '_moment']))))
        hold on
        title('ankle moment right')
    end
end
    