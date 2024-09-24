function [varargout] = run_track_sim(S,osim_path)
% --------------------------------------------------------------------------
% run_track_sim
%  This functions calls the external function, runs the optimal control 
%  problem and subfunction for each step in the simulation and postrocess
%  the data.

% Make sure casadi path is set up correctly. This needs to happen before
% adding a simulation to the batch
if ~isfield(S.solver,'CasADi_path')
    try
        S.solver.CasADi_path = casadi.GlobalOptions.getCasadiPath();
    catch
        error("Please add CasADi to the matlab search path, or pass the path " + ...
            "to your CasADi installation (top folder) to S.solver.CasADi_path.")
    end
elseif ~isempty(S.solver.CasADi_path) && ~isfolder(S.solver.CasADi_path)
    error("Unable to find the path assigned to S.solver.CasADi_path:" + ...
        " \n\t%s",S.solver.CasADi_path)
end

% Make sure folder to save results exists
OutFolder = S.subject.save_folder;
if ~isfolder(OutFolder)
    mkdir(OutFolder);
end

addpath([S.misc.main_path '\VariousFunctions']);
if isempty(S.post_process.result_filename)     
    cond = 1;
    ct = 1;
    while cond
        result_filename = [S.subject.name '_v' num2str(ct)];
        if ~isfile(fullfile(OutFolder,[result_filename '.mat']))
            cond = 0;
        end
        ct = ct+1;
    end
end

if nargout == 1
    varargout{1} = S.post_process.result_filename;
end

%% Start diary
t00 = tic;
Outname = fullfile(S.subject.save_folder,[S.post_process.result_filename '_log.txt']);
diary(Outname);
disp(' ')
disp(['Subject name: ' S.subject.name])
disp(' ')
disp(' ')

model_info.nq=6;
state_names=S.expdata.IKdata.colheaders(2:2:end);
for i=1:length(state_names)
    state_names_split=strsplit(state_names{i},'/');
    model_info.coord_names.all(i)=state_names_split(4);
end
model_info.jointi.rotations=1:3;
model_info.jointi.translations=4:6;

%% Creating casadi functions
addpath([S.misc.main_path '\CasadiFunctions']);
disp('Start creating CasADi functions...')
disp(' ')
t0 = tic;
[f_casadi] = createCasadiFunctions(S);
disp(' ')
disp(['...CasADi functions created. Time elapsed ' num2str(toc(t0),'%.2f') ' s'])
disp(' ')
disp(' ')

%% Formulating OCP
addpath([S.misc.main_path '\OCP'])
if ~S.post_process.rerun
    OCP_formulation(S,model_info,f_casadi);
    disp(' ')
    disp(' ')
end

%% PostProcessing
addpath([S.misc.main_path '\PostProcessing'])
disp('Start PostProcessing...')
disp(' ')
t0 = tic;
PostProcessing(S,model_info);
disp(' ')
disp(['...PostProcessing done. Time elapsed ' num2str(toc(t0),'%.2f') ' s'])
disp(' ')
disp(' ')

%% Conclude diary
disp(['Total time elapsed ' num2str(toc(t00),'%.2f') ' s'])
disp(' ')
disp(['Diary saved as ' Outname])
disp(' ')
disp(' ')
diary off