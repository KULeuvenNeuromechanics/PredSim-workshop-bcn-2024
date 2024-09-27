# Modelling neuromusculoskeletal deficits

Physics-based computer simulations that can predict the effect of treatments (e.g., bony and soft tissue correction, ankle-foot-orthoses) on gait in children with cerebral palsy (CP) have the potential to improve clinical decision-making. To this end, an important challenge is to accurately estimate patient-specific neuromusculoskeletal models.

## Example 1. Modelling msk impairments from clinical exam

The clinical exam is part of children's usual clinical care and is a comprehensive assesment of musculoskeletal functioning. In example 1. you will use (I.) manual muscle testing strength scores and (II.) passive Range of Motion (ROM) scores to personalize optimal muscle force and optimal muscle fiber length, respectively.

**Data:** Passive range of motion and muscle strength scores, provided in the subfolder 'Clinical Exam'/BCN_CP#:

### I. Strength
The strength is evaluated for the full active range of motion by manual muscle testing. 
 	
Strength scores indicate:
1. Evidence of slight contraction of the muscle but joint motion is not visible
2. Complete range of motion in gravity eliminated plane
3. Perfect motion against gravity 
4. Motion against gravity with some (moderate resistance)
5. Motion against gravity with maximal resistance

### II. Passive range of motion (pROM)
- Soleus (Silfversköld Test):
	The patient lies supine, with one knee flexed at 90°. The ankle is moved into maximal dorsiflexion. The maximal dorsiflexion angle is measured (evaluation of the length of soleus). The typical value for dorsiflexion in this situation is 20 to 30 °. (pROM_Ankledf 90_#side)
- Gastrocnemii (Silfversköld Test): 
	Tested in a similar position as Soleus. The knee is moved in full extension. At that time we can evaluate whether the length of the gastrocnemius limits the ankle motion, thereby decreasing the maximal dorsiflexion angle. The typical value for dorsiflexion in this situation is 10 to 20°. (pROM_Ankledf 0_#side)
- Hamstrings (popliteal angle - bilateral):
	The patient lies supine. The evaluated limb is flexed at the level of the hip. The contralateral limb is flexed. Starting from knee flexion, the knee of the evaluated limb is moved into maximum extension. The deficit until full knee extension is noted (as a negative angle: lack to full extension). The typical value for the bilateral test is -15° to -0°. (p_ROM_Poplbi_#side)

**Requirements:** OpenSim, CasADi.

**How to use the code:**
The code main_mskClinicalExam.m (PredSim-workshop-bcn-2024/Modelling neuromusculoskeletal deficits/Code/Example 1 - msk Clinical exam) guides the users through the estimation process. Users only have to edit the lines of code that are inbetween % ------ start edit ----- and % ----- end edit -----

**Running PredSim with estimated muscle-tendon parameters:**
Users have two options to run PredSim with the updated parameters.
1. Copy the structures `S.subject.muscle_strength line 58-line 70` and `S.subject.scale_MT_params - line 297-302` and paste in the %% Settings cell of PredSim main.m
2. Load the saved structures in the %% Settings cell of PredSim main.m.
 `S.subject.muscle_strength = load(fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name,'_muscle_strength.mat']));`
 `S.subject.scale_MT_params = load(fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name,'_MT_params.mat']));`

⚠️ NOTE: you also have to specify other subject specific parameters - see [Running PredSim](#running-predsim)
	
## Example 2. Modelling msk impairments from data-driven EMG torque relationships

Mechanical properties of muscles in children with CP are often altered. This code allows users to estimate the optimal fiber length, tendon slack length and tendon stiffness of each muscle such that the EMG-torque relationship across each joint can be satisfied. The code used for estimating these properties is available on https://github.com/KULeuvenNeuromechanics/MuscleRedundancySolver. This example will use data of walking gait and clinical tests to estimate the muscle-tendon properties for two children with CP. Data of two clinical tests are used. Instrumented passive Spasticity Assessment (IPSA) [Bar-On et al., 2013] and pendulum test [Fowler et al., 2000].

**Data:** Experiments record the EMG of major muscle groups is measured along with movement data and external forces for each trial. Movement data and external forces can then be used to run inverse kinematics and inverse dynamics to get the joint torques. These data along with the scaled OpenSim models for both patients is provided. 

**Requirements:** OpenSim, CasADi, clone of https://github.com/KULeuvenNeuromechanics/MuscleRedundancySolver repository.

**How to use the code:**

**Required inputs:**
In the code ParameterEstimation_BCN_workshop.m, specify the following:
1.	Subject name: either ‘CP1’ or ‘CP2’. 
2.	Path of your local directory of `PredSim-workshop-bcn-2024\S4 Modelling neuromusculoskeletal deficits`
3.	Path of the directory of MuscleRedundancySolver.
4.	Path of your CasADi folder
6.	A name to your analysis in Misc.AnalysisID

**Optional inputs:**
All the optional inputs are described on the https://github.com/KULeuvenNeuromechanics/MuscleRedundancySolver page and can be changed in the ParameterEstimation.m file.

**Outputs:**
The code will save the estimated muscle parameters, the data and settings in a .mat file. The code will also write the estimated optimal fiber length and estimated tendon slack length to the OpenSim model called 
`BCN_CP<1 or 2>_<Misc.AnalysisID>_paramEst.osim`. The code will also write the scaling factors of the optimal fiber length, tendon slack length and tendon stiffness in .mat files called `BCN_CP<1 or 2>_<Misc.AnalysisID>_paramEst_optimal_fiber_length_scale.mat`, `BCN_CP<1 or 2>_<Misc.AnalysisID>_paramEst_tendon_slack_length_scale.mat`, and `BCN_CP<1 or 2>_<Misc.AnalysisID>_paramEst_tendon_stiffness_scale.mat`.

**Running PredSim with estimated muscle-tendon parameters:**
Users have two options to run PredSim with the updated parameters.
1.	Use the updated OpenSim model (OpenSim model called `BCN_CP<1 or 2>_<Misc.AnalysisID>_paramEst.osim`) that has the estimated optimal fiber length and estimated tendon slack length already written in it, along with `BCN_CP<1 or 2>_<Misc.AnalysisID>_paramEst_tendon_stiffness_scale.mat` to set the S.subject.tendon_stiff_scale setting in PredSim.
2.	Use the .mat files corresponding to estimated scaling factors of optimal fiber length and tendon slack length (`BCN_CP<1 or 2>_<Misc.AnalysisID>_paramEst_optimal_fiber_length_scale.mat` and `BCN_CP<1 or 2>_<Misc.AnalysisID>_paramEst_tendon_slack_length_scale.mat`) to set the subject.scale_MT_param setting of PredSim, along with `BCN_CP<1 or 2>_<Misc.AnalysisID>_paramEst_tendon_stiffness_scale.mat` to set the S.subject.tendon_stiff_scale setting in PredSim.

To ensure that you ran the program correctly, we have provided some reference results for CP1. Your results after running PredSim using the model with personalized muscle tendon parameters using the default settings provided should match the reference result (https://github.com/KULeuvenNeuromechanics/PredSim-workshop-bcn-2024/blob/main/S4%20Modelling%20neuromusculoskeletal%20deficits/Code/Example%202%20-%20ParameterEstimation/ReferenceResult/PredSimResultForCP1.jpg). You can use the plotPredSimResults.m file provided in the Code folder to generate this result after your PredSim has converged. The section Results Analyses describes how to use the plotPredSimResults.m file.

## Example 3. Modelling neural impairments through muscle synergies

Muscle co-activation patterns derived from synergies might capture non-selective muscle control in children with CP and offer a way to include motor control deficits in predictive simulation workflows.

Scores from two clinical tests that evaluate motor control are provided in 'Clinical Exam' folder:
1. SCALE test
    Selective control assessment of the lower extremity was performed by the SCALE [Fowler et al., 2009].
    Normal selective voluntary motor control (SVMC) can be defined as the ability to perform isolated joint movement without using mass flexor⁄ extensor patterns or undesired movement at other joints, such as mirroring. The Selective Control Assessment of the Lower Extremity (SCALE) is a clinical tool developed to quantify SVMC in patients with CP. 
2. Selective test (scores explanation)
    0: no selective control, no (or minimal) contraction of some of the demanded muscles
    0.5: small contraction, but almost no motion, and/or a lot of co-contraction
    1: mild selective control, not all muscles working in a correct way, no smooth motion, with co-contraction (not always), limited range
    1.5 good co-contraction with correct muscles
    2: perfect control, perfect contraction with correct muscles

Results from the synergy analyses are provided in 'Models/BCN_CP#' folder:
Muscle synergies were extracted from the EMG signals of eight muscles using non-negative matrix factorisation. We selected the number of synergies that were needed to explain at least 90% of the variance accounted for (VAF) of the measured EMG. 
The eight measured muscles, for which the synergy analysis has been done, are: 'rect_fem', 'vasti_r', 'bifemsh_r', 'hamstrings_r', 'tib_ant_r', 'gastroc_r', 'soleus_r', 'glut_max_r'. The provided results in 'BCN_CP#_Syn.mat' are:
1. Number of synergies per leg (SynN.R and SynN.L)
2. Synergy weights per each muscle and synergy (SynW.R and SynW.L)
3. The VAF of only one synergy (VAF_1_syn.R and VAF_1_syn.L) during walking was computed as a measure of dynamic motor control [Schwartz et al., 2016].

In the example, two types of predictive simulations (additionally to the baseline simulation, without imposing synergies constraints) can be run:
1. Only imposing the number of synergies: The number of synergies per leg are imposed, by constraining all muscle activations to be controlled by a fixed number of synergies. This is done by adding a term in the cost function that minimises the difference between muscle activations in the optimisation and muscle actiovations reconstructed from the synergy activations and synergy weights. Additionally, this difference is constrained in an inequality constraint.
2. Imposing the number of synergies and tracking synergy weights : The patient-specific synergy weights (co-activation patterns) obtained from the synergy analysis are tracked fr the eight measured muscles. This is done by adding a term in the cost function that minimises the error between the synergy weights in the optimisation and the synergy weights from the synergy analysis.

## Running PredSim
Along with adjusting the parameters mentioned above, users will also need to specify the following settings for PredSim:
1. Initialize settings of the gait1018 model
2. Make sure to use the correct osim model (BCN_CP1_PredSim or BCN_CP2_PredSim) in main.m
3. Initial guess: Initial guess from inverse kinematics of 100% gait cycle can be used. The initial guesses are stored in `<path to PredSim-workshop-bcn-2024>\S4 Modelling neuromusculoskeletal deficits\Data\Data_<CP1 OR CP2>\<IK\<CP4_T0_10_IK_adjusted.mot for subject CP1 and CP16_T0_11_IK_adjusted.mot for subject CP2>`. 
4. While using initial guess base on inverse kinematics, it is advisable to adjust pelvis height of the initial guess `S.subject.adapt_IG_pelvis_y = 1;`.
5. Forward velocity (`S.misc.forward_velocity`): For CP1, the forward velocity is 1.1240 m/s. For CP2, the forward velocity is 1.0314 m/s.
6. Simulate a full gait cycle (`S.misc.gaitmotion_type = 'FullGaitCycle';`) instead of the default half gait cycle since the model is not symmetric.
7. Specify your casadi path

## Results Analyses:
After generating PredSim simulations, users can use the plotPredSimResults.m file in the Code folder to compare their results to IK results and result of PredSim run without any parameters estimated (Generic - no personalization). Users are allowed to compared multiple PredSim outputs at the same time. User can run multiple parameter estimations with varying settings and then the corresponding PredSim. This code can then be used to analyze how the predicted kinematics change with the different settings. The plotPredSimResults.m has the following required settings;
1.	subject: Subject name ‘CP1’ or ‘CP2’
2.	paramEstModelName: Names of the OpenSim models used to run the PredSim of each of the comparison.
3.	paramEstSuffix: PredSim outputs are followed by v<number> or job<number>. Please add this information here for each result
4.	modelLegend: name that the user want to add to identify each PredSim result
5.	path of your local directory of `PredSim-workshop-bcn-2024\S4 Modelling neuromusculoskeletal deficits`
6.	predsimResultsPath: path of PredSimResults folder where the PredSim results get stored


# References:

Fowler, E. G., Staudt, L. A., Greenberg, M. B., & Oppenheim, W. L. (2009). Selective Control Assessment of the Lower Extremity (SCALE): development, validation, and interrater reliability of a clinical tool for patients with cerebral palsy. Developmental Medicine & Child Neurology, 51(8), 607-614.

M.H. Schwartz, A. Rozumalski, K.M. Steele, ”Dynamic motor control is associated with treatment outcomes for children with cerebral palsy,” Dev. Medicine & Child Neurology, 58(11), 2016, 1139-1145

Bar-On, R., & Fiedeldey-Van Dijk, C. (2022). The Bar-On model and multifactor measure of human performance: Validation and application. Frontiers in psychology, 13, 872360.

Fowler, E. G., Nwigwe, A. I., & Ho, T. W. (2000). Sensitivity of the pendulum test for assessing spasticity in persons with cerebral palsy. Developmental medicine and child neurology, 42(3), 182-189.
