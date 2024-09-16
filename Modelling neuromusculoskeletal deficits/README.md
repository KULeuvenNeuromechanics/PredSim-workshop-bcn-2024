# Modelling neuromusculoskeletal deficits

Physics-based computer simulations that can predict the effect of treatments (e.g., bony and soft tissue correction, ankle-foot-orthoses) on gait in children with cerebral palsy (CP) have the potential to improve clinical decision-making. To this end, an important challenge is to accurately personalise patient-specific neuromusculoskeletal models.

## Example 1. Modelling msk impairments from clinical exam

The clinical exam is part of children's usual clinical care and is a comprehensive assesment of musculoskeletal functioning. In example 1. you will use (I.) passive Range of Motion (ROM) and (II.) to personalize optimal muscle force and optimal muscle fiber length, respectively.

Passive range of motion and muscle strength scores are provided in the 'Clinical Exam'/BCN_CP#:

I. Passive range of motion (pROM)
- Soleus (Silfversköld Test):
	The patient lies supine, with one knee flexed at 90°. The ankle is moved into maximal dorsiflexion. The maximal dorsiflexion angle is measured (evaluation of the length of soleus). The typical value for dorsiflexion in this 	situation is 20 to 30 °. (pROM_Ankledf 90_#side)
- Gastrocnemii (Silfversköld Test): 
	Tested in a similar position as Soleus. The knee is moved in full extension. At that time we can evaluate whether the length of the gastrocnemius limits the ankle motion, thereby decreasing the maximal dorsiflexion angle. The typical value 	for dorsiflexion in this situation is 	10 to 20°. (pROM_Ankledf 0_#side)
- Hamstrings (popliteal angle - bilateral):
	The patient lies supine. The evaluated limb is flexed at the level of the hip. The contralateral limb is flexed. Starting from knee flexion, the knee of the evaluated limb is moved into maximum extension. The deficit until full knee extension is noted (as a negative angle: lack to full extension). The typical value for the bilateral test is -15° to -0°. (p_ROM_Poplbi_#side)
	

II. Strength
	Strength is evaluated by manual muscle testing. The strength is evaluated for the full active range of motion. 
 	
	Strength scores:
	1. Evidence of slight contraction of the muscle but joint motion is not visible
	2. Complete range of motion in gravity eliminated plane
	3. Perfect motion against gravity 
	4. Motion against gravity with some (moderate resistance)
	5. Motion against gravity with maximal resistance

In this example, the code main_mskClinicalExam.m (PredSim-workshop-bcn-2024/Modelling neuromusculoskeletal deficits/Code/Example 1 - msk Clinical exam) .
will guide you through the estimation process. You only have to edit the lines of code that are inbetween % ------ start edit ----- and % ----- end edit -----
 

## Example 2. Modelling msk impairments from data-driven EMG torque relationships

(Description)



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


# References:

Fowler, E. G., Staudt, L. A., Greenberg, M. B., & Oppenheim, W. L. (2009). Selective Control Assessment of the Lower Extremity (SCALE): development, validation, and interrater reliability of a clinical tool for patients with cerebral palsy. Developmental Medicine & Child Neurology, 51(8), 607-614.

M.H. Schwartz, A. Rozumalski, K.M. Steele, ”Dynamic motor control is associated with treatment outcomes for children with cerebral palsy,” Dev. Medicine & Child Neurology, 58(11), 2016, 1139-1145