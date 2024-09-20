# Multi-objective optimal control and bilevel optimization for predictive simulations of walking

This project includes MATLAB scripts and functions designed to generate and analyze predictive simulations of walking of a planar five-link biped model. Each simulation is obtained by solving a multi-objective optimal control problem (MOCP), where the control policy is defined by a special weighted combination of selected cost functionals. Additionally, the project demonstrates an example of bilevel optimization used as an inverse optimal control method to identify the optimal control policy (specifically, the optimal weights for the selected cost functionals) that best replicates motion data from walking trials.

## Dependencies
This project relies on the following:
- Microsoft Windows 10 or 11.
- MATLAB (R2021a or later recommended) with *Optimization Toolbox™*, *Global Optimization Toolbox™*, *Parallel Computing Toolbox™*, and *CasADi*.
- The `walkingData.mat` file, which contains preprocessed gait data.

## Main scripts
The following three main scripts are provided, which can be run directly and independently, possibly in the order listed below. All of them include parallelization to speed up computations.
- **exploreParetoFront.m**:
    - Estimates the ideal and nadir objective vectors for the two objective functionals in `solveMOCP`, computes a user-specified number of points on the Pareto front, and visualizes the results.
- **bilevelOptimization.m**:
    - This is the main script for performing bilevel optimization. It identifies the optimal weights for the two objective functionals in `solveMOCP` to best match the walking kinematic data stored in the MAT-file `walkingData.mat`. Such optimal weights minimize the sum of the squared deviations (see `kinDeviation.m`) between the segment angles of the biped model, as recorded in `walkingData.mat`, and the corresponding angles predicted by `solveMOCP`.
- **bilevelOptimizationChallenge.m**:
    - Uses bilevel optimization to identify the optimal weights for the five objective functionals in `solveMOCPchallenge`, with the goal of minimizing the peak joint torque across the five torque motors at the joints. This kind of problem was originally presented to participants in the [ISB21 Workshop](https://github.com/antoinefalisse/ISB21-workshop) as a challenge.

## Functions
- **solveMOCP.m**:
    - Formulates and solves an MOCP for a predictive simulation of walking using a planar five-link biped model. Two objective functionals are implemented, namely the time integral of the squared 2-norm of the joint torques and that of the angular accelerations of the model segments. The two objectives are combined into a single objective using a weighted Chebyshev metric.
- **solveMOCPchallenge.m**:
    - Similar to `solveMOCP.m`, but includes five objective functionals, each of them being the time integral of the squared 2-norm of the *i*-th joint torque (provided by the ideal torque motor at the *i*-th joint of the model). The five objectives are combined into a single objective using the classical weighting method.
- **kinDeviation.m**:
    - Calculates the sum of the squared deviations between the segment angles of the biped model, read from `walkingData.mat`, and the corresponding angles obtained by solving a predictive MOCP, i.e., using `solveMOCP.m` with given weights.
	
The remaining `.m` functions included in this project are the original ones from the [ISB21 Workshop](https://github.com/antoinefalisse/ISB21-workshop) project.

## Acknowledgements
This project builds upon previous work by **Tom Van Wouwe**, **Antoine Falisse**, and **Gil Serrancoli**. Please see their GitHub project: [ISB21 Workshop](https://github.com/antoinefalisse/ISB21-workshop).

## Contact
If you have any questions or inquiries, please feel free to contact **Alessio Artoni** at alessio.artoni@unipi.it.