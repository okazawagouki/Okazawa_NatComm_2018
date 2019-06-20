# Code from Okazawa, Sha, Purcell, & Kiani (2018) Nat Com

This repository contains simulation code for the following paper:

Okazawa G, Sha L, Purcell BA, Kiani R (2018). [Psychophysical reverse correlation reflects both sensory and decision-making processes](https://www.nature.com/articles/s41467-018-05797-y). *Nature Communications* 9:3479

## Requirement

The code is tested under MATLAB R2018a (MathWorks).

## Main code

[DDM_Kernel_Simulation.m](./DDM_Kernel_Simulation.m): This code generates a psychophysical kernel based on drift diffusion model.

[RACE_Kernel_Simulation.m](./RACE_Kernel_Simulation.m): This code generates a psychophysical kernel based on competing accumulator (race) model.

## Script

[Simulation1_fixed_and_RT_task.m](./Simulation1_fixed_and_RT_task.m): The script shows psychophysical kernels based on drift diffusion model with/without the non-decision time under fixed duration/reaction time tasks (corresponds to Fig. 3 in the paper).

[Simulation2_urgency_and_race_model.m](./Simulation2_urgency_and_race_model.m): The script shows psychophysical kernels based on drift diffusion model with urgency and competing accumulator model with input correlation, reflective bound, and leak/mutual inhibition (corresponds to Fig. 6 in the paper).

## License

These files are being released openly in the hope that they might be useful but with no promise of support. If using them leads to a publication, please cite the paper.

The dataset is released under a CC-BY 4.0 license.

The code is released under a [BSD license](./LICENSE.md).
