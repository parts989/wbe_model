# wbe_model

Integrated mechanistic model describing dynamics of fecally shed signals into wastewater catchments

## Main Model Components and Description
- **Fecal Shedding Models** - Use fecal shedding data of SARS-CoV-2, PMMoV, and crAssphage as inputs. Generate interpolated shedding trajectories for SARS-CoV-2, static distributions for PMMoV and crAssphage
- **SLIR Model** - Generate model epidemiological data for an outbreak of SARS-CoV-2 to use as inputs for wastewater model
- **Wastewater Model** - Take the number of infected individuals and fecal shedding profiles to produce data for bulk wastewater measurements.


**Start with Basic model.R. Will produce a single iteration at a defined population** 
