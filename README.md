# multiple-contrast-count
R code for reproducing the simulations and figures as well as application code of the manuscript: Pigorsch, M., Hothorn, L. A., & Konietschke, F. (2025). Multiple Contrast Tests for Count Data: Small Sample Approximations and Their Limitations. Biometrical Journal, 67(6), e70098. https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.70098

* **What the project does:**
This is the R code used for the simulations in the manuscript "Multiple Contrast Tests for Count Data: Small Sample Approximations and Their Limitations.".
It additionally includes **application code** with an **update incorporating offsets** as well as the option to compute one-sided tests.
Moreover, additional simulation code examines simulation settings adjusted to the distribution of a count outcome referring to Naik, M. G., Budde, K., Koehler, K., Vettorazzi, E., Pigorsch, M., Arkossy, O., Stuard, S., Duettmann, W., Koehler, F. & Winkler, S. (2022). Remote patient management may reduce all-cause mortality in patients with heart-failure and renal impairment. Frontiers in Medicine, 9, 917466. https://www.frontiersin.org/journals/medicine/articles/10.3389/fmed.2022.917466/full

* **Why the project is useful:**
The simulation compares different methods for multiple contrast tests for count data. The application code allows users to use these methods for own data. 

* **How users can get started with the project:**
For application use the code in the Application folder, for reproducing the simulation results use the R files in the Simulation folder. 


The code was written by the first and last author. For any questions or comments,
please reach out to mareen.pigorsch@charite.de.

This folder contains the following data and files that can be used to reproduce the analysis and figures of the manuscript.
It contains three subfolders (Application, Results, Simulation) containing the following files:

**.\Application** \
A folder including code for application to reproduce application example and updated application code incorporating offsets as well as the option to compute one-sided tests. 

MCT_Count_Application.R  \
An R script including the function to execute the analyses for one data set and the code for the application used in the paper, 
creates Table 3 and 4 and Figure 4. 

MCT_Count_Application_function_boot.R \
An R script including the function to execute solely the Nonparametric Bootstrap for an application data set, including option to use offsets and deciding between one-sided and two-sided tests.

MCT_Count_Application_function_with_offset.R \
An R script including the function to execute all methods examined for an application data set, update including option to use offsets.

**.\Results** \
\Results_Simulation_Excel \
MCT_Count_typeI_3groups.xlsx \
An Excel file including the simulation results for the type-I error for 3 groups. 
		
MCT_Count_typeI_4groups.xlsx \
An Excel file including the simulation results for the type-I error for 4 groups. 

MCT_Count_power_anyall_4groups.xlsx \
An Excel file including the simulation results for any and all pairs power for 4 groups. 

MCT_Count_power_3groups.xlsx \
An Excel file including the simulation results for the global power for 3 groups. 

MCT_Count_power_4groups.xlsx \
An Excel file including the simulation results for the global power for 4 groups. 

MCT_Count_typeI_TIM-HF2.xlsx \
An Excel file including the simulation results for the type-I error for 4 groups for simulation settings referring count outcome of TIM-HF2 trial. 

\Results_Plots \
A folder including .png-files with the plotted simulation results. 
	
\Results_Tables \
A folder including Tables 1, 2, A1 and S1-S5 as Excel files. 

MCT_Count_plot_3.R \
An R script to create the figures for the simulation study regarding 3 groups, shown in the supplementary material. 

MCT_Count_plot_4.R \
An R script to create the figures for the simulation study regarding 4 groups, creates Figure 1-3 and A1-A5. 

MCT_Count_Powertable_anyall_4.R \
An R script to create the tabeles for the simulation study for any and all pairs power, creates Table S1-S4 in supplementary material. 

MCT_Count_Tables.R \
An R script to create Table 1, 2, A1 and S5. 

**.\Simulation** \
MCT_Count_functions.R \
An R script that contains all the functions used for the simulation study to compare different analysis strategies for count data in context of multiple contrast tests (for 3 and 4 groups).	

MCT_Count_sim_3groups_Power.R \
An R script that contains the simulation setup and calls the simulation function for power regarding 3 groups. 	The simulation is conducted on the high-performance compute (HPC) cluster from the Berlin Institute of Health after including the functions from MCT_Count_functions.R. 

MCT_Count_sim_3groups_TypeI.R \
An R script that contains the simulation setup and calls the simulation function for type-I error regarding 3 groups. 	The simulation is conducted on HPC cluster from the Berlin Institute of Health after including the functions from 	MCT_Count_functions.R. 		

MCT_Count_sim_4groups_Power.R \
An R script that contains the simulation setup and calls the simulation function for power regarding 4 groups. The simulation is conducted on the HPC cluster from the Berlin Institute of Health after including the functions from MCT_Count_functions.R.		

MCT_Count_sim_4groups_TypeI.R \
An R script that contains the simulation setup and calls the simulation function for type-I error regarding 4 groups. 	The simulation is conducted on the HPC cluster from the Berlin Institute of Health after including the functions from 	MCT_Count_functions.R.		

MCT_Count_sim_TypeI_TIM-HF2.R \
An R script that contains the function and simulation setup for type-I error regarding 4 groups for the count outcome referring to TIM-HF2. The simulation is conducted on the HPC cluster from the Berlin Institute of Health.	

MCT_Count_sim_example_reproduction.R \
An R script that contains examples to check reproducibility without repeating the whole simulation study, one for a setting regarding type-I error and one for a setting regarding power.  

MCT_Count_sim_parallel_without_cluster.R \
An R script that contains an alternative code calling the simulations for 4 groups without using a cluster computer but a parallelization with function future_lapply() instead of slurm_apply(), examplary for type-I error for 4 groups.
