# ContMM_RT_BoundDMP
Code and data related to the paper Efficiently bounding deadline miss probabilities of Markov chain real-time tasks,
extended from Continuous-Emission Markov Models for Real-Time Applications: Bounding Deadline Miss Probabilities.

Link to paper: https://link-springer-com.ep.bib.mdh.se/article/10.1007/s11241-024-09431-7


## How to reproduce the main results of the paper from the data in this repository

### Installation of a virtual machine

The analysis has been performed with an Ubuntu virtual machine, with at least 1GB RAM.
A virtual machine from osboxes can be used.

Install git
```console
user@host:~$sudo apt install git
````
Install R, packages and dependencies
```console
user@host:~$sudo apt install install libcurl4-openssl-dev 
user@host:~$sudo apt install -no-install-recommends rbase
user@host:~$R
>install.packages(‘ggplot2’)
>install.packages(‘forecast’)
>install.packages(‘modules’)
>install.packages(‘depmixS4’)
>install.packages(‘data.tree’)
>q()
````

Install scipy, numpy
```console
user@host:~$sudo apt install python3-numpy
user@host:~$sudo apt install python3-scipy
````

Clone this repository
```console
user@host:~$git clone https://github.com/annafriebe/ContMM_RT_BoundDMP.git
````
Clone PROSITool
```console
user@host:~$git clone https://bitbucket.org/luigipalopoli/prositool.git
````
### Reproducing figs 16, 17 and 18

```console
user@host:~$cd ~/ContMM_RT_BoundDMP/fig
user@host:~/ContMM_RT_BoundDMP/fig$Rscript seqHistAcorrFigs.R
````
This produces figures 16, 17 and 18 from the timing data in ContMM_RT_BoundDMP/data/traces_reports_csv/control_time_full.csv

Fig. 16 is reproduced in ContMM_RT_BoundDMP/fig/trace_sequence.png (and .eps)
Fig. 17  is reproduced in ContMM_RT_BoundDMP/fig/trace_density.png (and .eps)
Fig 18 is reproduced in ContMM_RT_BoundDMP/fig/ggplot_acorr.png (and .eps)

### Fitting the HMM to the timing data
```console
user@host:~$cd ~/ContMM_RT_BoundDMP/contMMfit
user@host:~/ContMM_RT_BoundDMP/contMMfit$Rscript fitHMM.R
````
This fits a HMM to the timing data in ContMM_RT_BoundDMP/data/traces_reports_csv/control_time.csv

The fitted model is described by the files in ContMM_RT_BoundDMP/data/models/cont. 
Eq. 41 is reproduced in ContMM_RT_BoundDMP/data/models/cont/transitionMatrix.txt
The mean and standard deviation columns of Table 3 are found in nanoseconds in the first and second columns of ContMM_RT_BoundDMP/data/models/cont/normalParams.txt. Multiplying those values with 1e-6 to convert to milliseconds gives the resulting parameters.
The Stationary Prob. Column of Table 3 is found in the ContMM_RT_BoundDMP/data/models/cont/stationaryDistr.txt file.

### Simulation with independence assumption
```console
user@host:~$cd ~/ContMM_RT_BoundDMP/eval/Ind
user@host:~/ContMM_RT_BoundDMP/eval/Ind$Rscript DMPEstCBSInd.R
````
This estimates the DMP with an independence assumption, by randomly reordering the data in 
ContMM_RT_BoundDMP/data/traces_reports_csv/control_time.csv before entering it into a CBS simulation.
The result is output in ContMM_RT_BoundDMP/data/ind/result.csv.

### Simulation with continuous model
```console
user@host:~$cd ~/ContMM_RT_BoundDMP/eval/Sim-Cont
user@host:~/ContMM_RT_BoundDMP/eval/Sim-Cont$Rscript DMPSimContMM.R
````
This estimates the DMP from the fitted HMM by generating execution times that are entered into a CBS simulation.
The resulting files are saved to ContMM_RT_BoundDMP/data/sim_cont.

### Simulation with merged continuous model
```console
user@host:~$cd ~/ContMM_RT_BoundDMP/eval/Sim-Cont
user@host:~/ContMM_RT_BoundDMP/eval/Sim-Cont$Rscript DMPSimContMM_merged.R
````
This simulation is used to obtain the starting beta values for the merged model below.


### Evaluating the DMP bound with the 8-state model
```console
user@host:~$cd ~/ContMM_RT_BoundDMP/eval/Bound
user@host:~/ContMM_RT_BoundDMP/eval/Bound$python3 DMPAccum8State.py
````
This runs the accumulation process for the DMP bound, using the fitted model in ContMM_RT_BoundDMP/data/models/cont/
The resulting files are saved to ContMM_RT_BoundDMP/data/dmp_bound/pzwl_dmp_control_bound_XXXX.csv.
Computation times for bounds are saved to ContMM_RT_BoundDMP/data/pzwl_dmp_control_5time.csv and ContMM_RT_BoundDMP/data/pzwl_dmp_control_10time.csv.

### Evaluating the DMP bound with the merged 2-state model
```console
user@host:~$cd ~/ContMM_RT_BoundDMP/eval/Bound
user@host:~/ContMM_RT_BoundDMP/eval/Bound$python3 DMPAccum8StateMerged.py
````
This runs the accumulation process for the DMP bound, using a merged 2-state model in ContMM_RT_BoundDMP/data/models/cont/
The resulting files are saved to ContMM_RT_BoundDMP/data/dmp_bound/pzwl_dmp_control_merged_bound_XXXX.csv.
Computation times for bounds are saved to ContMM_RT_BoundDMP/data/pzwl_dmp_control_merged_5time.csv and ContMM_RT_BoundDMP/data/pzwl_dmp_control_merged_10time.csv.

### PROSIT
```console
user@host:~$cp ContMM_RT_BoundDMP/eval/PROSIT/Makefile prositool
user@host:~$cd prositool
user@host:~/prositool$make
user@host~/prositool$cd ~
user@host:~$cp prositool/cli_solver ContMM_RT_BoundDMP/eval/PROSIT/
user@host:~$cd ContMM_RT_BoundDMP/eval/PROSIT/
user@host:~/cd ContMM_RT_BoundDMP/eval/PROSIT$./run_cli_solver.sh
````
Runs the CLI solver from PROSITool with the model in ContMM_RT_BoundDMP/data/models/discr. The results are saved to ContMM_RT_BoundDMP/data/prosit/.

### Generating fig 19
```console
user@host:~$cd ~/ContMM_RT_BoundDMP/fig
user@host:~/ContMM_RT_BoundDMP/fig$python3 genEvalFig.py
````
Fig. 19 is reproduced in eval_fig.png (and .eps) from files in ContMM_RT_BoundDMP/data.

Table 4 contains times stored when calculating the bounds. The first table row is taken from the last rows of ContMM_RT_BoundDMP/data/pzwl_dmp_control_merged_5time.csv.
The second table row is taken from the last rows of pzwl_dmp_control_merged_10time.csv, the third from pzwl_dmp_control_5time.csv and the last row from pzwl_dmp_control_10time.csv.

## How to run the Furuta pendulum control tests on a Raspberry pi

### Installation

The evaluation has been performed on a Raspberry Pi Model 3B+, with a kernel patched with PREEMPT_RT fetched from https://github.com/kdoren/linux/releases/
The used kernel is the 32-bit kernel 5.15.65-rt49-v7l+

Clone the repository. Build the programs in test_program/controller and test_program/simulator

### Running the test program
Traces are obtained by running test_program/runSimulatorControllerTrace.sh

Tests are run with the various test_program/runSimulatorControllerBW_XX.sh files.




