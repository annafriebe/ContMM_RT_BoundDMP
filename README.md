# ContMM_RT_BoundDMP
Code and data related to the paper Continuous-Emission Markov Models for Real-Time Applications: Bounding Deadline Miss Probabilities.

An extended version of the paper is also available here.

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
### Reproducing figs 11 and 12

```console
user@host:~$cd ~/ContMM_RT_BoundDMP/fig
user@host:~/ContMM_RT_BoundDMP/fig$Rscript seqHistAcorrFigs.R
````
This produces figures 11 and 12 from the timing data in ContMM_RT_BoundDMP/data/traces_reports_csv/control_time_full.csv

Fig. 11 is reproduced in ContMM_RT_BoundDMP/fig/trace_sequence.png (and .eps)
Fig. 12 a)  is reproduced in ContMM_RT_BoundDMP/fig/trace_density.png (and .eps)
Fig 12 b) is reproduced in ContMM_RT_BoundDMP/fig/ggplot_acorr.png (and .eps)

### Fitting the HMM to the timing data
```console
user@host:~$cd ~/ContMM_RT_BoundDMP/contMMfit
user@host:~/ContMM_RT_BoundDMP/contMMfit$Rscript fitHMM.R
````
This fits a HMM to the timing data in ContMM_RT_BoundDMP/data/traces_reports_csv/control_time_full.csv

The fitted model is described by the files in ContMM_RT_BoundDMP/data/models/cont. 
Eq. 32 is reproduced in ContMM_RT_BoundDMP/data/models/cont/transitionMatrix.txt
The mean and standard deviation columns of Table IV are found in nanoseconds in the first and second columns of ContMM_RT_BoundDMP/data/models/cont/normalParams.txt. Multiplying those values with 1e-6 to convert to milliseconds gives the resulting parameters.
The Stationary Prob. Column of Table IV is found in the ContMM_RT_BoundDMP/data/models/cont/stationaryDistr.txt file.

### Simulation with independence assumption

### Simulation with continuous model

### Evaluating the DMP bound
```console
user@host:~$cd ~/ContMM_RT_BoundDMP/eval/Bound
user@host:~/ContMM_RT_BoundDMP/eval/Bound$python3 DMPAccum8State.py
````
This runs the accumulation process for the DMP bound, using the fitted model in ContMM_RT_BoundDMP/data/models/cont/
The resulting files are saved to ContMM_RT_BoundDMP/data/dmp_bound

### PROSIT

### Generating fig 13
```console
user@host:~$cd ~/ContMM_RT_BoundDMP/fig
user@host:~/ContMM_RT_BoundDMP/fig$python3 genEvalFig.py
````
Fig. 13 is reproduced in eval_fig.png (and .eps) from files in ContMM_RT_BoundDMP/data.

## How to run the Furuta pendulum control tests on a Raspberry pi

### Installation

### Running the test program





