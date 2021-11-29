The entire workflow is:

1. Run R --no-save --args < data.R to generate data file into the data folder;
2. Go to pml folder and run autosub.joint.sh script to run the first step of jointGES;
3. Go to indivi folder and run autosub.subset.sh script to run the second step of jointGES, also, run autosub.sep.sh to generate results of separate estimation;
4. Run R --no-save --args < indivi_eval.R to get SHD figures;
5. Go to roc folder and run autosub.subset.sh and autosub.sep.sh to get results for ROC figure;
6. Run R --no-save --args < tpfig.R to get ROC figures.
