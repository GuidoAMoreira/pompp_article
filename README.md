This is the code to reproduce the results in the paper Presence-Only through Marked Point Process.

The code can take very long to run. The pompp R package version used is in the pompp_0.0.0.9000 tarball (source package with compilation).

## Artificial data

The steps to reproduce the results in the "Artificial data" section are:

1. Set simulated folder as the working directory.
2. Backup data and chains folders.
3. Empty data and chains folders contents.
4. Run simulator.R file. May be run multiple times in parallel sessions. It handles the data generation and storing by itself (tracking done in the files stored in the data folder. If the folder is not cleared, it will consider that all simulations are finished).
5. Repeat steps 1-4 for the contents of misspec folder inside the simulated folder.
6. File results.R contains the code to generate figures in the paper. It requires the contents of the data and chains folders (including the ones in the misspec folder) for exact replicas.

## Steps to reproduce the application

The data that support the findings of this study are available from the authors but restrictions apply to the availability of these data, which were used under license from Portuguese fisheries authority “DGRM - Direcção Geral dos Recursos Marinhos, Segurança e Serviços Marinhos” for the current study, and so are not publicly available. Data are available from the paper's third author upon reasonable request and with permission from the DGRM.

The steps to reproduce the results in the "Application" section are:

1. Set fishery folder as the working directory.
2. Run all code in the runModel.R file.
3. File readMCMC.R contains the code to generate figures in the paper.
