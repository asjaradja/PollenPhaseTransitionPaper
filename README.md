# PhaseDiagramCalculations
This is a repository of the code essential for the [TITLE] paper. It contains the following files and folders:
	* PhaseDiagramCalculations: All code for building the phase diagram and finding kinetically arrested states
		* GradientDescent_1l: Gradient descent for 1 el value used to find equilibrium states
	⁃	GradientDescent_2l: Gradient descent for 2 el values used to find equilibrium states
	⁃	SimulatedAnnealing_1l: Simulated annealing for 1 el value used to find equilibrium states
	⁃	SimulatedAnnealing_2l: Simulated annealing for 2 el values used to find equilibrium states
	⁃FiPyPhaseFieldCrystalSphere.ipynb: Finite volume solver used to find kinetically arrested states
	•	TraitReconstruction
	⁃	3State_WavelengthModel
	⁃	SpermatophyteTree.tre: Tree file
	⁃	SpermatophyteTree_wavelengthhypothesis.txt: Data file
	⁃	2State_EquilibriumModel
	⁃	SpermatophyteTree.tre: Tree file
	⁃	SpermatophyteTree_equilibriumhypothesis.txt: Data file
	⁃	NexusFile_Spermatophytes.tree: Nexus file for ancestral trait reconstruction

Each of the folders contains the following source code:
	•	Makefile
	•	main.cpp: main program
	•	Header files from Numerical Recipes [1] necessary for gradient descent and simulated annealing:
	⁃	nr.h
	⁃	nrtypes.h
	⁃	nrtypes_nr.h
	⁃	nrutil.h
	⁃	nrutil_nr.h
	•	surface_pattern_1l.py/surface_pattern_2l.py: plots the surface pattern onto a sphere
	•	wigxjpf.h (for Simulated Annealing): header file for the evaluation of wigner 3j coefficients

REQUIREMENTS
In order to successfully run main.cpp you need to have the following installed:
	•	Numpy, scipy, FiPy, and matplotlib (Python)
	•	wigxjpf library and program for evaluation of wigner 3j coefficients (http://fy.chalmers.se/subatom/wigxjpf/)
	•	Eigen: a c++ template library for linear algebra (http://eigen.tuxfamily.org/)
	•	BayesTraits software (http://www.evolution.rdg.ac.uk/BayesTraitsV3.0.1/BayesTraitsV3.0.1.html)

USAGE
1. To compile and generate an executable file, “make” is sufficient. This creates an executable file that 
	can be run simply with: ./main.
	The executable generates several text files:
	⁃	T_H.txt: 2 columns, annealing temperature and lowest Hamiltonian value at every simulated annealing step
	⁃	T_dH.txt: 2 columns, annealing temperature and Hamiltonian derivative value at every simulated annealing step
	⁃	endcms.txt: the complex coefficients (cms) for the final pattern
	⁃	cms.txt: the complex coefficient (cms) for the pattern at every annealing step
	⁃	parameters.txt: Contains information on the parameters of the simulation (el_not, lambda3, etc.), the starting cos, the eigenvalues, the lowest point derivative value and the lowest hamiltonian value
	⁃	time.txt: simulation time

2. surface_pattern_1(2)l.py reads in two input files, endcms.txt and parameters.txt, and runs:
	python surface_pattern_1(2)l.py endcms.txt parameters.txt
    You can show the plot with “plt.show()” or save the figure with “plt.savefig([file name])”

3. Run FiPyPhaseFieldCrystalSphere.ipynb in jupyter to find the kinetically arrested states

4. To run the ancestral trait reconstruction, run the following command with the included tree and data            
    files:
	./BayesTraitsV3 TreeFile DataFile


[1] Press, W.H., Teukolsky, S.A., Vetterling, W.T., & Flannery, B.P. Numerical Recipes: The Art of Scientific Computing. Cambridge [Cambridgeshire]: Cambridge University Press (1986)
