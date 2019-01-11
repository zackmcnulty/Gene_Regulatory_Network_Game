### How to read an Antimony model
* `//` marks comments.
* Constants and species are declared at the beginning.
* Reactions are formatted as the following, with one reaction each line
  * Reaction_Name: [Chemical => Reaction ] ; [Rate equation] ;
* Species and parameters are initialized at the very end.  
# Assumptions with the Antimony model
* This model is **incomplete** and **missing information**.
* There is an external input on one gene, called `INPUT`, that affects the system.
* Any rate equation parts missing info has been replaced with `0.1` to allow the model to run in Tellurium. 
* If any gene has two inputs affecting it, and the information about one of them is missing, the reaction with the known input is written as a modified mass action reaction.

### ASSUMPTIONS
* You will be given rough, biologically reasonable ranges for the values of each parameter
* There are no outside players; all genes involved are already present
* There can be at most two regulators for each gene
* A single protein cannot act on both regulatory sites of another gene
* The number of regulatory sites a protein can act on is not limited
* INPUT does not act on any genes besides the one it is currently attached to
* A protein might act on nothing
* self-feedback loops are allowed
* The true network is fully connected
* Nothing can act on INPUT
* Your experimental data will have some error/noise
* All of the currently provided connections are correct
* Concerning perturbations: Perturbations have a minimum of 0% perturbation (exclusive), for obvious reasons. There is no upper limit, meaning that perturbations can be over 100%. In studies, scientists have upregulated genes by over thousand-fold.
* The following are parameter ranges that you can expect the parameters of the model to be in.

| Parameter | Value      |
|-----------|------------|
| d_p       | 0.01-0.03  |
| d_m       | 0.5-2.0    |
| L         | 0.01-0.03  |
| Vm        | 0.5-2.0    |
| a_p       | 0.05-0.15  |
| H         | 2-8        |
| K         | 0.01-0.03  |

### ADVICE
* Biotapestry (http://www.biotapestry.org/) has a convenient CSV format for building models (see the tutorials on their website
	for more information). If you write code to convert this CSV format to an antimony string,
	you will save yourself a lot of time trying to simulate potential connections. We also wrote
	some helper methods that can quickly add/remove connections from the CSV format.
* Generating test models (which you can then "break") are a great way to test your code. From these
 	 models, you can generate your own data, add noise, or do whatever you want to see how
	 robust your testing procedure is.