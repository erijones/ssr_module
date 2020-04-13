Supplemental Code for the Steady-State Reduction (SSR) and SSR-generated Parameter Change (SPARC) Methods
============================

This repository stores supplemental code for the papers:

+ "Steady State Reduction of generalized Lotka-Volterra systems in the
  microbiome," by Eric Jones and Jean Carlson, published in Phys. Rev. E in
  2019, (called the SSR paper) and
+ "Control of ecological outcomes through deliberate parameter changes in a
  model of the gut microbiome", by Zipeng Wang, Eric Jones, Joshua Mueller, and
  Jean Carlson, published in Physical Review E in 2020 (called the SPARC
  paper).

The code in `steady_state_reduction_example.py` generates Fig. 2 of the SSR
paper.  Steady State Reduction (SSR) approximates a high-dimensional gLV system
by two-dimensional (2D) subspaces that each obey gLV dynamics and span a pair
of high-dimensional steady states. Therefore, the process of SSR determines the
2D parameters that best approximate in-plane dynamics of the original
high-dimensional system.

The code in `SPARC_example.py` generates Fig. 3 of the SPARC paper. The control
framework SPARC (SSR-generated Parameter Change) alters the steady state
outcome in bistable regions of high-dimensional gLV systems. SPARC operates by
deliberately changing interaction parameters of the high-dimensional model:
once the parameter change is made, the dynamics of an initial condition that
was initially flowing towards some steady state A will be altered so that it
now flows towards some steady state B.  Here it is applied to an 11-dimensional
generalized Lotka-Volterra (gLV) model fit to gut microbial data from a C.
difficile infection (CDI) mouse experiment.

`steady_state_reduction_example.py`
-----------------------------------

The structure of this code is, briefly:

1) Import the parameters that characterize a high-dimensional system (growth
   rates \rho and interactions K), and store them in the class `Params`. In this
   code, as in the paper, we use the 11-dimensional parameters fit by Stein et al.
   (PLOS Comp Bio, 2013). To apply SSR to your own high-dimensional system, create
   your own `Params` class specific to your data.

2) Generate a pair of steady states to use in SSR. This code finds steady
   states of the high-dimensional model with `get_all_stein_steady_states`, but to
   apply SSR to a system you are interested in, identify the pair of
   high-dimensional steady states of your particular system you would like to
   study.
3) Generate the 2D parameters with SSR. Parameters of this 2D system are stored
   in the class `ssrParams`; generate this class using the parameters and steady
   state of your interest.
4) Plot and compare (in-plane) trajectories in the high-dimensional and
   SSR-reduced systems, saved in `SSR_demo_1.pdf`.
5) Plot the separatrix of the SSR-reduced 2D system, saved in `SSR_demo_2.pdf`.
6) Plot the in-plane separatrix of the high-dimensional system, saved in
   `SSR_demo_2.pdf`. Depending on the spatial resolution of the generated
   high-dimensional in-plane separatrix, this simulation can take a while-- use
   the `load_data` parameter to save a generated separatrix to the disk.

To use SSR yourself, generate your own `Params` class, say `p`, and generate the pair of
steady states you intend to use for SSR, say `ssa` and `ssb`. Then generate new
ssrParams with `s = ssrParams(p, ssa, ssb)`. You can check the validity of SSR
as implemented in your own system with the three `plot_*` commands.

`SPARC_example.py`
-----------------------------------

In total, this code:

1) imports the SSR formulae by Jones and Carlson (Phys. Rev. E 2019)
   (referred to as the SSR paper) and gLV system parameters fit by Stein et al.
   (PLOS Comp. Bio. 2013);
2) compresses the high-dimensional parameters into 2D SSR-generated
   parameters and generate trajectories in both systems starting at the initial
   condition (0.5,0.5);
3) finds a 2D parameter change that alters the steady-state outcome of this
   initial condition and plots the resulting new trajectory in the reduced 2D
   system; and
4) generates a high-dimensional parameter modification that corresponds to
   the 2D parameter change and plots the corresponding trajectory.

To apply the SPARC method to your own gLV model, in the main function at the
end of the file, create your own `Params` container class with your own growth
rates `rho`, interaction matrix `K`, and labels `labels` with the command `my_p
= ssr.Params([labels, rho, K])`. Name the two steady states that define your
bistable region of interest `ssa` and `ssb`, and compute the SSR-reduced
parameter set with the command `my_s = ssr.ssrParams(my_p, ssa, ssb)`. Finally, to
generate the four-panel plot for your own gLV model starting from an initial
condition `my_IC`, use the command `plot_all_panels(my_p, my_s, my_IC)`.

Please email ewj@physics.ucsb.edu with questions or bugs.
