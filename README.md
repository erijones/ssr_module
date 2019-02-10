Steady State Reduction (SSR)
============================

This repository stores supplemental code for the paper "Steady State Reduction
of generalized Lotka-Volterra systems in the microbiome," by Eric Jones and
Jean Carlson, (to be) published in Phys. Rev. E in 2019. The code in
`steady_state_reduction_example.py` generates Fig. 2 of the paper. Steady State
Reduction (SSR) approximates a high-dimensional gLV system by two-dimensional
(2D) subspaces that each obey gLV dynamics and span a pair of high-dimensional
steady states. Therefore, the process of SSR determines the 2D parameters that
best approximate in-plane dynamics of the original high-dimensional system.

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
5) Plot the separatrix of the SSR-reduced 2D system, saved in `SSR_demo_2.pdf`
6) Plot the in-plane separatrix of the high-dimensional system, saved in
`SSR_demo_2.pdf`

To use SSR yourself, generate your own `Params` class, say `p`, and generate the pair of
steady states you intend to use for SSR, say `ssa` and `ssb`. Then generate new
ssrParams with `s = ssrParams(p, ssa, ssb)`. You can check the validity of SSR
as implemented in your own system with the three `plot_*` commands.

Please email ewj@physics.ucsb.edu with questions or bugs.
