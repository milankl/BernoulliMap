# Generalised Bernoulli Map with various number formats

Revisiting BM Boghosian, PV Coveney, and H Wang. *A new pathology in the Simulation of Chaotic Dynamical Systems on Digital Computers*, **Adv. Theory Simul. 2019**, 2, 1900125, DOI: [10.1002/adts.201900125](https://dx.doi.org/10.1002/adts.201900125)

![beta2](plots/beta2.png?raw=true "Bernoulli map with beta=2")

**Fig. 1** The generalised Bernoulli map fᵦ for β=2 with different number formats. (a) Float32, (b) Float16, (c) Posit32 and (d) LogFixPoint16. The Bernoulli map f₂ does not introduce arithmetic rounding errors in (a-c), such that stochastic rounding has no impact. The logfix subtraction in f₂ introduces rounding errors in (d) that prevent the stalling at 0 in (a-c).

![invariant measures](plots/inv_measures.png?raw=true "Invariant measures of the Bernoulli map")

**Fig. 2** The invariant measures of the generalised Bernoulli map for β = (a) 32, (b) 43, (c) 54 and (d) 65 calculated with different number formats and rounding modes: 64 and 32-bit floating-point numbers (Float64,Float32), 32-bit posits and Float32 with stochastic rounding. The invariant measures are obtained from long integrations of the Bernoulli map from many different initial conditions x₀ in [0,1). Histograms of x with bin width 0.025 and the analytical invariant measures are normalised.

![error](plots/error_invariant.png?raw=true "Error in invariant measures")

**Fig. 3** Error in the invariant measures for the Bernoulli map. The invariant measures h*(x) simulated with different number formats are quantified from the histograms in Fig. 2. The errors are calculated with respect to the analytical invariant measures h(x). The standard deviations of the error distributions are denoted with 𝜎. For Float16 and LogFixPoint16 parts of the error distribution is outside of the axis limits.

![bifurcation](plots/bifurcation.png?raw=true "Bifurcation")

**Fig. 4** Bifurcation of h(x) with changing β as simulated with various number formats. (a) Analytical bifurcation from Hofbauer's h(x), (b-f) Bifurcation obtained from the histograms of x simulated with various number formats: (b) Float64, (c) Float32, (d) Posit32, (e) Float32 + stochastic rounding, (f) Float16.

