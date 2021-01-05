# Generalised Bernoulli Map with various number formats

![beta2](plots/beta2.png?raw=true "Bernoulli map with beta=2")

**Fig. 1** The generalised Bernoulli map 𝑓_𝛽 for 𝛽=2 with different number formats. (a) Float32, (b) Float16, (c) Posit32 and (d) LogFixPoint16. The Bernoulli map 𝑓_2 does not introduce arithmetic rounding errors in (a-c), such that stochastic rounding has no impact. The logfix subtraction in 𝑓_2 introduces rounding errors in (d) that prevent the stalling at 0 in (a-c).

![invariant measures](plots/inv_measures.png?raw=true "Invariant measures of the Bernoulli map")

**Fig. 2** The invariant measures of the generalised Bernoulli map for β = (a) 32, (b) 43, (c) 54 and (d) 65 calculated with different number formats and rounding modes: 64 and 32-bit floating-point numbers (Float64,Float32), 32-bit posits and Float32 with stochastic rounding. The invariant measures are obtained from long integrations of the Bernoulli map from many different initial conditions 𝑥0 in [0,1). Histograms of 𝑥 with bin width 0.025 and the analytical invariant measures are normalised.

![error](plots/error_invariant.png?raw=true "Error in invariant measures")

**Fig. 3** Error in the invariant measures for the Bernoulli map. The invariant measures ℎ∗(𝑥) simulated with different number formats are quantified from the histograms in Fig. 2. The errors are calculated with respect to the analytical invariant measures ℎ(𝑥). The standard deviations of the error distributions are denoted with 𝜎.