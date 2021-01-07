# Generalised Bernoulli Map with various number formats

Revisiting BM Boghosian, PV Coveney, and H Wang. *A new pathology in the Simulation of Chaotic Dynamical Systems on Digital Computers*, **Adv. Theory Simul. 2019**, 2, 1900125, DOI: [10.1002/adts.201900125](https://dx.doi.org/10.1002/adts.201900125)

M Kloewer, 2021.
The corresponding notebook can be found [here](https://github.com/milankl/BernoulliMap/blob/master/src/Bernoulli_map.ipynb).

## The generalised Bernoulli map and its invariant measures

The generalised Bernoulli map is a one-variable chaotic system

```
xⱼ₊₁ = fᵦ(xⱼ) = βxⱼ mod 1
```
with a dense attractor for many choices of β>1. Its invariant measures are (Hofbauer, 1978)

```
hᵦ(x) = C∑ᴺⱼ₌₀ β ⁻¹ θ(1ⱼ - x)
```
with `xⱼ = fʲᵦ(xⱼ)` and `1ⱼ = fʲᵦ(1)`. The Heaviside function is `θ(x)` and `C` a normalisation constant.
This series is evaluated with `N→∞`, but in practice `N=100` sufficies to calculate the invariant measure
at much higher precision than can be obtained from the numerical solution of the Bernoulli map.

## The Bernoulli map with β=2

For β=2 the Bernoulli map is poorly represented with most binary arithmetics, such as floating-point number or 
[posits](https://github.com/milankl/SoftPosit.jl).
The multiplication with β acts as a left-bitshift in binary formats with linearly-spaced significand. 
The `mod 1` operation maps `x ∈ [1,2)` to `[0,1)`, which is for floats and posits an error-free subtraction with 1 as
for every `x ∈ [1,2)` the corresponding `x-1` is perfectly representable in the same number number format.
Consequently, an initial condition `x₀` at finite precision (i.e. `0`s exist implicitly for the unrepresented significant bits)
will yield `xₙ = 0` after `n` time steps as `βx` sets one signficant bit per iteration to zero, starting from the least significant.
Analytically, any rational initial condition `x₀` will show the same behaviour, but irrational `x₀` will remain chaotic in time.
As floats or posits only represent rational numbers, their attractors will eventually collapse to `x=0` too.

This effect is shown in Fig. 1 for floats (32-bit single-precision Float32 and 16-bit half-precision Float16) and posits (32-bit Posit32).
[Stochastic rounding](https://github.com/milankl/StochasticRounding.jl) as an alternative to the deterministic 
round to nearest (as in the IEEE-754 standard) has no impact here for Float32 or Float16, as all operations in the Bernoulli map are rounding error-free.

![beta2](plots/beta2.png?raw=true "Bernoulli map with beta=2")

**Fig. 1** The generalised Bernoulli map fᵦ for β=2 with different number formats. (a) Float32, (b) Float16, (c) Posit32 and (d) LogFixPoint16 starting from random initial conditions in [0,1). The Bernoulli map f₂ does not introduce arithmetic rounding errors in (a-c), such that stochastic rounding has no impact. The logfix subtraction in f₂ introduces rounding errors in (d) that prevent the stalling at 0 in (a-c).



![invariant measures](plots/inv_measures.png?raw=true "Invariant measures of the Bernoulli map")

**Fig. 2** The invariant measures of the generalised Bernoulli map for β = (a) 32, (b) 43, (c) 54 and (d) 65 calculated with different number formats and rounding modes: 64 and 32-bit floating-point numbers (Float64,Float32), 32-bit posits and Float32 with stochastic rounding. The invariant measures are obtained from long integrations of the Bernoulli map from many different initial conditions x₀ in [0,1). Histograms of x with bin width 0.025 and the analytical invariant measures are normalised.

![error](plots/error_invariant.png?raw=true "Error in invariant measures")

**Fig. 3** Error in the invariant measures for the Bernoulli map. The invariant measures h*(x) simulated with different number formats are quantified from the histograms in Fig. 2. The errors are calculated with respect to the analytical invariant measures h(x). The standard deviations of the error distributions are denoted with 𝜎. For Float16 and LogFixPoint16 parts of the error distribution is outside of the axis limits.

![bifurcation](plots/bifurcation.png?raw=true "Bifurcation")

**Fig. 4** Bifurcation of h(x) with changing β as simulated with various number formats. (a) Analytical bifurcation from Hofbauer's h(x), (b-f) Bifurcation obtained from the histograms of x simulated with various number formats: (b) Float64, (c) Float32, (d) Posit32, (e) Float32 + stochastic rounding, (f) Float16. Normlisation is applied with max(hᵦ(x)).

