# Generalised Bernoulli Map with various number formats

Revisiting BM Boghosian, PV Coveney, and H Wang. *A new pathology in the Simulation of Chaotic Dynamical Systems on Digital Computers*, **Adv. Theory Simul. 2019**, 2, 1900125, DOI: [10.1002/adts.201900125](https://dx.doi.org/10.1002/adts.201900125)

M Kloewer, 2021.
The corresponding notebook with details on the numerical implementation can be found [here](https://github.com/milankl/BernoulliMap/blob/master/src/Bernoulli_map.ipynb).

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
round to nearest (as in the IEEE-754 standard) has no impact here for Float32 or Float16, as all operations in the Bernoulli map
are rounding error-free.

![beta2](plots/beta2.png?raw=true "Bernoulli map with beta=2")

**Fig. 1** The generalised Bernoulli map fᵦ for β=2 with different number formats. (a) Float32, (b) Float16, (c) Posit32
and (d) LogFixPoint16 starting from random initial conditions in `[0,1)`. The Bernoulli map f₂ does not introduce arithmetic
rounding errors in (a-c), such that stochastic rounding has no impact. The logfix subtraction in f₂ introduces rounding errors
in (d) that prevent the collapse of the attractor to 0 in (a-c).

Interestingly, the logarithmic-fixed point number format [LogFixPoint16](https://github.com/milankl/StochasticRounding.jl) does not 
show the collapse of the attractor after (at most) `n` time steps, which corresponds to the number of significant/fraction bits in
the respective number format. Logfixs are rounding error-free in multiplication, but not in addition, such that the `mod 1` operation
of the Bernoulli map yields both `0` or `1` in the least signficant bit, preventing the attractor from collapse. Consequently, the
rounding error for logfixs preserves the chaotic property of the analytical Bernoulli map.

## Invariant measures of the Bernoulli map with varying β

The Hofbauer invariant measure `hᵦ(x)` of the Bernoulli map are compared with numerical estimations when simulating the system with various number 
formats. The Bernoulli map is calculated with an ensemble of `Nens = 1000` members (each starting from a different random initial condition)
for `Nsteps = 10,000` time steps, discarding a spin-up of 100 time steps. The sample size for `x` is therefore about `10,000,000`, smaller ensembles
or fewer time steps were not found to yield robust statistics. The invariant measures are estimated as histograms `h*(x)` over `x` with a bin-width of 0.025 and normalised to probability 1, `∫₀¹ h*(x) dx = 1`.

The free parameter β is changed to be 3/2, 4/3, 5/4 or 6/5, such that the multiplication with β in the Bernoulli map is not rounding error-free.
The estimated invariant measures resemble the analytical one from Hofbauer, but differ due to rounding errors. Note that the discrete analytical invariant measures do not coincide with the chosen bins, such that the discrete steps are not correctly represented - an artifact from the methodology.

![invariant measures](plots/inv_measures.png?raw=true "Invariant measures of the Bernoulli map")

**Fig. 2** The invariant measures of the generalised Bernoulli map for β = (a) 3/2, (b) 4/3, (c) 5/4 and (d) 6/5 
calculated with different number formats and rounding modes: 64 and 32-bit floating-point numbers (Float64,Float32),
32-bit posits and Float32 with stochastic rounding. The invariant measures are obtained from long integrations of the
Bernoulli map from many different initial conditions `x₀` in `[0,1)`. Histograms of x with bin width 0.025 and the
analytical invariant measures are normalised.

Other than for β=2, where the attractor entirely collapses, the misrepresentation of the invariant measure here is only slight.
The question arises whether this is improved with higher precision or with other number formats. The error between the
analytical invariant measure and the numerical estimation is shown in Fig. 3. 

![error](plots/error_invariant.png?raw=true "Error in invariant measures")

**Fig. 3** Error in the invariant measures for the Bernoulli map. The invariant measures h*(x) simulated with different number formats are quantified from the histograms in Fig. 2. The errors are calculated with respect to the analytical invariant measures h(x). The standard deviations of the error distributions are denoted with 𝜎. For Float16 and LogFixPoint16 parts of the error distribution is outside of the axis limits.

![bifurcation](plots/bifurcation.png?raw=true "Bifurcation")

**Fig. 4** Bifurcation of h(x) with changing β as simulated with various number formats. (a) Analytical bifurcation from Hofbauer's h(x), (b-f) Bifurcation obtained from the histograms of x simulated with various number formats: (b) Float64, (c) Float32, (d) Posit32, (e) Float32 + stochastic rounding, (f) Float16. Normlisation is applied with max(hᵦ(x)).

