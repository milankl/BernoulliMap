# Hofbauer invariant measure h with C=1
heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))
hofbauer(x,β,N=100) = heaviside(1-x) + [β^-j for j in 1:N]'*heaviside.(bernoulli_map(1e0,β,N)[2:end] .- x)