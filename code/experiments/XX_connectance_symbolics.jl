using Symbolics

# Probability of sharing species
# Probability of rewiring
@variables p q

# Connectance of the two networks, c2 = a*c, 0 <= c <= c2 <= 1, a >= 1
@variables a c

# Links unique to N1 due to turnover
turnover1 = c - p^2*c

# Links unique to N2 due to turnover
turnover2 = a*c - p^2*c*a

# Total unique links due to turnover
turnover = turnover1 + turnover2

# Shareable links
shareable = p^2*c

# Common shared links
common = (1-q) * c * p^2

# Rewired shared links
rewired1 = p^2*c-common
rewired2 = p^2*a*c-common
rewired = rewired1 + rewired2

# Partitions as expressions
A = simplify(common)
U = simplify(turnover)
S = simplify(rewired)
