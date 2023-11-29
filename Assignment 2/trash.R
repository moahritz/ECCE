# set prarameters

S0 = 100
mu = 0.1
sigma = 0.2
T = 1 #final date expressed in years
dt = 1/252 # daily frequency
nsteps = floor(T/dt) 
dW = sqrt(dt)*rnorm(nsteps)
S_euler = numeric(nsteps)
S_euler[1] = 50
r = 0.05
#path simulation using the euler approximation

for (i in 2:nsteps){
  S_euler[i] = S_euler[i-1] + r * S_euler[i-1] *dt + sigma * S_euler[i-1] * dW[i-1]
}


W = cumsum(dW)
S_exact = S0 * exp((r - 0.5 * sigma^2) * dt + sigma * W)

S_exact[252]

S0 * exp((r - 0.5 * sigma^2) + sigma * W[252])
