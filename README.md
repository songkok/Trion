# Trion
This is a Julia numerical package for calculating the spectrum and wavefunction of exciton and trion in two-dimensional materials.

## Tutorial
This package includes: 
- "Trion.jl", for calculating the exciton/trion bound state properties. 
- "Nonlinearity.jl", for evaluating the nonlinearity arise from exciton-exction and trion-trion interactions. 

### One exciton/trion boundstate calculation ("Trion.jl")

To calculate the exciton bound state, we call 'spectrum' in "Trion.jl" as follow
```julia
exciton = spectrum([W],[me mh],alpha,Q,N,N0,L0)
# W is the screening dielectric constant in 2D (see later text)
# me = conduction band mass, mh = valence band mass [ in the unit of free electron mass]
# alpha: 1 = ground state, 2,3,4,... = excited states
# Q = [Qx Qy] is the exciton momentum.
# N is the cutoff parameter for the basis set.
# N0 and L0 are the initial guess for the bound state calculation. 
# Typically, N0 < 3 (initial size of the basis set) and L0 is about the Bohr radius.
```
In this package, the 2D potential is defined as

$$ V(q) = \frac{1}{\epsilon_q} \frac{2\pi e^2}{q}\quad\quad \text{  (cgs unit)}$$

where $\epsilon_q$ is a 2D momentum-dependent screening dielectric constant (dimensionless). 
We note that, $\epsilon_q=1$ in vacuum. For the commonly used Keldysh potential, $\epsilon_q=1+r_\ast q$ where $r_\ast$ is the screening parameter. 

For input in the code,  
$$W(q)=\frac{1}{\epsilon_q}.$$

Once the calculation is done. One can retrive the bound state properties by the object 'exciton' as follow
```julia
exciton.energy  # boundstate energy [ eV ]
exciton.r       # electron and hole average distance [sqrt(<r^2>), Angstrom] 
exciton.A       # exciton wavefunction expanded in the Hermite functions basis [ [nx ny] => C_{nx, ny} ]
```
The expansion of the exciton wavefunction (in momentum space) is
$$\psi_X(k_x,k_y)=\sum_{0\leq n_x+n_y \leq N}C_{n_x,n_y} \varphi_{n_x}(k_x\lambda)\varphi_{n_y}(k_y\lambda)$$
where $\varphi_n(k)=\sqrt{\frac{\lambda}{\pi^{1/2} 2^n n!}}H_n(k)e^{-\frac{1}{2}k^2}$ with Hermite polynomial $H_n(k)$ and optimal basis length $\lambda$ that contains in
```julia
exciton.lambda  # optimal basis length [Angstrom]
```

Similarly, for calculating Trion bound state. We use the following
```julia
trion = spectrum([W12 W1h W2h],[m1 m2 mh],alpha,Q,N,N0,9.0)
# m1, m2 = masses for the two particles with SAME charges. 
# mh = mass for the particle with charges different from m1, m2 particles.  [in the unit of free electron mass]
# W12 = interaction between the particles with mass m1 and mass m2
# W1h = interaction between the particles with mass m1 and mass mh
# W2h = interaction between the particles with mass m2 and mass mh
```
The expansion of the trion wavefunction in momentum ($k$) space is
$$\psi_T(k_{1x},k_{1y},k_{2x},k_{2y})=\sum_{0\leq n_{1x}+n_{1y}+n_{2x}+n_{2y} \leq N}C_{n_{1x},n_{1y},n_{2x},n_{2y}} \varphi_{n_{1x}}(k_{1x}\lambda)\varphi_{n_{1y}}(k_{1y}\lambda) \varphi_{n_{2x}}(k_{2x}\lambda')\varphi_{n_{2y}}(k_{2y}\lambda').$$
We note that $\lambda=\lambda'$ if $m_1=m_2$.  

### Nonlinearity calculation ("Nonlinearity.jl")
```julia
gX(W,exciton,Ncut) # exciton
gT(W,trion,Ncut)   # trion
# W = e-e. e-h, h-h interaction (e=electron, h=hole)
# 'exciton' and 'trion' are the object that are returned by 'spectrum' function.
# Ncut is the number of basis function that will be use in nonlinearity calculation.
```
We note that the code currently only support the nonlinearity calculation for monolayer. Also, the nonlinearity calculation for the trion wavefunction with $\lambda\neq\lambda'$ is not supported too. 

### Example: exciton and trion in TMDC monolayer

```julia
include("Trion.jl")
include("Nonlinearity.jl")
#=== Model parameters =====================#
m1=0.38 #(1st conduction band mass [per free electron mass])
m2=0.38 #(2nd conduction band mass [per free electron mass] *for trion only)
mh=0.44 #(Valance band mass [per free electron mass])
epsilon=1          #(substrate dieletric constant)
r0=40.0/epsilon    #(screen length for monolayer Keldysh potential [Ang])
Q=[0 0]            #(Total momentum [Ang^(-1)])

r1=40.0            #(screen length for Layer 1 [Ang])
r2=45.0            #(screen length for Layer 2 [Ang])
L=6.48             #(Layer 1 and 2 interlayer distance [Ang])

N0=3    # basis size (small) for optimizing the basis length
N=9     # basis size for convergent calculation
l0=[10] # initial guess of the basis length [Ang]. input [lambda1 lambda2] will active the optimization with 2 lengths. 
alpha=1 # 1 = ground state, 2,3,4,... = excited states

#V=[W(e1<->e2), W(e1<->h), W(e2<->h)] (trion)
#V=[W(e<->h)] (exciton)
#where The default potential are:
# W=VKel (monolayer Keldysh interaction)
VKel(q)=1/(epsilon*(1+r0*q))
# W=Vintra1, Vintra2 (intralayer interactions in layer 1,2) 
Vintra1(q)=(1+r2*q*(1-exp(-2*L*q)))/((1+r1*q)*(1+r2*q)-r1*r2*q^2*exp(-2*L*q))
Vintra2(q)=(1+r1*q*(1-exp(-2*L*q)))/((1+r1*q)*(1+r2*q)-r1*r2*q^2*exp(-2*L*q))
# W= Vinter (interlayer interaction)
Vinter(q)=exp(-q*L)/((1+r1*q)*(1+r2*q)-r1*r2*q^2*exp(-2*L*q))

#=== Exiton/Trion Energy & Wavefunction =====================#
exciton=spectrum([VKel],[m1 mh],alpha,Q,N,N0,9.0)
trion=spectrum([VKel VKel VKel],[m1 m2 mh],alpha,Q,N,N0,9.0)


#=== Nonlinearity calculation ===============================#
#exciton
 gX(VKel,exciton,5)
#trion
 gT(VKel,trion,10)

```




## References:

1. [KW Song, S. Chiavazzo, I. A. Shelykh, O. Kyriienko, Attractive trion-polariton nonlinearity due to Coulomb scattering](https://arxiv.org/abs/2204.00594)
    
1. [KW Song, S. Chiavazzo, I. A. Shelykh, O. Kyriienko, Theory for Coulomb scattering of trions in 2D materials](https://arxiv.org/abs/2207.02660)

