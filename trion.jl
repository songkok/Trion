include("Hamiltonian.jl")
include("nonlinearity.jl")

#=== Model parameters =====================#
m1=0.38/7.61994776 #(1st conduction band mass [eV Ang^2])
m2=0.38/7.61994776 #(2nd conduction band mass [eV Ang^2] *for trion only)
mh=0.44/7.61994776 #(Valance band mass [eV Ang^2])
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
#W= Vintra1, Vintra2 (intralayer interactions in layer 1,2) 
#W= Vinter (interlayer interaction) 
#W=VK (monolayer Keldysh interaction)
#define your own potental: W=V(q)

#=== Exiton/Trion Energy & Wavefunction =====================#
#trion=[spectrum([Vintra1 Vintra1 Vintra1],alpha,Q,N,N0,l0) spectrum([Vintra1 Vinter Vinter],alpha,Q,N,N0,l0) spectrum([Vinter Vintra1 Vinter],alpha,Q,N,N0,l0) spectrum([Vinter Vinter Vintra1],alpha,Q,N,N0,l0) spectrum([Vinter Vinter Vintra2],alpha,Q,N,N0,l0)  spectrum([Vinter Vintra2 Vinter],alpha,Q,N,N0,l0)  spectrum([ Vintra2 Vinter Vinter],alpha,Q,N,N0,l0) spectrum([ Vintra2 Vintra2 Vintra2],alpha,Q,N,N0,l0)]
#exciton=[spectrum([Vintra1],alpha,Q,N,N0,l0) spectrum([Vinter],alpha,Q,N,N0,l0) spectrum([Vintra2],alpha,Q,N,N0,l0)]
exciton=spectrum([VK],alpha,Q,N,N0,9)
trion=spectrum([VK VK VK],alpha,Q,N,N0,9)

#=== Nonlinearity calculation ===============================#
#exciton
gX(Theta1,exciton,N,3)
#trion
gT(Theta1,trion,N,3)
