"""
```Julia
boundstate( energy::Float, r::Float, A::Dict, lambda::Float, N::Int, order::Vec )
```

A data type that store the bound state properties.

energy = binding energy [eV]
r      = root mean square distances, sqrt(<(r_e-r_h)^2>) for exciton, 
        [sqrt(<(r_1-r_2)^2>), sqrt(<(r_1-r_h)^2>), sqrt(<(r_2-r_h)^2>)] for trion [Ang]
A      = wavefunction expanded in Hermite basis.
lambda = optimized basis length [Ang]
N      = Cutoff parameter for basis
order  = vector, rank the basis according to their weight in the wavefunction basis expansion

```
"""
struct boundstate
   energy
   r
   A
   lambda
   N
   order
end

"""
```Julia
spectrum(V::Vec{function}, mass::Vec{Float}, alpha::Int, Q::Vec{Float}, N::Int, N0::Int, l0::Float)
```

Diagonalizing the Wannier equation for exciton/trion.

Inputs:
   
V     = [W] for exciton,     [W12, W13, W23] for trion, where W is the q-dependent 2D dielectric constant (see README.md)
mass  = [me mh] for exciton, [m1 m2 mh]      for trion
alpha = 1 is ground state, 2,3,4,... are excited states
Q     = [Qx Qy] Trion total momentum [ Ang^{-1} ]
N     = cutoff parameter of the basis set
N0    = cutoff parameter of the basis set, use for optimizing the basis length lambda (N0~3 is suffiient) 
l0    = Initial guess for basis length, approximate to Bohr radius [Ang] 

Returns: 

exciton::boundstate or trion::boundstate which inlcude the boundstate properties (see struc boundstate)

```
"""

function spectrum(V, mass, alpha, Q, N, N0, l0)
   if length(V)==3
      m1, m2, mh = mass/7.61994776
      ET(ind0,l1,l2)=eigvals(HT(ind0,N0,V[1],V[2],V[3],m1, m2, mh,Q[1],Q[2],l1,l2))[alpha]
      if length(l0)==2
         ET2(l)=ET(fTind(N0),l[1],l[2])
         spec=optimize(ET2,l0,BFGS())
         lambda=[spec.minimizer[1] spec.minimizer[2]]
      else
         ET1(l)=ET(fTind(N0),l,l) 
         spec=optimize(ET1,l0[1],10*l0[1])
         lambda=[spec.minimizer spec.minimizer]
      end
      ind=fTind(N)
      result=eigen(HT(ind,N,V[1],V[2],V[3], m1, m2, mh, Q[1],Q[2],lambda[1],lambda[2]))
   else
      me, mh = mass/7.61994776
      EX(ind0,l)=eigvals(HX(ind0,N0,V[1],me, mh,Q[1],Q[2],l))[alpha]
      EX0(l1)=EX(fXind(N0),l1)
      spec=optimize(EX0,l0[1],10*l0[1]) 
      lambda=[spec.minimizer]
      ind=fXind(N)
      result=eigen(HX(ind,N,V[1],me, mh ,Q[1],Q[2],spec.minimizer))
   end
   A=hcat(ind,result.vectors[:,1])
   wavefunc=[A[i,:] for i in 1:length(ind)]
   sort!(wavefunc, by = x -> -abs(x[2]))
   A=Dict(wavefunc)
   function r2f(A, ind, i, di)
      return sum([ A[n]*((n[i] >= 2 ? -sqrt(n[i]*(n[i]-1)/4)*A[n-di] : 0)+(2*n[i]+1)/2*A[n]-(n[i]+2<=N && sum(n)+2<=N ? sqrt((n[i]+1)*(n[i]+2)/4)*A[n+di] : 0)) for n in ind]) 
   end
   if length(V)==3
      r13=lambda[1]*sqrt(r2f(A,ind,1,[2 0 0 0])+r2f(A,ind,2,[0 2 0 0]))
      r23=lambda[2]*sqrt(r2f(A,ind,3,[0 0 2 0])+r2f(A,ind,4,[0 0 0 2]))
      r=[sqrt(r13^2+r23^2-2*lambda[1]*lambda[2]*sum([A[n]*(-(n[1]>=1 && n[3]>=1 ? sqrt(n[1]*n[3]/2^2)*A[n-[1 0 1 0]] : 0) - (sum(n)+2<=N ? sqrt((n[1]+1)*(n[3]+1)/4)*A[n+[1 0 1 0]] : 0)+
      (n[3]>=1 ? sqrt((n[1]+1)*n[3]/2^2)*A[n+[1 0 -1 0]] : 0)+(n[1]>=1 ? sqrt(n[1]*(n[3]+1)/4)*A[n+[-1 0 1 0]] : 0)+
      -(n[2]>=1 && n[4]>=1 ? sqrt(n[2]*n[4]/2^2)*A[n-[0 1 0 1]] : 0)-( sum(n)+2<=N ? sqrt((n[2]+1)*(n[4]+1)/4)*A[n+[0 1 0 1]] : 0)+
      ( n[4]>=1 ? sqrt((n[2]+1)*n[4]/2^2)*A[n+[0 1 0 -1]] : 0)+( n[2]>=1 ? sqrt(n[2]*(n[4]+1)/4)*A[n+[0 -1 0 1 ]] : 0)) for n in ind])),r13,r23 ]
   else
      r=lambda[1]sqrt(r2f(A,ind,1,[2 0])+r2f(A,ind,2,[0 2]))
   end
   return boundstate(result.values[1],r,A,lambda,N,[ O[1] for O in wavefunc])
end