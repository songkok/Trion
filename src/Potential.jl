
"""
```Julia
W(V::function, a::Float, lambda::Float, Nmax::Int)
```

Calculating the matrix elements for the potential for Nonlinearity calculation

Inputs:

V    = Interaction between particles (see definition of 2D dielectric constant in README.md )
a    = the factor in the exponent, exp(-a*q^2), this exponential factor appear in the Hermite basis expansion (see Ref[2] in README.md).
Nmax = Cutoff parameter that determines the size of the basis set ( see fTind )
lambda = Basis length need to be optimized by bound state energy ( see README.md )

Returns:

An interacting matrix = (4*Nmax+1)x(4*Nmax+1) matrix{Float}

```
"""

# Keldysh potential matrix elements (Hermite basis)
function W(V, a, lambda, Nmax)
   Hv=zeros(4*Nmax+1,4*Nmax+1)
   for bx in 0:(4*Nmax)
      for by in bx:(4*Nmax)
         Hv[bx+1,by+1]=((bx%2==0 && by%2==0) ? 2*pi*2*quadgk(q -> V(q/lambda)*exp(-a*q^2+(bx+by)*log(q)+log(beta((1+bx)/2,(1+by)/2))), 0, Inf )[1] : 0)
         Hv[by+1,bx+1]=Hv[bx+1,by+1]
      end
   end
   return Hv    
end

"""
```Julia
Wp(V::function, a::Float, Nmax::Int)
```

Calculating the matrix elements for the potential for general purposes.

Inputs:

V    = Interaction between particles (see definition of 2D dielectric constant in README.md )
a    = the factor in the exponent, exp(-a*q^2), this exponential factor appear in the Hermite basis expansion (see Ref[2] in README.md).
Nmax = Cutoff parameter that determines the size of the basis set ( see fTind )

Returns:

An interacting matrix = (4*Nmax+1)x(4*Nmax+1) matrix{Float}

```
"""
function Wp(V, a, Nmax)
   Hv=zeros(4*Nmax+1,4*Nmax+1)
   for bx in 0:(4*Nmax)
      for by in bx:(4*Nmax)
         Hv[bx+1,by+1]=((bx%2==0 && by%2==0) ? 2*pi*2*quadgk(q -> V(q)*exp(-a*q^2+(bx+by)*log(q)+log(beta((1+bx)/2,(1+by)/2))), 0, Inf )[1] : 0)
         Hv[by+1,bx+1]=Hv[bx+1,by+1]
      end
   end
   return Hv    
end


