using SpecialFunctions
using QuadGK
using LinearAlgebra
using Optim

struct boundstate
   energy
   r
   lambda
   A
   order
end

F( l, m, n)=sum([ sum([ (l+m-2*(r+s))==n ? sqrt(factorial(n)/(2^(2*r+l)*factorial(m)))*binomial(m,Int((l+m-n)/2)-r)*factorial(l)/(factorial(r)*factorial((Int((l-m+n)/2)-r))) : 0 for s in 0:minimum([m,l-2r])]) for r in 0:Int(floor(l/2))])

function fTind(Nmax)
   index=[]
   for i in 0:Nmax 
      for j in 0:Nmax
         for k in 0:Nmax 
            for l in 0:Nmax
               if i+j+k+l<=Nmax
                  push!(index,[i j k l])
               end
            end 
         end
      end 
   end
   return index  
end

function fXind(Nmax)
   index=[]
   for i in 0:Nmax 
      for j in 0:Nmax
         if i+j<=Nmax
            push!(index,[i j])
         end
      end 
   end
   return index  
end

# Keldysh potential matrix elements (Hermite basis)
function W(V, a, lambda, Nmax)
   Hv=zeros(3*Nmax+1,3*Nmax+1)
   for bx in 0:(3*Nmax)
      for by in bx:(3*Nmax)
         #Hv[bx+1,by+1]=((bx%2==0 && by%2==0) ? 2*pi*2*beta((1+bx)/2,(1+by)/2)*1/lambda*quadgk(q -> V(q/lambda)*q^(bx+by)*exp(-a*q^2), 0, Inf )[1] : 0)
         Hv[bx+1,by+1]=((bx%2==0 && by%2==0) ? 2*pi*2*quadgk(q -> V(q/lambda)*exp(-a*q^2+(bx+by)*log(q)+log(beta((1+bx)/2,(1+by)/2))), 0, Inf )[1] : 0)
         Hv[by+1,bx+1]=Hv[bx+1,by+1]
      end
   end
   return Hv    
end

function Wp(V, a, Nmax)
   Hv=zeros(3*Nmax+1,3*Nmax+1)
   for bx in 0:(3*Nmax)
      for by in bx:(3*Nmax)
         #Hv[bx+1,by+1]=((bx%2==0 && by%2==0) ? 2*pi*2*beta((1+bx)/2,(1+by)/2)*1/lambda*quadgk(q -> V(q/lambda)*q^(bx+by)*exp(-a*q^2), 0, Inf )[1] : 0)
         Hv[bx+1,by+1]=((bx%2==0 && by%2==0) ? 2*pi*2*quadgk(q -> V(q)*exp(-a*q^2+(bx+by)*log(q)+log(beta((1+bx)/2,(1+by)/2))), 0, Inf )[1] : 0)
         Hv[by+1,bx+1]=Hv[bx+1,by+1]
      end
   end
   return Hv    
end

function VK(bx,by,lambda) 
   global r0
   lambdar0=lambda/r0
   return 2*pi*2*beta((1+bx)/2,(1+by)/2)*(-lambdar0)^(bx+by)*(lambdar0/2*exp(-(lambdar0/2)^2)* (pi*erfi(lambdar0/2)-expinti(( lambdar0/2)^2))-sum([(-lambdar0/2)^(-j)*gamma((1+j)/2) for j in 0:(bx+by-1)]))
end

function HT(Ind, Nmax, V12, V13, V23, m1, m2, mh, Qx, Qy, lambda1, lambda2)
   Nind=length(Ind)
   H=zeros(Nind,Nind)
   k0=13.605692*2*0.529177/(2*pi)^2
   t1=(1/lambda1)
   t2=(1/lambda2)
   lbar=sqrt(lambda1^2+lambda2^2)
   l1=lambda1/lbar
   l2=lambda2/lbar
   KE(t,mc,Q,m,n)=t*(t*(1/mc+1/mh)/2*F(2,m,n)-(Q/mh)*F(1,m,n))
   g(m,n,s,lambda)=(lambda)^(n-s)*(-lambda)^(m-s)*binomial(n,s)*binomial(m,s)*(2^(s-(n+m)/2)*factorial(s)/(sqrt(factorial(m))*sqrt(factorial(n))))#*exp(s*log(2)+log(factorial(s))-1/2*((n+m)*log(2)+log(factorial(m)+log(factorial(n)))))
   W12=W(V12,0.25,lbar,Nmax)
   W13=W(V13,0.25,lambda1,Nmax)
   W23=W(V23,0.25,lambda2,Nmax)
   for i in 1:Nind
      for j in i:Nind
         n=Ind[i]
         m=Ind[j]
         #kinetic terms
         h0=((n==m ? (Qx^2+Qy^2)/(2*mh) : 0) 
         + (n[1:3]==m[1:3] ? KE(t2,m2,Qy,m[4],n[4]) : 0) 
         + (n[1:2]==m[1:2] && n[4]==m[4] ? KE(t2,m2,Qx,m[3],n[3]) : 0) 
         + (n[3:4]==m[3:4] && n[1]==m[1] ? KE(t1,m1,Qy,m[2],n[2]) : 0) 
         + (n[2:4]==m[2:4] ? KE(t1,m1,Qx,m[1],n[1]) : 0)
         + (n[1]==m[1] && n[3]==m[3] ? t1*t2/mh*F(1,m[2],n[2])*F(1,m[4],n[4]) : 0 )
         + (n[2]==m[2] && n[4]==m[4] ? t1*t2/mh*F(1,m[1],n[1])*F(1,m[3],n[3]) : 0 ) )
         #interacting terms
         sm=map(minimum, [vcat(n,m)[:,k] for k in 1:length(n)])
         b=m+n
         hint = sum([
            sum([
               sum([
                  sum([
                     ((b[1]+b[3])%2==0 && (b[2]+b[4])%2==0 ? 1/lbar*W12[b[1]+b[3]-2*(s1+s3)+1,b[2]+b[4]-2*(s2+s4)+1]*g(m[1],n[1],s1,l1)*g(m[2],n[2],s2,l1)*g(n[3],m[3],s3,l2)*g(n[4],m[4],s4,l2) : 0) 
                  for s1 in 0:sm[1]]) 
               for s2 in 0:sm[2]])
            for s3 in 0:sm[3]])
         for s4 in 0:sm[4]])

         hint += sum([sum([-(n[1]==m[1] && n[2]==m[2] && b[3]%2==0 && b[4]%2==0 ? 1/lambda2*W23[b[3]-2*s3+1,b[4]-2*s4+1]*g(m[3],n[3],s3,1)*g(m[4],n[4],s4,1) : 0) for s3 in 0:sm[3]]) for s4 in 0:sm[4]])
         hint += sum([sum([-(n[3]==m[3] && n[4]==m[4] && b[1]%2==0 && b[2]%2==0 ? 1/lambda1*W13[b[1]-2*s1+1,b[2]-2*s2+1]*g(m[1],n[1],s1,1)*g(m[2],n[2],s2,1) : 0) for s1 in 0:sm[1]]) for s2 in 0:sm[2]])
         
      H[i,j]=h0+k0*hint
      H[j,i]=h0+k0*hint
      end
   end
   return H
end

function HX(Ind, Namx, V, m1, mh, Qx, Qy, lambda)
   Nind=length(Ind)
   H=zeros(Nind,Nind)
   k0=13.605692*2*0.529177/(2*pi)^2
   t=(1/lambda)
   KE(t,mc,Q,m,n)=t*(t*(1/mc+1/mh)/2*F(2,m,n)-(Q/mh)*F(1,m,n))
   g(m,n,s,l)=(l)^(n-s)*(-l)^(m-s)*binomial(n,s)*binomial(m,s)*(2^(s-(n+m)/2)*factorial(s)/(sqrt(factorial(m))*sqrt(factorial(n))))
   V1=W(V,0.25,lambda,Namx)
   for i in 1:Nind
      for j in i:Nind
         n=Ind[i]
         m=Ind[j]
         h0=((n==m ? (Qx^2+Qy^2)/(2*mh) : 0) 
         + (n[1]==m[1] ? KE(t,m1,Qy,m[2],n[2]) : 0) 
         + (n[2]==m[2] ? KE(t,m1,Qx,m[1],n[1]) : 0) )

         sm=map(minimum, [vcat(n,m)[:,k] for k in 1:length(n)])
         b=m+n
         hint=sum([sum([-(b[1]%2==0 && b[2]%2==0 ? 1/lambda*V1[b[1]-2*s1+1,b[2]-2*s2+1]*g(m[1],n[1],s1,1)*g(m[2],n[2],s2,1) : 0) for s1 in 0:sm[1]]) for s2 in 0:sm[2]])
   
      H[i,j]=h0+k0*hint
      H[j,i]=h0+k0*hint
      end
   end
   return H
end

function spectrum(V, mass, n,Q,N,N0,l0)
   if length(V)==3
      m1, m2, mh = mass/7.61994776
      ET(ind0,l1,l2)=eigvals(HT(ind0,N0,V[1],V[2],V[3],m1, m2, mh,Q[1],Q[2],l1,l2))[n]
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
      EX(ind0,l)=eigvals(HX(ind0,N0,V[1],me, mh,Q[1],Q[2],l))[n]
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
   return boundstate(result.values[1],r,lambda,A,[ O[1] for O in wavefunc])
end
