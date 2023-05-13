using SpecialFunctions
using QuadGK

include("Hamiltonian.jl")

function Theta1(a,bx,by,lambda)
    global r0
   lambdar0=lambda/r0
    return 2*pi*2*beta((1+bx)/2,(1+by)/2)*lambda/2*
        (-lambdar0)^(bx+by)*(lambdar0/2*exp(-a*(lambdar0)^2)*
        (pi*erfi(sqrt(a)*lambdar0)-expinti(a*( lambdar0)^2))-1/sqrt(a)*
        sum([(-sqrt(a)*lambdar0)^(-j)*gamma((1+j)/2) for j in 0:(bx+by-1)]))
end

function H4(h, s1, s2, t, c)
  m1, m2, n1, n2 = h
  g = sum(h)
  return sqrt(pi/2)*(-1)^(m2+n2-t)*(2.0)^(s1+s2-t)*binomial(m1, s1)*binomial(n1, s1)*binomial(m2, s2)*binomial(n2, s2)*
         exp(logfactorial(s1)+logfactorial(s2)+logfactorial(g-2*(s1+s2))-(logfactorial(t)+logfactorial(g-2*(s1+s2+t))))*(2.0*c)^(g-2*(s1+s2+t)) 
end

Nc(m)=1/sqrt(sqrt(pi)*factorial(m)*2^m)

G0(m,n,s,r)=(-1)^(Int(m-s))*binomial(m, s)*binomial(n, r)
f0(m, n, s) =G0(m,n,s,s)*exp(logfactorial(s)-(m+n-2*s)*log(2)/2-logfactorial(n)/2-logfactorial(m)/2)

js(m1, m2, n1p, n2p, m1p, m2p, n1, n2, s1, s2, s1p, s2p) = ((pi/2)*((-1)^(Int(m2 + n2p + m2p + n2))*2^(s1 + s2 + s1p + s2p)*
       exp(logfactorial(s1)+logfactorial(s2)+logfactorial(s1p)+logfactorial(s2p)))/(2^((m1 + m2 + n1p + n2p + m1p + m2p + n1 + n2)/2 - (s1 + s2 + s1p + s2p)))*
       binomial(m1, s1)*binomial(n1p,s1)*binomial(m2, s2)*binomial(n2p, s2)*binomial(m1p, s1p)*binomial(n1, s1p)*binomial(m2p, s2p)*binomial(n2, s2p))

function J0(m1, m2, n1p, n2p, m1p, m2p, n1, n2) 
    S=sum(m1 + m2 + n1p + n2p + m1p + m2p + n1 + n2)
    return ( S%2 == 0 ? sum([sum([sum([sum([js(m1, m2, n1p, n2p, m1p, m2p, n1, n2, s1, s2, s1p, s2p)*
                        (2^(S/2 - (s1 + s2 + s1p + s2p))*sqrt(pi)/gamma((1-S)/2 + (s1 + s2 + s1p + s2p))) 
         for s1 in 0:minimum([m1, n1p])]) for s2 in 0:minimum([m2, n2p])]) for s1p in 0:minimum([m1p, n1])]) for s2p in 0:minimum([m2p, n2])]) : 0)
end

function J(m1x, m2x, n1px, n2px, m1px, m2px, n1x, n2x, m1y, m2y, n1py, n2py, m1py, m2py, n1y, n2y, gx, gy, Tab, c) 
    Sx = m1x + m2x + n1px + n2px + m1px + m2px + n1x + n2x;
    Sy = m1y + m2y + n1py + n2py + m1py + m2py + n1y + n2y;
    return ((gx + Sx)% 2 == 0 && (gy + Sy)% 2 == 0 ? 
        sum([sum([sum([sum([    
        sum([sum([sum([sum([js(m1x, m2x, n1px, n2px, m1px, m2px, n1x, n2x, s1x, s2x, s1px, s2px)*js(m1y, m2y, n1py, n2py, m1py, m2py, n1y, n2y, s1y, s2y, s1py, s2py)*
            sum([sum([(-1)^(Int( tx + ty))*exp((logfactorial(Sx - 2*(s1x + s2x + s1px + s2px))+logfactorial(Sy - 2*(s1y + s2y + s1py + s2py)))-
                    (logfactorial(tx)+logfactorial(ty)+logfactorial(Sx - 2*(s1x + s2x + s1px + s2px + tx))+logfactorial(Sy - 2*(s1y + s2y + s1py + s2py + ty)))
                    )*c^(Sx + Sy - 2*(s1x + s2x + s1px + s2px - tx + s1y + s2y + s1py + s2py - ty))*
                    Tab[ gx + Sx - 2*(s1x + s2x + s1px + s2px + tx) + 1, gy + Sy -  2*(s1y + s2y + s1py + s2py + ty) + 1]
            for tx in 0:Int(floor(Sx/2 - (s1x + s2x + s1px + s2px))) ]) for ty in 0:Int(floor(Sy/2 - (s1y + s2y + s1py + s2py))) ])
        for s1x in 0:minimum([m1x, n1px])]) for s2x in 0:minimum([m2x, n2px])]) for s1px in 0:minimum([m1px, n1x])]) for s2px in 0:minimum([m2px, n2x])])
        for s1y in 0:minimum([m1y, n1py])]) for s2y in 0:minimum([m2y, n2py])]) for s1py in 0:minimum([m1py, n1y])]) for s2py in 0:minimum([m2py, n2y])]) : 0)
end

function ExI(m, n, mp, np, T12, T38)
   mx, my = m
   nx, ny = n
   mpx, mpy = mp
   npx, npy = np
   hx = [nx, mx, npx, mpx]
   hy = [ny, my, npy, mpy]
   gx = sum(hx)
   gy = sum(hy)
   Imn=[0,0]
   
    if (gx%2 == 0 && gy%2 == 0)
        for spx in 0:minimum([mx, mpx])
            for spy in 0:minimum([my, mpy])                    
                Imn +=[sum([sum([ sum([sum([H4(hx, sx, spx, tx, 1/2)*H4(hy, sy, spy, ty, 1/2)*T12[gx-2*(sx+spx+tx)+1, gy-2*(sy+spy+ty)+1] 
                    for tx in 0:Int(floor(gx/2 -(sx+spx))) ]) for ty in 0:Int(floor(gy/2-(sy+spy))) ]) for sx in 0:minimum([nx, npx]) ]) for sy in 0:minimum([ny, npy]) ]),
                    -sum([sum([sum([sum([(-1)^(nx+ny-rx-ry)*binomial(nx,rx)*binomial(ny,ry)*binomial(npx,rpx)*binomial(npy,rpy)*sum([sum([ sum([sum([ H4([rx, mx, rpx, mpx], sx, spx, tx, 1/4)*H4([ry, my, rpy, mpy], sy, spy, ty, 1/4)*T38[gx-2*(sx+spx+tx)+1, gy-2*(sy+spy+ty)+1] 
                    for tx in 0:Int(floor((mx+rx+mpx+rpx)/2 - (sx+spx))) ]) for ty in 0:Int(floor((my+ry+mpy+rpy)/2 - (sy+spy))) ]) for sx in 0:minimum([rx,rpx]) ]) for sy in 0:minimum([ry,rpy]) ]) for rx in 0:nx ])  for ry in 0:ny]) for rpx in 0:npx]) for rpy in 0:npy]) ]
            end
        end
    end
    return -2*prod(map(Nc, vcat(hx, hy)))*Imn
end

function TrI(m, n, mp, np, T1, T12, T14, T34, T38, T516)
m1x, m1y, m2x, m2y = m
m1px, m1py, m2px, m2py = mp
n1x, n1y, n2x, n2y = n
n1px, n1py, n2px, n2py = np
h1x = [n1x, m1x, n1px, m1px]
h1y = [n1y, m1y, n1py, m1py]
h2x = [n2x, m2x, n2px, m2px]
h2y = [n2y, m2y, n2py, m2py]
hx = [m2x, m1px, n1px, n2x]
hy = [m2y, m1py, n1py, n2y]
g1x = sum(h1x)
g1y = sum(h1y)
g2x = sum(h2x)
g2y = sum(h2y)
gx = sum(hx)
gy = sum(hy)

JRmn=4*sum([sum([sum([sum([sum([sum([sum([sum([
    ((g1x - s1x - r1px - r1x - s1px)%2==0 && (g1y - s1y - r1py - r1y - s1py)%2==0 ?
            (-1)^(Int(n1x + m1px + n1y + m1py - r1x - s1px - r1y - s1py))*
            G0(m1x, n1px, s1x, r1px)*G0(m1y, n1py, s1y, r1py)*G0(m1px, n1x, s1px, r1x)*G0(m1py, n1y, s1py, r1y)*
            J0(s1x, m2x, r1px, n2px, s1px, m2px, r1x, n2x)*J0(s1y, m2y, r1py, n2py, s1py, m2py, r1y, n2y)*
            T12[g1x - s1x - r1px - r1x - s1px + 1, g1y - s1y - r1py  - r1y - s1py + 1] : 0)
        for s1x in 0:m1x]) for r1px in 0:n1px]) for s1y in 0:m1y]) for r1py in 0:n1py]) for r1x in 0:n1x]) for s1px in 0:m1px]) for r1y in 0:n1y]) for s1py in 0:m1py])+
        J(m1x, m2x, n1px, n2px, m1px, m2px, n1x, n2x, m1y, m2y, n1py, n2py, m1py, m2py, n1y, n2y, 0, 0, T14, 1)       
        
IRmn = ((m1x + n1x + m2px + n2px + gx)%2 == 0 && (m1y + n1y + m2py + n2py + gy)%2 == 0 ?
    sum([sum([sum([sum([
    sum([sum([
        sum([sum([sum([sum([
            4*f0(m1x, n1x, s1x)*f0(m1y, n1y, s1y)*f0(m2px, n2px, s2px)*f0(m2py, n2py, s2py)*H4(hx, s2x, s1px, tx, 1/2)*H4(hy, s2y, s1py, ty, 1/2)*
            T1[ m1x + n1x + m2px + n2px + gx - 2*(s1x + s2px + s2x + s1px + tx)+1, m1y + n1y + m2py + n2py  + gy - 2*(s1y + s2py + s2y + s1py + ty)+1]
        for s2px in 0:minimum([m2px, n2px]) ]) for s2py in 0:minimum([m2py, n2py]) ]) for s1x in 0:minimum([m1x, n1x]) ]) for s1y in 0:minimum([m1y, n1y]) ])   
    for tx in 0:Int(floor(gx/2 - s2x - s1px)) ]) for ty in 0:Int(floor(gy/2 - s2y - s1py)) ]) 
    for s2x in 0:minimum([m2x, n1px]) ]) for s2y in 0:minimum([m2y, n1py]) ]) for s1px in 0:minimum([m1px, n2x]) ]) for s1py in 0:minimum([m1py, n2y]) ]) : 0) + 
     (m2px == n2px && m2py == n2py && m1x == n1x && m1y == n1y && (m2x + m1px + n2x + n1px)%2== 0 && (m2y + m1py + n2y + n1py%2) == 0 ?
     sum([sum([sum([sum([
        sum([sum([ H4(hx, s2x, s1px, tx, 1/2)*H4(hy, s2y, s1py, ty, 1/2)*T12[ gx - 2*(s1px + s2x + tx)+1, gy - 2*(s2y + s1py + ty)+1]
        for tx in 0:Int(floor(gx/2 - s2x - s1px)) ]) for ty in 0:Int(floor(gy/2 - s2y - s1py)) ])
    for s2x in 0:minimum([m2x, n1px]) ]) for s2y in 0:minimum([m2y, n1py]) ]) for s1px in 0:minimum([m1px, n2x]) ]) for s1py in 0:minimum([m1py, n2y]) ]) : 0)

JAmn = -2*sum([sum([sum([sum([
            G0(m1x, n1px, s1x, r1px)* G0(m1y, n1py, s1y, r1py)* 
                J(s1x, m2x, r1px, n2px, m1px, m2px, n1x, n2x, s1y, m2y, r1py, n2py, m1py, m2py, n1y, n2y, (m1x + n1px - s1x - r1px), (m1y + n1py - s1y - r1py), T516, -(1/2))
            for s1x in 0:m1x]) for r1px in 0:n1px]) for s1y in 0:m1y]) for r1py in 0:n1py])-  
        2*sum([sum([sum([sum([ 
            G0(m1px, n1x, s1px, r1x)*G0(m1py, n1y, s1py, r1y)*(-1)^(Int(m1py + n1y - s1py - r1y + m1px + n1x - s1px - r1x))*
                J(m1x, m2x, n1px, n2px, s1px, m2px, r1x, n2x, m1y, m2y, n1py, n2py, s1py, m2py, r1y, n2y, (m1px + n1x - s1px - r1x), (m1py + n1y - s1py - r1y), T516, 1/2)
            for r1x in 0:n1x]) for s1px in 0:m1px]) for r1y in 0:n1y]) for s1py in 0:m1py])

IAmn = (-2*(m2px == n2px && m2py == n2py && (m1x + n1x + gx)%2 == 0 && (m1y + n1y + gy)%2 == 0 ?
    sum([sum([
        sum([sum([sum([sum([
            sum([sum([f0(m1x, n1x, s1x)*f0(m1y, n1y, s1y)*H4(hx, s2x, s1px, tx, 1/2)*H4(hy, s2y, s1py, ty, 1/2)* 
                T34[ m1x + n1x  + gx - 2*( s1x + s2x + s1px + tx) + 1, m1y + n1y + gy - 2*(s1y + s2y + s1py + ty) + 1]
            for tx in 0:Int(floor(gx/2 - s2x - s1px)) ]) for ty in 0:Int(floor(gy/2 -  s2y -  s1py)) ])
        for s2x in 0:minimum([m2x, n1px]) ]) for s2y in 0:minimum([m2y, n1py]) ]) for s1px in 0:minimum([m1px, n2x]) ]) for s1py in 0:minimum([m1py, n2y]) ])
    for s1x in 0:minimum([m1x, n1x]) ]) for s1y in 0:minimum([m1y, n1y]) ]) : 0) - 
    2*(m1x == n1x && m1y == n1y && (gx + m2px + n2px)%2 == 0 && (gy + m2py + n2py)%2 == 0 ? 
    sum([sum([
        sum([sum([sum([sum([
            sum([sum([f0(m2px, n2px, s2px)*f0(m2py, n2py, s2py)*H4(hx, s2x, s1px, tx, 1/2)*H4(hy, s2y, s1py, ty, 1/2)*
                T34[ m2px + n2px + gx - 2*(s2px + s2x + s1px + tx) + 1, m2py + n2py + gy - 2*(s2py + s2y + s1py + ty) + 1]
            for tx in 0:Int(floor(gx/2 - s2x - s1px)) ]) for ty in 0:Int(floor(gy/2 -  s2y -  s1py)) ])
        for s2x in 0:minimum([m2x, n1px]) ]) for s2y in 0:minimum([m2y, n1py]) ]) for s1px in 0:minimum([m1px, n2x]) ]) for s1py in 0:minimum([m1py, n2y]) ])
    for s2px in 0:minimum([m2px, n2px]) ]) for s2py in 0:minimum([m2py, n2py]) ]) : 0)
    )

    return (-1)^(Int(g2x+g2y))*prod(map(Nc, vcat(vcat(h1x,h2x), vcat(h1y,h2y))))*([JRmn,JAmn])-prod(map(Nc, vcat(hx, hy)))*[IRmn,IAmn]
end

function gX(V,exciton,N,Ncut)
    k0=13.605692*2*0.529177/(2*pi)^2
    Amn=[exciton.A[i][2] for i in 1:Ncut]
    m=[exciton.A[i][1] for i in 1:Ncut]
    lambda=exciton.lambda[1]
    SumI=[]
    for n1 in 1:Ncut
        for n2 in 1:Ncut
            for n3 in 1:Ncut
                for n4 in 1:Ncut
                    push!(SumI,[n1, n2, n3, n4])
                end
            end
        end
    end
    sort!(SumI,by=sum)
    T12= lambda/2*W(V, 0.5, lambda, N)
    T38= lambda/2*W(V, 3.0/8.0, lambda, N)
    #[repulsive,attractive]
    return  (2*pi)^2*k0/100*sum([Amn[i[3]]*Amn[i[2]]*Amn[i[1]]*Amn[i[4]]*ExI(m[i[3]],m[i[1]],m[i[4]],m[i[2]],T12,T38) for i in SumI]) 
end

function gT(V,trion,N,Ncut)
    k0=13.605692*2*0.529177/(2*pi)^2
    Amn=[trion.A[i][2] for i in 1:Ncut]
    m=[trion.A[i][1] for i in 1:Ncut]
    lambda=trion.lambda[1]
    SumI=[]
    for n1 in 1:Ncut
        for n2 in 1:Ncut
            for n3 in 1:Ncut
                for n4 in 1:Ncut
                    push!(SumI,[n1, n2, n3, n4])
                end
            end
        end
    end
    sort!(SumI,by=sum)
    T1= lambda/2*W(V, 1.0, lambda, N)
    T12= lambda/2*W(V, 0.5, lambda, N)
    T14= lambda/2*W(V, 0.25, lambda, N)
    T34= lambda/2*W(V, 3.0/4.0, lambda, N)
    T38= lambda/2*W(V, 3.0/8.0, lambda, N)
    T516= lambda/2*W(V, 5.0/16.0, lambda, N)
    #[repulsive,attractive]
    return  (2*pi)^2*k0/100*sum([Amn[i[3]]*Amn[i[4]]*Amn[i[1]]*Amn[i[2]]*TrI(m[i[1]],m[i[2]],m[i[3]],m[i[4]], T1, T12, T14, T34, T38, T516) for i in SumI]) 
end


#
