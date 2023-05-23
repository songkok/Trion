"""
```Julia
fXind( Nmax::Int )
```
Indexing the basis set for exciton Hamiltonian matrix construction. 

Inputs: 

Nmax = cutoff parameter that determines the size of the basis set

Returns: 

a vector of index = [ [0 0], [1 0], [0 1], .... , [i j] ...., [0 Nmax] ] 
                    where i + j <= Nmax

```
"""
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

"""
```Julia
fTind( Nmax::Int )
```
Indexing the basis set for trion Hamiltonian matrix construction. 

   Nmax = cutoff parameter that determines the size of the basis set

Returns a vector of index: [ [0 0 0 0], [1 0 0 0], [0 1 0 0], .... , [i j k l] ...., [0 0 0 Nmax] ] 
                           where i + j + k + l <= Nmax

```
"""
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

