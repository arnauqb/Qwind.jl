export searchsortednearest, countsignchanges

function searchsortednearest(a,x)
   idx = searchsortedfirst(a,x)
   if (idx==1); return idx; end
   if (idx>length(a)); return length(a); end
   if (a[idx]==x); return idx; end
   if (abs(a[idx]-x) < abs(a[idx-1]-x))
      return idx
   else
      return idx-1
   end
end

function countsignchanges(array::Vector{Float64})
   counter = 0 
   if length(array) == 0
      return 0
   end
   current_sign = sign(array[1])
   for elem in array
      if sign(elem) != current_sign
         current_sign = sign(elem)
         counter += 1
      end
   end
   return counter
end
