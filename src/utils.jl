export searchsorted_nearest, searchsorted_first, countsignchanges, remove_close_elements

function searchsorted_nearest(a,x)
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

function searchsorted_first(a,x, direction=1)
   idx = searchsortedfirst(a,x)
   if (idx>length(a)); return length(a); end
   if a[idx] == x; return idx; end
   if (idx==1); return 1; end
   return idx - (direction + 1) / 2
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

function remove_close_elements(args...; digits=4)
    ret = [[array[1] for array in args]]
    ret = hcat(ret...)
    for i in 2:length(args[1])
        noinsert = 0
        toinsert = zeros(length(args))
        for (j, array) in enumerate(args)
            element = trunc(array[i], digits=digits)
            if element in ret[j,:]
                noinsert += 1
            end
            toinsert[j] = element
        end
        if noinsert < length(args)
            ret = hcat(ret, toinsert)
        end
    end
    return [ret[i,:] for i in 1:length(args)]
end
