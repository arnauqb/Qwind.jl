using Distances
using BenchmarkTools

function my_distance(dp, d1, d2)
    ret = d2^2 - ((d1^2 - dp^2 - d2^2) / (2dp))^2
    return sqrt(ret)
end

function true_distance(x1, y1, x2, y2, p1, p2)
    a = abs((y2-y1) * p1 - (x2-x1)*p2 + x2*y1 - y2*x1)
    b = sqrt((y2-y1)^2 + (x2-x1)^2)
    return a/b
end

p1 = 175
p2 = 22

x1 = 168.45030367578497
y1 = 4.593011820974516

x2 = 168.45509110166472
y2 = 4.6009526868356145

dp = evaluate(Euclidean(), [x1,y1], [x2,y2])
d1 = evaluate(Euclidean(), [x1,y1], [p1,p2])
d2 = evaluate(Euclidean(), [x2,y2], [p1,p2])

@test my_distance(dp, d1, d2) ≈ true_distance(x1, y1, x2, y2, p1, p2)
