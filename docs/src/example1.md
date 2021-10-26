# Example 1: Mdot dependence with radius

This example illustrates how one can change the normalised accretion rate as a function of radius.

Let us load our example model

```julia
using Qwind
model = Model("configs/config_example.yaml")
```

since we specified `mdot=0.5`, the mdot grid is flat at 0.5:

```julia
using PyPlot
plt.plot(model.rad.disk_grid, model.rad.mdot_grid)
```

Let us now run the model for a couple of iterations

```julia
iterations_dict = Dict()
run!(model, iterations_dict, start_it=1, n_iterations=2)
```

We can plot the last iteration,

```julia
streamlines = iterations_dict[2]["streamlines"]
fig, ax = plt.subplots()
for sl in streamlines
    ax.plot(sl.r, sl.z)
end
```

At this point, imagine we would like to change our mdot grid to have a different value with radius. Let's suppose we want it to go from 0.2 to 0.5, we can do:

```julia
n_disk = length(model.rad.disk_grid)
new_mdot_grid = range(0.2, 0.5, length=n_disk)
model.rad.mdot_grid .= new_mdot_grid
```

And check it has updated properly

```julia
plt.plot(model.rad.disk_grid, model.rad.mdot_grid)
```

Now we can run one more iteration

```julia
run!(model, iterations_dict, start_it=3, n_iterations=1)
```

and check the new results

```julia
streamlines = iterations_dict[3]["streamlines"]
fig, ax = plt.subplots()
for sl in streamlines
    ax.plot(sl.r, sl.z)
end
```
