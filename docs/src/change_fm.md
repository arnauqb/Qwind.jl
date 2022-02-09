# Using a different force multiplier

In this section we explain how to change the algorithm that computes the force multiplier. The computation needs to take 2 arguments: effective optical depth `t` and ionisation parameter `xi`.



We first need to import `Qwind` and the `compute_force_multiplier` function



```julia
using Qwind
import Qwind: compute_force_mulitplier
```



We now need to define a new `type` to dispatch our calculation of the force multiplier to.

```julia
struct FMNewFlag <: FMInterpolationFlag end
```

and then we can define our new method:

```julia
function compute_force_multiplier(t, ionization_parameter, mode::FMNewFlag)
    return 1.0
end
```

Now we need to specify the new mode in the parameters that are passed to the model. For now the way to do it is:

```julia
parameters = Parameters("../configs/config_example.yaml")
parameters = change_parameter(parameters, :fm_interp_method_flag, FMNewFlag()) # For now new modes cannot be specified in the .yaml files
```

and then we just initialise and run the model

```julia
model = Model(parameters)
iterations_dict = run!(model)
```

And to make sure it worked we can plot the force multiplier for one streamline:

```julia
using Plots
integrator = iterations_dict[2]["integrators"][10]
fm = integrator.p.data[:fm]
plot(fm)
```

