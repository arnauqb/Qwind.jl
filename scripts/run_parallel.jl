using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind

model = Model("paper2/plot_nn.yaml");

iterations_dict = Dict();
do_iteration!(model, iterations_dict, it_num=1);

do_iteration!(model, iterations_dict, it_num=2);
