% run_poopsi
cd ~/Research/oopsi/pop-oopsi/code/JoVo/
clear

if 0
    load('../../data/JoVo/sim1_data.mat')   % load data
else
    pop_sim
end
smc = pop_oopsi(Cell,V);                % do everything
omegahat = plot_poopsi(V,smc,Cell{1}.P.omega)        % plot results