min_t = 0; max_t = 100;  % want to end up w/ 1 second at 204800samp/s
Fs = 204800;
dt = (max_t - min_t)/Fs;
tstep = min_t:dt:max_t - dt;
IC = randn(1,3);
[t, xyz] = ode23(@lorenz_ode, tstep, IC);
t = t/100;