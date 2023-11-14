function a_moon = moon_pert(time, r)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DA TOGLIERE DALLA CARTELLA ONE DRIVE, MOON PERTURBATION INSERITA NELL ODE_2BP_PERTURBED %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_moon = astroConstants(20);

day = time/(60*60*24);

[r_moon] = ephMoon(day);
r_moon = r_moon';

%%

r_sc_moon = r_moon - r;
r_sc_moon_norm = norm(r_sc_moon);
r_moon_norm = norm(r_moon);


a_moon = mu_moon * ( r_sc_moon ./ r_sc_moon_norm^3 - r_moon ./ r_moon_norm^3);


end

