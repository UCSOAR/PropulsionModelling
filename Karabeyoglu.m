clear
clc

%% Set Parameters

mu_g = 6.5E-6;%Pa-s  %Average gaseous viscosity in port -> guess

%Densities
rho_f = 930;%kg/m^3  %Solid fuel density -> sourced from material properties database
rho_g = 4;%kg/m^3  %Average gaseous density in port -> guess

%Ratio of radiative to convective heat transfer -> reasonable range given by Development and testing of Paraffin-based Hybrid Rocket Fuels
Q_dot_r = 0.15;
Q_dot_c = 1;

%Motor conditions
B = 4.7;  %Blowing parameter -> sourced from Evaluation of the Homologous Series of Normal Alkanes as Hybrid Rocket Fuels
G = 200;%kg/m^2-s  %Oxidizer mass flux -> sourced from Dracarys
z = 0.01;%m  %Position along fuel grain -> variable parameter

%Temperatures
T_g = 1673;%K  %Average gas phase temperature -> guess, assumed to be equal to paraffin flame temperature
T_v = 731;%K  %Vaporization temperature -> sourced from NIST
T_m = 341;%K  %Melting temperature -> sourced from NIST
T_a = 290;%K  %Initial temperature -> set parameter

%Specific heats
C_s = 2.1;%kJ/kg-K  %Solid specific heat -> sourced from material properties database
C_l = 2.1;%kJ/kg-K  %Liquid specific heat -> sourced from material properties database

%Heats and effective heats
L_v = 1600;%kJ/kg  %Heat of vaporization -> sourced from NIST
L_m = 220;%kJ/kg  %Heat of fusion -> sourced from NIST

h_m = L_m + C_s*(T_m - T_a);%kJ/kg  %Effective heat of melting
h_e = h_m + C_l*(T_v - T_m);%kJ/kg  %Effective heat of gas stream

%Empirical entrainment parameters
a_ent = 1E-11;  %Entrainment coefficient -> guess based on Development and testing of Paraffin-based Hybrid Rocket Fuels
alpha_hat = 1;  %Dynamic pressure exponent -> reasonable estimate from Evaluation of the Homologous Series of Normal Alkanes as Hybrid Rocket Fuels
beta_hat = 1;  %Thickness exponent -> reasonable estimate from Evaluation of the Homologous Series of Normal Alkanes as Hybrid Rocket Fuels

%% Initial Values

%Guess values for regression
%For some reason these initial values are the only way solutions converge
r_dot = 1;%m/s  %Total regression rate
r_dot_v = 1;%m/s  %Regression rate due to fuel vaporization
r_dot_ent = 0;%m/s  %Regression rate due to droplet entrainment

prev = 1;%m/s  %Previous regression rate, for convergence purposes

%Convergence values
converged = false;
converge = 1;
converge_change = 0;

%Iteration and divergence counters
i=1;
solution_resets = 0;

%% Equations (all sourced from Development and testing of Paraffin-based Hybrid Rocket Fuels)

%Blowing correction coefficients
C_B1 = 2/(2+ 1.25*B^0.75);
C_B2 = (1.25*B^0.75)/(2+ 1.25*B^0.75);

%Classical regression rate (no entrainment)
r_dot_cl = (0.03*mu_g^0.2/rho_f)*(1+ Q_dot_r/Q_dot_c)*B*C_B1*G^(0.8)*z^(-0.2);%m/s

%Roughness parameter
F_r = 1+ (14.1*rho_g^0.4)/(G^0.8*((T_g/T_v)^0.2));

%Ratios of effective heat of gasifications for entrainment and vaporization
R_hv = (C_l*(T_v - T_m))/(h_e + L_v);
R_he = h_m / (h_e + L_v);

%% Computation loop
while converged == false

    %Divergence correction
    %For some reason only r_dot_v can be corrected; correcting any other variables fail the solution
    if r_dot_v < 0
        r_dot_v = 0;
        solution_resets = solution_resets + 1;
    end

    %Convergence verification
    %Tracks both convergence of r_dot between iterations and the change between the convergences (sometimes r_dot doesn't change by happenstance when the solution hasn't stabilized)
    if converge < 0.001 && (converge_change < .1 || isnan(converge_change)) 
        converged = true;
        disp('Solution converged.')
        
        %Solution sanity check, must be positive and entirely real - exits otherwise
        %Solution error printout
        if imag(r_dot_ent + r_dot_v + r_dot) ~= 0 || r_dot < 0
            disp('Solution error, solution is negative or imaginary.')
            disp([newline,'Solution stats:'])
            disp(['Number of iterations: ',num2str(length(iterations))])
            disp(['Last convergence: ',num2str(converge)])
            disp(['Change in convergence between last 2 iterations: ',num2str(converge_change)])
            disp(['Number of corrected divergences: ',num2str(solution_resets)])     
            break
        end

        %Successful solution printout
        disp(['Computed regression rate (in mm/s): ',num2str(r_dot*1000)])
        disp([newline,'Solution stats:'])
        disp(['Number of iterations: ',num2str(length(iterations))])
        disp(['Last convergence: ',num2str(converge)])
        disp(['Change in convergence between last 2 iterations: ',num2str(converge_change)])
        disp(['Number of corrected divergences: ',num2str(solution_resets)])
    end

    %Program max iteration limit and fail printout
    if i>=300000
        disp(['Solution failed to converge within ', num2str(i), ' iterations.' ])
        disp([newline,'Solution stats:'])
        disp(['Last convergence: ',num2str(converge)])
        disp(['Change in convergence between last 2 iterations: ',num2str(converge_change)])
        disp(['Number of corrected divergences: ',num2str(solution_resets)])

        break
    end

    prev = r_dot;
   
    %Calculating all three regressions
    r_dot_v = -(R_he+R_hv*(r_dot_v/r_dot))*r_dot_ent + F_r*(0.03*mu_g^(0.2)/rho_f)*(1+ Q_dot_r/Q_dot_c)*B*((C_B1)/(C_B1+C_B2*(r_dot_v/r_dot_cl)^0.75))*G^(0.8)*z^(-0.2);
    r_dot_ent = a_ent * G^(2*alpha_hat)/r_dot^(beta_hat);
    r_dot = r_dot_v + r_dot_ent;

    %Array tracking all three regressions to monitor solution stability
    iterations(i,:) = [r_dot, r_dot_ent, r_dot_v];
    i=i+1;
    
    %Saving previous convergence and calculating current convergence
    converge2 = converge;
    converge = abs(r_dot - prev)/abs(prev);

    %Calculating change in convergences
    converge_change = abs(converge - converge2)/abs(converge2);

end
