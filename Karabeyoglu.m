function r_dot = Karabeyoglu(G, z)
    %% Set Parameters
    
    if nargin < 1E-10
        z = 0.01;
        G = 500;
    end

    mu_l = 3.067E-3;%Pa-s  %Fuel liquid viscosity
    mu_g = 7.376E-5;%Pa-s  %Average gaseous viscosity in port -> from RPA

    %Densities
    rho_f = 930;%kg/m^3  %Solid fuel density (fuel+matrix average) -> sourced from material properties database
    rho_l = 750;%kg/m^3  %Liquid fuel density (fuel+matrix average)
    rho_g = 2.252;%kg/m^3  %Average gaseous density in port -> from RPA
    
    %Ratio of radiative to convective heat transfer -> reasonable range given by Development and testing of Paraffin-based Hybrid Rocket Fuels
    Q_dot_r = 0.15;
    Q_dot_c = 1;
    
    %Motor conditions
    B = 4.7;  %Blowing parameter -> sourced from Evaluation of the Homologous Series of Normal Alkanes as Hybrid Rocket Fuels
    %G = 700;%kg/m^2-s  %Oxidizer mass flux -> sourced from Dracarys
    %z = 0.01;%m  %Position along fuel grain -> variable parameter
    
    %Temperatures
    T_g = 2378.5;%K  %Average gas phase temperature -> from RPA
    T_v = 731;%K  %Vaporization temperature -> sourced from NIST
    T_m = 341;%K  %Melting temperature -> sourced from NIST
    T_a = 290;%K  %Initial temperature -> set parameter
    
    %Specific heats
    C_s = 2100;%kJ/kg-K  %Solid specific heat -> sourced from material properties database
    C_l = 2100;%kJ/kg-K  %Liquid specific heat -> sourced from material properties database
    
    %Heats and effective heats
    L_v = 160000;%kJ/kg  %Heat of vaporization -> sourced from NIST
    L_m = 22000;%kJ/kg  %Heat of fusion -> sourced from NIST
    
    h_m = L_m + C_s*(T_m - T_a);%kJ/kg  %Effective heat of melting
    h_e = h_m + C_l*(T_v - T_m);%kJ/kg  %Effective heat of gas stream
    
    %Empirical entrainment parameters
    a_ent = 1.4E-10;  %Entrainment coefficient -> guess based on Development and testing of Paraffin-based Hybrid Rocket Fuels
    %a_ent = 6.6554E-13*(0.8*rho_l)/(mu_l*rho_f);
    alpha_hat = 1;  %Dynamic pressure exponent -> reasonable estimate from Evaluation of the Homologous Series of Normal Alkanes as Hybrid Rocket Fuels
    beta_hat = 1;  %Thickness exponent -> reasonable estimate from Evaluation of the Homologous Series of Normal Alkanes as Hybrid Rocket Fuels
    
    %% Initial Values
    
    %Guess values for regression
    r_dot_0 =  [0.001;%m/s  %Total regression rate
                0.001;%m/s  %Regression rate due to fuel vaporization
                0];%m/s  %Regression rate due to droplet entrainment
    
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
    
    %% Computation
    %Calculating all three regressions

    equations = @(r_dot) [
        F_r*(0.03*mu_g^(0.2)/rho_f)*(1+ Q_dot_r/Q_dot_c)*B*((C_B1)/(C_B1+C_B2*(r_dot(2)/r_dot_cl)^0.75))*G^(0.8)*z^(-0.2)-r_dot(2)-(R_he+R_hv*(r_dot(2)/r_dot(1)))*r_dot(3);
        a_ent*G^(2*alpha_hat)/r_dot(1)^(beta_hat)-r_dot(3);
        r_dot(2)+r_dot(3)-r_dot(1)];
        
    options = optimoptions('lsqnonlin', 'FunctionTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10);
    solution = lsqnonlin(equations, r_dot_0, [0,0,0],[1,1,1], options);

    r_dot = solution(1);

    disp(solution)
    
end
