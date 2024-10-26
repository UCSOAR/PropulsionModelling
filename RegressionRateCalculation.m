
function [RegressionRate] = RegressionRateCalculation(Zeta_v, G)
    
    % Check number of input arguments

    if nargin < 1E-10
        Zeta_v = 0.10;
        G = 700;
        %Rho_l = 850;
    end
    
    % RegressionRateCalculation(0.1, 700, 850)
    
    %% Input Variables
    
    %Zeta_v = 0.1; % lattice volume fraction

    h_v_ent = 200; % entrainment vaporized specific enthalpy [J/kg]
    
    h_v = ((230-180)/2 + 180); % vaporized specific enthalpy [J/kg] % between 180 - 230 kJ/kg
   
    Mu_g = 6.5E-5; % dynamic gas viscosity [Pa·s] % From simulation software 7.376E-5
    
    Rho_l = 850; % liquid density [kg/m3] % density of light liquid paraffin as 0.83-0.86 g/mL and heavy liquid paraffin as 0.86-0.89 g/mL
     
    B = 4.7; % Blowing Parameter 7 or 14
    
    C_h = 2; % Stanton number with blowing (from paper) 
    
    C_ho = 2+1.25*B^0.75; % Stanton number without blowing
    
    %G = 700; % oxidizer flux [kg/(m2·s)]
    
    z = 0.065; % Axial location [m] % Colin 5.4.1 0.065 
    
    Q_dot_r = 0.15; % Q_dot_r/Q_dot_C ~ 0.1 to 0.3 % Radiative heat flux [J/(m2·s)]
    
    Q_dot_c = 1; % Convective heat flux [J/(m2·s)]
    
    h_g_wax = 1300*1000; % Gaseous specific enthalpy of wax [J/kg]
    
    h_g_lat = 1300*1000; % Gaseous specific enthalpy of lattice [J/kg]
    
    %u_e = 1.22; % velocity of oxidizer gas parallel to fuel surface [m/s]
    
    %u_b = 1; % velocity of flame parallel to fuel surface [m/s] % Karabeyoglu states that Ue/Ub ~ 1.22
    
    %delta_h = 1000; % enthalpy difference between the flame and the surface,
    
    %h_g = 1800; % heat of gasification
    
    alpha = 1.5; % thermal diffusivity [m2/s] or absorptivity
    
    beta = 1.5; % entrainment liquid layer exponent
    
    T_b_fu = 2220; % Flame zone temperature of lattice augmented fuel K
    
    T_p_fu = 850; % Peak surface temperature of lattice augmented fuel K
   
    delta_h_g_prf = 1600; % Heat of vaporization paraffin [KJ/Kg]
    
    delta_h_g_fu = 2300; % Heat of vaporization of lattice augmented fuel [KJ/Kg] % pure ABS in N2O is 2.3MJ/Kg or 3.23 MJ/kg
    
    delta_h_g_prf_fu = (delta_h_g_prf * (1-Zeta_v) + delta_h_g_fu * (Zeta_v)); % Weighted heat of vaporization of wax + ABS lattice depending on infill percentage

    T_b_prf = 2220; % Flame zone temperature of paraffin [K]
    
    T_p_prf = 518; % Peak surface temperature of paraffin [K] % paraffin temp of vaporization 518 K from NIST
    
    rho_wax = 900; % Density of Wax [kg/m3] % approx density of solid paraffin
    
    rho_lat = 1100; % Density of lattice [kg/m3] % approx density of ABS
    
    a_ent = 8E-14; % entrainment regression rate coefficient [m/s]
    
    %R_dot = 1;
    R_dot_v = 1;
    R_dot_ent = 1;
    
    %% Beginning of Calculations

    % Enthalpy modification term that accounts for entrainment of wax at temperatures below the vaporization temperature,
    % R_he + R_hv(rv/r), has been omitted in favor of the original Rh.
    % equation 5.30
        
    R_h = h_v_ent / h_v ;
    
    % R_rho_hv is a non-dimensional ratio that captures the density and effective heat of gasification ratios for the wax to
    % the lattice equation 5.39
        
    R_rho_hv = (1 - Zeta_v) * (h_g_wax / h_g_lat) * (rho_wax / rho_lat) ;
        
    %Since the lattice heat of gasification differs from that of paraffin, the influence of its addition needs to be
    %corrected for by modifying the blowing parameter (B) equation 5.34
        
    %B = (u_e * delta_h) / (u_b * h_g) ;
        
    % CB is the correction term that is a function of the lattice fraction and regression rate
    % Subscripts b and p denote temperatures associated with the flame zone and the peak surface temperature, respectively. 
    % Subscripts prf and fu refer to neat paraffin and the lattice-augmented
    % fuels. equation 5.36
        
    %C_B = ((T_b_fu - T_p_fu) * delta_h_g_prf) / ((T_b_prf - T_p_prf) * delta_h_g_fu);

    C_B = ((T_b_fu - T_p_fu) * delta_h_g_prf) / ((T_b_prf - T_p_prf) * delta_h_g_prf_fu);
    
    %% Non-Linear Equations

    tolerance = 1E-5; % Set the tolerance level 0.882
    difference_tolerance = inf; % Initialize difference to infinity
    prev_value = 0; % Initialize previous value
    
    %array = [];    
        
    while difference_tolerance >= tolerance
        
        R_dot = R_dot_v + R_dot_ent;
    
        R_dot_ent = (a_ent * (G ^ (2 * alpha))) / (R_dot ^ beta);
        
        R_dot_v = ((((0.03*(Mu_g^0.2))/Rho_l) * (C_h/C_ho) * B * C_B * (G^0.8) * (z^-0.2) * (1 + (Q_dot_r/Q_dot_c))) - (R_dot_ent * (Zeta_v + (R_rho_hv * R_h))))/(Zeta_v + R_rho_hv) ;
        
        
        % Messing with the order of computations changes the outputs:
    
        % Rdot,RdotEnt,RdotV = ok
        % Rdot,RdotV,RdotEnt = not ok (complex)
        % RdotEnt,Rdot,RdotV = not ok (complex)
        % RdotEnt,RdotV,Rdot = ok
        % RdotV,RdotEnt,Rdot = not ok (complex)
        % RdotV,Rdot,RdotEnt = ok
    
        % Order must follow RdotV --> Rdot --> RdotEnt (at any starting point)
        % To work properly and give non-complex solution   
        
        % Setting the current value and checking the weighted difference between them
    
        current_value = R_dot;
    
        difference_tolerance = abs(current_value - prev_value);
    
        prev_value = current_value;
    
        %array(end + 1) = current_value;
    end
    
    %RegressionRate = ['Predicted Regression Rate: ', num2str(current_value*1000), ' mm/s'];
    
    RegressionRate = (current_value*1000);
end

% disp(array);
% disp(RegressionRate)
% plot(array);