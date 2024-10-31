clear
clc

%% Set Parameters
B = 4.7; %- %Blowing parameter -> sourced from Evaluation of the Homologous Series of Normal Alkanes as Hybrid Rocket Fuels
mu_g = 6.5E-6; %Pa-s %Average gaseous viscosity in port -> guess

rho_f = 930; %kg/m^3 %Solid fuel density -> sourced from material properties database
rho_g = 4; %kg/m^3 %Average gaseous density in port -> guess

%Ratio for later use, reasonable estimate in range given by Karabeyoglu
Q_dot_r = 0.15;
Q_dot_c = 1;

G = 700; %kg/m^2-s -> sourced from Dracarys
z = 0.01; %m -> variable parameter

T_g = 2200; %K %Average gas phase temperature -> guess, assumed to be equal to paraffin flame temperature
T_v = 731; %K %Vaporization temperature -> sourced from NIST
T_m = 341; %K %Melting temperature -> sourced from NIST
T_a = 290; %K %Initial temperature -> set parameter

C_s = 2.1; %kJ/kg-C %Solid specific heat -> sourced from material properties database
C_l = 2.1; %kJ/kg-C %Liquid specific heat -> sourced from material properties database

L_v = 1600; %kJ/kg %Heat of vaporization -> sourced from NIST
L_m = 220; %kJ/kg %Heat of fusion -> sourced from NIST

h_m = L_m + C_s*(T_m - T_a); %kJ/kg %Effective heat of melting
h_e = h_m + C_l*(T_v - T_m); %kJ/kg %Effective heat of gas stream

a_ent = 1E-11; %- %Entrainment coefficient -> guess based on Karabeyoglu
alpha_hat = 1; %- %Dynamic pressure exponent -> sourced from Evaluation of the Homologous Series of Normal Alkanes as Hybrid Rocket Fuels
beta_hat = 1; %- %Thickness exponent -> sourced from Evaluation of the Homologous Series of Normal Alkanes as Hybrid Rocket Fuels

%% Initial Values
r_dot = 1;
r_dot_v = 1;
r_dot_ent = 0;

converge = 1;
prev = 1;
i=1;

%% Equations

C_B1 = 2/(2+ 1.25*B^0.75)

C_B2 = (1.25*B^0.75)/(2+ 1.25*B^0.75)

r_dot_cl = (0.03*mu_g^0.2/rho_f)*(1+ Q_dot_r/Q_dot_c)*B*C_B1*G^(0.8)*z^(-0.2)

F_r = 1+ (14.1*rho_g^0.4)/(G^0.8*((T_g/T_v)^0.2))

R_hv = (C_l*(T_v - T_m))/(h_e + L_v)

R_he = h_m / (h_e + L_v)

while converge > 0.0001

    if r_dot_v < 0
        r_dot_v = 0;
    end

    r_dot_v = -(R_he+R_hv*(r_dot_v/r_dot))*r_dot_ent + F_r*(0.03*mu_g^(0.2)/rho_f)*(1+ Q_dot_r/Q_dot_c)*B*((C_B1)/(C_B1+C_B2*(r_dot_v/r_dot_cl)^0.75))*G^(0.8)*z^(-0.2);
    
    r_dot_ent = a_ent * G^(2*alpha_hat)/r_dot^(beta_hat);

    r_dot = r_dot_v + r_dot_ent

    a(i,:) = [r_dot, r_dot_ent, r_dot_v];
    i=i+1;
    converge = abs(r_dot - prev)/abs(prev);
    prev = r_dot;
end