function [] = Changing_Infill_and_G_on_RegressionRate(Lower_Infill, Upper_Infill, Lower_G, Upper_G, Num_Steps)

% Infill_G_and_RhoFuel_VS_RegressionRate(0.1, 0.1, 700, 250, 1000)

Infill_StepSize = (Upper_Infill - Lower_Infill) / Num_Steps;

G_StepSize = (Upper_G - Lower_G) / Num_Steps;

%Rho_Liq_StepSize = (Upper_Rho_liquid - Lower_Rho_liquid) / Num_Steps;

x_array = [];

y_array = [];

for i = 0 : Num_Steps

    Regress_Rate = RegressionRateCalculation ((Infill_StepSize*i + Lower_Infill), (G_StepSize*i + Lower_G));

    x_array = [x_array, i];

    y_array = [y_array, Regress_Rate];


end 
figure;
plot(x_array, y_array, 'b-');
hold on;
plot(x_array, y_array, 'bo','MarkerSize', 2, 'MarkerFaceColor', 'b');
%yticks(0:0.5:1); 
xlabel('Step Number');
ylabel('Calculated Regression Rate (mm/s)');
title(sprintf('Regression Rate Change\n %% Infill from [%g to %g]\n G from [%g to %g]',Lower_Infill, Upper_Infill, Lower_G, Upper_G));