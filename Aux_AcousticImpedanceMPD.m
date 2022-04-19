function [z] = Aux_AcousticImpedanceMPD(rho,co,eta,dh,sigma,t,omega)
    %%%%%%%%%% ========== Aux_AcousticImpedanceMPD ========== %%%%%%%%%% 
    %%%%% ----- Overview ----- %%%%%
    %   Function to compute the characteristic acoustic impedance of micro-
    %   perforated layer
    %
    %%%%% ----- Input ----- %%%%%
    % rho   -> Air density [kg/m^3] 
    % co    -> Sound speend in air [m/s]
    % eta   -> Air viscosity [Pa*s]
    % dh    -> Perforation diameter [m]
    % sigma -> Equivalent porosity
    % t     -> Micro-perforated layer thickness
    % omega -> Angular frequency
    %
    %%%%% ----- Output ----- %%%%%
    % z -> Characteristic impedance 
    % =====================================================================
    
    K = dh*sqrt((omega*rho)/(4*eta));
    re_z = ((32*eta*t)/(sigma*rho*co*dh^2))*( sqrt(1+K^2/32) + (sqrt(2)*K*dh)/(32*t));
    imag_z = ((omega*t)/(sigma*co))*( 1 + 1/sqrt(9+0.5*K^2) + 0.85*dh/t );
    z = re_z + 1i*imag_z;
    
end