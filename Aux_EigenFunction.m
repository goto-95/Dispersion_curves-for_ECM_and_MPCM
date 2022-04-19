function [Fb] = Aux_EigenFunction(rm,rc,ko,zeta,krB)
    %   ===============================================================================
    %%%%%%%%%% ========== Aux_EigenFunction ========== %%%%%%%%%% 
    %%%%% ----- Description ----- %%%%%
    %   Function to calculate the value of the charateristic equation for 
    %   the expansion chamber domain with a disspative micro-perforated layer  
    %
    %%%%% ----- Arguments ----- %%%%%
    % rm     -> Internal micro-perforated duct diameter [m]
    % rc     -> Outer duct diameter [m]
    % ko     -> Total air wavenumber [1/m]
    % zeta   -> Specific acoustic impedance
    % krB    -> Radial wave mode for the expansion chamber with a
    %           micro-perforated layer
    %
    %%%%% ----- Output ----- %%%%%
    % Fb -> Transcendetal characteristic equation for dissipative
    %       micro-perforated layers
    %  
    % ===============================================================================
    
    aux1 = besselj(0,krB*rm) + 1i*zeta*(krB/ko).*besselj(1,krB*rm);
    aux2 = besselj(1,krB*rm).*bessely(1,krB*rc) - besselj(1,krB*rc).*bessely(1,krB*rm);
    aux3 = besselj(0,krB*rm).*bessely(1,krB*rc) - besselj(1,krB*rc).*bessely(0,krB*rm);
    
    Fb = aux1.*aux2 - besselj(1,krB*rm).*aux3;
end
    
    
    
    
    
    
    
    
    
    
    
    