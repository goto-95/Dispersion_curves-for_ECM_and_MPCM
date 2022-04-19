function [kout] = Aux_WavenumberComp(ko,kin)
    %   ===============================================================================
    %%%%%%%%%% ========== Aux_AlphaCoef ========== %%%%%%%%%% 
    %%%%% ----- Description ----- %%%%%
    %   Function to compute correctly the wavenumber for any wave mode type
    %
    %%%%% ----- Input ----- %%%%%
    % ko -> Absolute wavenumber 
    % kr -> Radial wavemode obtained from the characteristic equation
    %
    %%%%% ----- Output ----- %%%%%
    % kz -> Axial wavenumber
    % 
    % ===============================================================================
    
    kout = zeros(1,length(kin));
    for cont=1:length(kin)
        if real(ko)>=real(kin(cont))
            kout(cont) = sqrt(ko^2-kin(cont)^2);
        else
            kout(cont) = sqrt(ko^2-kin(cont)^2);
        end
    end
end
    