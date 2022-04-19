function [kroot, value] = Aux_RootsEigen(rm,rc,ko,zeta,kr_prev,Nmax,tol,delta)
    %   ===============================================================================
    %%%%%%%%%% ========== EigenFunction ========== %%%%%%%%%% 
    %%%%% ----- Description ----- %%%%%
    % Computes the eigenroots of the characteristic equation the for expansion chamber 
    % domain with a disspative micro-perforated layer  
    %
    %%%%% ----- Arguments ----- %%%%%
    % rm     -> Internal micro-perforated duct diameter [m]
    % rc     -> Outer duct diameter [m]
    % ko     -> Total air wavenumber [1/m]
    % zeta   -> Specific acoustic impedance
    % Nmax   -> Maximuum number of iterations
    % tol    -> Error tolerance
    % delta  -> Step of iteration on modified secant method 
    %
    %%%%% ----- Answer ----- %%%%%
    % kroots  -> Eigenvalues axial wavenumber 
    % value   -> Value of the function at the root candidate (For inspection)
    % ===============================================================================
    
    Nite = 0;
    value = 0;
    res=1;
    
    while res==1
        aux1 = Aux_EigenFunction(rm,rc,ko,zeta,kr_prev);
        num = (1+delta)*kr_prev - (1-delta)*kr_prev;
        den = Aux_EigenFunction(rm,rc,ko,zeta,(1+delta)*kr_prev) - ...
            Aux_EigenFunction(rm,rc,ko,zeta,(1-delta)*kr_prev);
        kroot = kr_prev - aux1*num/den;
        
        value = abs(Aux_EigenFunction(rm,rc,ko,zeta,kroot));
        Nite = Nite+1;
        kr_prev = kroot;
        if value < tol
            res= 0;
        elseif Nite >= Nmax
            res = 0;
        end
    end
end