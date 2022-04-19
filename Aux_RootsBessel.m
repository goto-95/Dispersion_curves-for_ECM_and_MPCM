function [kroot, value] = Aux_RootsBessel(kr_prev,Nmax,tol,delta)
    %   ===============================================================================
    %%%%%%%%%% ========== EigenFunction ========== %%%%%%%%%% 
    %%%%% ----- Description ----- %%%%%
    %   Computes the Bessel roots for the expansion chamber domain without
    %   any dissipative layer
    %
    %%%%% ----- Arguments ----- %%%%%
    % kr_prev -> Previous radial roots
    % Nmax    -> Maximuum number of iterations
    % tol     -> Error tolerance
    % delta   -> Step of iteration on modified secant method 
    %
    %%%%% ----- Answer ----- %%%%%
    % kroots  -> Eigenvalues axial wavenumber
    % value   -> Value of the function at the root candidate (For inspection)
    % ===============================================================================
    
    Nite = 0;
    value = 0;
    res=1;
    
    while res==1
        aux1 = besselj(1,kr_prev);
        num = (1+delta)*kr_prev - (1-delta)*kr_prev;
        den = besselj(1,(1+delta)*kr_prev) - besselj(1,(1-delta)*kr_prev);
        kroot = kr_prev - aux1*num/den;
        
        value = abs(besselj(1,kroot));
        Nite = Nite+1;
        kr_prev = kroot;
        if value < tol
            res= 0;
        elseif Nite >= Nmax
            res = 0;
        end
    end
end