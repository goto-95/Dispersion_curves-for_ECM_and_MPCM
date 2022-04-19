function [alpha_n, beta_n] = Aux_BesselJZeros(rm,rc,ko,zeta,alpha,beta,Nmax,tol,delta,indroots)
    %   ===============================================================================
    %%%%%%%%%% ========== EigenFunction ========== %%%%%%%%%% 
    %%%%% ----- Description ----- %%%%%
    %   Function for solve the charateristic equation for the conventional
    %   expansion chamber muffler and for the micro-perforated chamber
    %   muffler
    %
    %%%%% ----- Inputs ----- %%%%%
    % rm       -> Internal micro-perforated duct diameter [m]
    % rc       -> Expansion chamber diameter [m]
    % ko       -> Total air wavenumber [1/m]
    % zeta     -> Specific acoustic impedance
    % alpha    -> Radial eignemodes for the expansion chamber domain
    % beta     -> Radial eigenmodes for the micro-perforated chamber domain
    % Nmax     -> Maximuum number of iterations
    % tol      -> Error tolerance
    % delta    -> Step of iteration on modified secant method
    % indroots -> Roots index of the alpha and beta vectors
    %
    %%%%% ----- Answer ----- %%%%%
    % alpha_n  -> Eigenvalues axial wavenumber for the expansion chamber
    % beta_n   -> Eigenvalues axial wavenumber for the micro-perforated chamber
    % ===============================================================================
    
    ind1 = indroots(1,:);
    ind2 = indroots(2,:);
    Nroots = length(ind1);

    alpha_n =  zeros(Nroots,1);
    beta_n = alpha_n;
    for cont=1:Nroots  
        [alpha_n(cont),~] = Aux_RootsBessel(alpha(ind1(cont)),Nmax,tol,delta);
        [beta_n(cont),~] = Aux_RootsEigen(rm,rc,ko,zeta,beta(ind2(cont)),Nmax,tol,delta);
    end
end
    