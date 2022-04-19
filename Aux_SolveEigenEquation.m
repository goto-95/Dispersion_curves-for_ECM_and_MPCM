function [krB_n, Fzeros] = Aux_SolveEigenEquation(rho,co,eta,dh,sigma,t,rm,rc,...
    freq,Nmax,Nz,delta,tol,krB_n)
    %   ===============================================================================
    %%%%%%%%%% ========== Aux_SolveEigenFunction ========== %%%%%%%%%% 
    %%%%% ----- Descrição ----- %%%%%
    %   Function for solve the charateristic equation 
    %
    %%%%% ----- Arguments ----- %%%%%
    % rho   -> Air density [kg/m^3] 
    % co    -> Sound speed in air [m/s]
    % eta   -> Air viscosity [Pa*s]
    % dh    -> Diameter of the hole [m]
    % sigma -> Equivalent porosity
    % t     -> Perforated layer thickness [m]  
    % rm    -> Inner perforated duct diameter [m]
    % rc    -> Outer duct diameter [m]
    % freq  -> Frequency vector [Hz]
    % Nmax  -> Maximuum number of iterations
    % Nz    -> Number of steps through the fake-complex-function process (zeta/Nz)
    % delta -> Step of iteration on modified secant method 
    % tol   -> Error tolerance
    %
    %%%%% ----- Answer ----- %%%%%
    % krB_n  -> Radial wavenumber 
    % Fzeros -> Value of the function at the root candidate (For inspection)
    % ===============================================================================
    
    % Computing roots for z==0 e f=max
    f_max = max(freq);
    Nroots = length(krB_n);
    roots =  krB_n;
    
    % Computing roots for z~=0 e f=fmax
    roots_prev = roots;
    omega_max = 2*pi*f_max;
    ko_max = omega_max/co;
    [zeta] = Aux_AcousticImpedanceMPD(rho,co,eta,dh,sigma,t,omega_max);
    for n=0:Nz-1
        % Fiber wavenumber for R~=0
        aux = (n+1)/Nz;
        for cont=1:Nroots  
            [roots(cont), ~] = Aux_RootsEigen(rm,rc,ko_max,aux*zeta,...
                roots_prev(cont),Nmax,tol,delta);
            roots_prev(cont) = roots(cont);
        end
    end
    
    % Computing roots for others frequencies (decrementing)  
    krB_n = zeros(Nroots,length(freq));
    Fzeros = zeros(Nroots,length(freq));
    
    krB_n(:,end) = roots;
    Fzeros(:,end) = Aux_EigenFunction(rm,rc,ko_max,zeta,krB_n(:,end));
    for n=length(freq)-1:-1:1
        omega = 2*pi*freq(n);
        ko = omega/co;
        
        [zeta] = Aux_AcousticImpedanceMPD(rho,co,eta,dh,sigma,t,omega);
        for cont=1:Nroots
            [krB_n(cont,n), ~] = Aux_RootsEigen(rm,rc,ko,zeta,krB_n(cont,n+1),Nmax,tol,delta);
        end
        Fzeros(:,n) = Aux_EigenFunction(rm,rc,ko,zeta,krB_n(:,n));
    end
end 