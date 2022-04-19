                %%%%%%%%%% ========== DispersionDiagram ========== %%%%%%%%%%
% =============================================================================================
% Author: Adriano M. Goto
% Date: 18/01/2021
% ---------------------------------------------------------------------------------------------
% Code overview: Code for compute the wavenumber for higher order wave modes 
% ---------------------------------------------------------------------------------------------
%
% =============================================================================================

%%%%% ---------- Cleaning ---------- %%%%%
    clear;
    close all;
    clc

%%%%% ---------- Expansion chamber data ---------- %%%%%    
    Data_Silencer;

%%%%% ---------- Parameters of numerical solution of the characteristic equation ---------- %%%%%     
    Nmax = 20;      % Maximuum number of iterations        
    delta = 1e-2;   % Secant method step 
    tol = 1e-14;    % Error tolerance    
    Nmode = 2;      % Number of extra wave modes
    Nz = 1000;      % Number of step over the impedance value
    
%%%%% ---------- Finding roots ---------- %%%%%       
    alpha = 0.1:0.01:100;
    krB = 0.1:1:1000;
    Fa = besselj(1,alpha);
    Fb = Aux_EigenFunction(rm,rc,2*pi*max(freq)/co,0,krB);
    [Nra,indra] = Aux_VerifyRoots(Fa);
    [Nrb,indrb] = Aux_VerifyRoots(Fb);
    [alpha_n,krB] = Aux_BesselJZeros(rm,rc,2*pi*max(freq)/co,0,alpha,krB,Nmax,tol,...
        delta,[indra(1:Nmode+1);indrb(1:Nmode+1)] );
    [krB, Fzeros] = Aux_SolveEigenEquation(rho,co,eta,dh,sigma,t,rm,rc,...
    freq,Nmax,Nz,delta,tol,krB);
    alpha_n = [0;alpha_n(1:end-1)];


%%%%% ---------- Frequency loop ---------- %%%%% 
    kxA = zeros(Nmode+1,length(freq));
    kxB = kxA;
    krA_ecm = alpha_n/rc;
    for cont=1:length(freq)
        % Wavenumbers
        ko = 2*pi*freq(cont)/co;
        ko = ko*(1-1i*eta/2);
        kxA(:,cont) = Aux_WavenumberComp(ko,krA_ecm).';
        kxB(:,cont) = Aux_WavenumberComp(ko,krB(:,cont)).';
  
    end

    
%%%%% ---------- Plotting ---------- %%%%
    %%% ---------- Wavenumbers curves for ECM ---------- %%%    
        figure(1);
        subplot(211)
        plot(freq,real(kxA*Lc),'k','linewidth',2);
        hold on;
        plot(freq,imag(kxA*Lc),'b','linewidth',2);
        grid on;
        xlabel('Frequency [Hz]'); ylabel('\Im[k_{II,n}L_c]  \Re[k_{II,n}L_c]');
    
    %%% ---------- Wavenumbers curves for MPCM ---------- %%%
        subplot(212)
        plot(freq,real(kxB*Lc),'k','linewidth',2);
        hold on;
        plot(freq,imag(kxB*Lc),'b','linewidth',2);
        grid on;
        xlabel('Frequency [Hz]'); ylabel('\Im[k_{II,n}L_c]  \Re[k_{II,n}L_c]');
   
    