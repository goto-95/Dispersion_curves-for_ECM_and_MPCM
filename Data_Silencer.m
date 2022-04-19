            %%%%% ---------- Data_Silencer ---------- %%%%%
% =============================================================================== % 
% Acoustical silencer parameters 
% =============================================================================== %

%%%%% ----- Frequency range ----- %%%%%   
    fo   = 0.1;       % Initial frequency [Hz]
    fmax = 3e3;       % Final frequency [Hz]
    df = 2.5;         % Frequency step [Hz]
    freq = fo:df:fmax; 
    omega = 2*pi*freq;

%%%%% ----- Properties ----- %%%%%
    eta = 1.8e-5;   % Air viscosity [Pa*s]
    rho = 1.204;    % Air density [kg/m^3]
    co  = 343.3;    % Sound speed in air [m/s]
       
%%%%% ----- Acoustical silencer parameters ----- %%%%%  
    Lc = 100e-3;     % Chamber length [m]
    Ld = 68e-3;      % Inlet/outlet duct length [m]
    Dd = 60e-3;      % Inlet/outlet duct diameter [m]
    rd = Dd/2;       % Inlet duct radius [m]
    Dm = 60e-3;      % Internal micro-perforated duct diameter [m] 
    rm = Dm/2;       % Internal micro-perforated duct radius [m] 
    Dc = 150e-3;     % Expansion chamber diamter [m]
    rc = Dc/2;       % Expansion chamber radius [m]
    Sc = pi*rc^2;    % Cross-section chamber area [m^2]
    Sd = pi*rd^2;    % Cross-section duct area [m^2]
    SL = 2*pi*rm*Lc; % Lateral internal duct area [m^2]

%%%%% ----- Microperforation parameters ----- %%%%%       
    dh = .80e-3;        % Perforation diameter [m]
    Ap = 0.25*pi*dh^2;  % Perforation area [m^2]
    t = 3e-3;           % Micro-perforated duct thickness [m]
    Nf = 36*10;         % Number of perforations
    sigma = Nf*Ap/SL;   % Equivalent porosity
    


