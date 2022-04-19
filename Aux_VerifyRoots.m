function [Nroots,indroots] = Aux_VerifyRoots(funvec)
    % =====================================================================
    %%%%%%%%%% ========== Aux_VerifyRoots ========== %%%%%%%%%% 
    %%%%% ----- Overview ----- %%%%%
    %   Function to identify the position of potential  roots accoroding to
    %   the change of function signal
    %
    %%%%% ----- Input ----- %%%%%
    % funvec -> Function vector
    %
    %%%%% ----- Output ----- %%%%%
    % Nroots   -> Number of finded roots
    % indroots -> Roots index according to the function vector provided
    % =====================================================================
    
    Nroots = 0;
    indroots = zeros(1,1);
    for cont=1:length(funvec)-1
        f1 = sign(funvec(cont));
        f2 = sign(funvec(cont+1));

        if f1 ~= f2
            Nroots = Nroots + 1;
            indroots(Nroots) = cont;
        end
    end
end
    