%Возвращает полный потенциал решетки для отображения на графике
function U = fullpotential(type)
global N_GaAs
global N_AlAs
global U_left
global N_period
global U_elstatic
global N_surface
global N_well
if strcmp('gamma', type)
    %Это гамма электроны
    U_GaAs = 0;
    U_AlAs = 1;    
else
    U_GaAs = 0.3;
    U_AlAs = 0;   
end
U = [];
Uc = U_left*ones(1, N_surface);
U_elstatic_part = U_elstatic(1:N_surface);
U = [U Uc + U_elstatic_part];
%Счетчик текущего шага
N_count = N_surface;
for i=1:N_period/2
    Uc = U_GaAs*ones(1, N_GaAs);
    U_elstatic_part = U_elstatic(N_count + 1:N_count + N_GaAs);
    N_count = N_count + N_GaAs;
    U = [U Uc + U_elstatic_part];
    Uc = U_AlAs*ones(1, N_AlAs);
    U_elstatic_part = U_elstatic(N_count + 1:N_count + N_AlAs);
    U = [U Uc + U_elstatic_part];
    N_count = N_count + N_AlAs;
end
Uc = U_GaAs*ones(1, N_well);
U_elstatic_part = U_elstatic(N_count + 1:N_count + N_well);
U = [U Uc + U_elstatic_part];
U = symm_array(U);


