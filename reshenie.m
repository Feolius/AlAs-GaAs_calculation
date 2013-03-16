%Задает граничные условия и вызывает решение Шредингера для каждого
%значения энергии
function [psi, E_ans]=reshenie(type)
global N_GaAs
global N_AlAs
global U_left
global N
global N_period
global U_elstatic
global N_surface
global m
global N_well
global U
global j
%Масса свободного электрона в г
m_e = 9.109*10^(-28);
%Шаг по расстоянию
h = 1/N;
%Задаем разрыв зоны проводимости и эффективную массу в слоях
if strcmp('gamma', type) > 0
    %Это гамма электроны
    U_GaAs = 0;
    U_AlAs = 1;
    m_GaAs = 0.068 * m_e;
    m_AlAs = 0.15 * m_e;
else
    U_GaAs = 0.3;
    U_AlAs = 0;
    m_GaAs = 1.3 * m_e;
    m_AlAs = 1.1 * m_e;
end
%Задаем максимум сканирования по энергии
if U_GaAs ~= 0
    E_max = U_GaAs;
else
    E_max = U_AlAs;
end;
%Число шагов по энергии и величина шага
nE = 10000;
dE = E_max/nE;
%Счетчики итераций и решений
j = 1;
n = 1;
E = 0;
% if strcmp('X', type)
% hold on;
% x = linspace(0, 1, N);
% axis([0 N -5 5]);
% psi_temp1 = zeros(1, N);
% h1=plot(x, psi_temp1,'g');
% set(h1,'EraseMode','xor','MarkerSize',18);
% end;
while E < 1
    
    %Сюда складируем решения в течении одной итерации
    psi_temp1 = zeros(1, N);
    psi_temp2 = zeros(1, N);
    %Задаем граничные условия на значение и производную слева
    psi_temp1(1) = 10^-20;
    psi_temp2(1) = psi_temp1(1)*sqrt(U_left-E);
    %Решаем Шредингера для каждого слоя п/п. На каждом шаге вычисляются
    %куски потенциала от разрыва зоны проводимости и от носителей. В каждом
    %из таких кусков одна лишняя точка в начале для сшивки. Затем на
    %каждом слое решается Шредингер с граничными условиями, обеспечивающими
    %сшивку вф между слоями.
    
    Uc = U_left*ones(1, N_surface);
    U_elstatic_part = U_elstatic(1:N_surface);   
    U = Uc + U_elstatic_part;
    m = m_AlAs;    
    [psi1, psi2] = runge_kutt(psi_temp1(1), psi_temp2(1), N_surface - 1, E);       
    psi_temp1(2:N_surface) = psi1;
    psi_temp2(2:N_surface) = psi2;
    %Счетчик текущего шага
    N_count = N_surface;    
    %Решаем Шредингера для каждого слоя п/п, но для числа
    %перидов равного половине + 1
    %Число периодов, для которого будем решать
    N_period_sol = N_period/2;
    for i = 1:N_period_sol
        %Решение для слоя GaAs
        Uc = U_GaAs*ones(1, N_GaAs + 1);
        U_elstatic_part = U_elstatic(N_count:N_count + N_GaAs);
        U = Uc + U_elstatic_part;
        m = m_GaAs;
        [psi1, psi2] = runge_kutt(psi_temp1(N_count), psi_temp2(N_count), N_GaAs, E);
        psi_temp1(N_count + 1:N_count + N_GaAs) = psi1;
        psi_temp2(N_count + 1:N_count + N_GaAs) = psi2;
        N_count = N_count + N_GaAs;
        %Решение для слоя AlAs
        Uc = U_AlAs*ones(1, N_AlAs + 1);
        U_elstatic_part = U_elstatic(N_count:N_count + N_AlAs);
        U = Uc + U_elstatic_part;
        m = m_AlAs;
        [psi1, psi2] = runge_kutt(psi_temp1(N_count), psi_temp2(N_count), N_AlAs, E);
        psi_temp1(N_count + 1:N_count + N_AlAs) = psi1;
        psi_temp2(N_count + 1:N_count + N_AlAs) = psi2;
        N_count = N_count + N_AlAs;
    end
    %Решаем подзадачу для ямы
    Uc = U_GaAs*ones(1, N_well + 1);
    U_elstatic_part = U_elstatic(N_count:N_count + N_well);
    U = Uc + U_elstatic_part;
    m = m_GaAs;
    [psi1, psi2] = runge_kutt(psi_temp1(N_count), psi_temp2(N_count), N_well, E);
    psi_temp1(N_count + 1:N_count + N_well) = psi1;
    psi_temp2(N_count + 1:N_count + N_well) = psi2;    
    if j > 1
        %Если значение функции или производной поменяло знак в середине, значит это
        %решение
        psians = [];
        if w1*psi_temp1(fix(N/2) + 1) < 0
            psians = symm_array(psi_temp1);
        elseif w2*psi_temp2(fix(N/2) + 1) < 0
            psians = symm_array(psi_temp1);
        end
        [s1, s2] = size(psians);
        if s2 ~= 0
            square = 0;
            for g = 1:N
                square = square + psians(g)^2*h;
            end
            psi(n, :) = psians/square^0.5;
            E_ans(n) = E;
            if n == 2 && strcmp(type, 'gamma')
                break
            elseif n == 2 && strcmp(type, 'X')
                break
            end            
            n = n + 1;
        end
    end
%     if strcmp('X', type)
%      set(gca,'YLim',[min(psi_temp1) max(psi_temp1)])
%          set(h1,'YData',psi_temp1);
%          pause(0.00001);
%     end;
    w1 = psi_temp1(fix(N/2) + 1);
    w2 = psi_temp2(fix(N/2) + 1);    
    j = j + 1;
    E = E + dE;
    
end


