%Задает начальные данные
clearvars
global fi1
global t
t = [];
%Общий потенциал решетки
global lattice_U
%lattice_U = [];
%Количество точек для апроксиммации одного нанометра (надо чтоб было кратно
% 10)
N_apr = 10;
%Количество полных периодов решетки должно быть кратно двум
global N_period
N_period = 40;
%Толщина слоя GaAs в нм
l_GaAs = 2.3;
%Количество точек для апроксимации слоя GaAs
global N_GaAs
N_GaAs = l_GaAs * N_apr;
%Толщина слоя AlAs в нм
l_AlAs = 1.1;
%Количество точек для апроксимации слоя AlAs
global N_AlAs
N_AlAs = l_AlAs * N_apr;
%Толщина ямы в нм
l_well = 46;
global N_well
N_well = l_well*N_apr;
%Длина боковых отступов в нм
l_surface = 10;
%Количество точек для аппроксимации боковых отступов
global N_surface
N_surface = l_surface*N_apr;
%Длина всей решетки в нм
global l_lattice
l_lattice = 2*(N_period/2*(l_GaAs + l_AlAs) + l_well/2 + l_surface);
%Всего шагов на решетку
global N
N = l_lattice*N_apr;
%Начальный потенциал без учета разрыва дна зоны проводимости
%Разрыв зоны проводимости между GaAs и AlAs в эв
global delta_U
delta_U = 1;
%Левая стенка в эв
global U_left
U_left = 3;
%Находим положение дельта слоя примеси. Он будет располагаться в 12 периоде
%в слое GaAs
%Начало слоя
global N_start_delta
N_start_delta = 11*(l_GaAs + l_AlAs)*N_apr + l_surface*N_apr + 10;
global U_elstatic
U_elstatic = zeros(1, N);
nd = zeros(1, N);
for i = N_start_delta + 1:N_start_delta + 3
    nd(i) = 4.2*10^15/(3*1/N);
end
x = linspace(0, 1, N);
% [psi, E] = reshenie();
% psi1 = psi(1, :);
% n1 = psi1.^2*8.4*10^15;
% n = nd - n1;
% n = symm_array(n);
% U_elstatic = puasson(n);
for i=1:80
    %Решаем Шредингера
    [psi, E] = reshenie();
    E_analize(2*i-1) = E(1);
    E_analize(2*i) = E(2);    
    psi1 = psi(1, :);
    psi2 = psi(2, :);
    %Задаем распределение концентрации в подзонах
    n1 = psi1.^2*4.4*10^15;
    n2 = psi2.^2*4*10^15;
    %Полное распределение концентрации
    n = nd - n1 - n2;
    n = symm_array(n);
    U_elstatic = puasson(n);
end
plot(x, U_elstatic);



