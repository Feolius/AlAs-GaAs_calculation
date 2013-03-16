%Задает начальные данные
clearvars
global j
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
l_well = 36;
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
%U_elstatic = zeros(1, N);
nd = zeros(1, N);
for i = N_start_delta + 1:N_start_delta + 3
    nd(i) = 4.65*10^15/(3*1/N);
end
x = linspace(0, 1, N);
%Пробуем одну итерацию
% [psi, E] = reshenie('gamma');
% psi1_gamma = psi(1, :);
% psi2_gamma = psi(2, :);
% n1_gamma = psi1_gamma.^2*4.4*10^15;
% n2_gamma = psi2_gamma.^2*4*10^15;
% n = nd - n1_gamma - n2_gamma;
% n = symm_array(n);
% U_elstatic = puasson(n);
% Вставим добавочку для Х электронов
% koeff = 0.003;
% num_of_points = 50;
% U_add_max = koeff*num_of_points;
% U_add = ones(1, N)*U_add_max;
% for i = 1:num_of_points
%     U_add(483 - num_of_points + i) = U_add_max - koeff*i;
% end;
% for i = 1:num_of_points
%     U_add(483 + i) = U_add(483) + koeff*i;
% end;
% U_add = symm_array(U_add);
[psi, E] = reshenie('X');
psi_x = psi(1, :);




% n1 = psi1.^2*8.4*10^15;
% n = nd - n1;
% n = symm_array(n);
% U_elstatic = puasson(n);
% for i=1:5
%     %Решаем Шредингера для гамма электронов
%     [psi, E] = reshenie();
%     E_analize_gamma(2*i-1) = E(1);
%     E_analize_gamma(2*i) = E(2);    
%     psi1_gamma = psi(1, :);
%     psi2_gamma = psi(2, :);
%     %Решаем Шредингера для Х электронов
%     delta_U = -0.3;
%     [psi, E] = reshenie();
%     E_analize_X(i) = E(1);   
%     psi1_X = psi(1, :);    
%     %Задаем распределение концентрации в подзонах
%     n1_gamma = psi1_gamma.^2*4.4*10^15;
%     n2_gamma = psi2_gamma.^2*4*10^15;
%     n1_X = psi1_X.2*0.9;
%     %Полное распределение концентрации
%     n = nd - n1;
%     n = symm_array(n);
%     U_elstatic = puasson(n);
% end
plot(x, psi);



