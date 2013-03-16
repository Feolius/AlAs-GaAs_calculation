%������ ��������� ������
clearvars
global j
%lattice_U = [];
%���������� ����� ��� ������������� ������ ��������� (���� ���� ���� ������
% 10)
N_apr = 10;
%���������� ������ �������� ������� ������ ���� ������ ����
global N_period
N_period = 40;
%������� ���� GaAs � ��
l_GaAs = 2.3;
%���������� ����� ��� ������������ ���� GaAs
global N_GaAs
N_GaAs = l_GaAs * N_apr;
%������� ���� AlAs � ��
l_AlAs = 1.1;
%���������� ����� ��� ������������ ���� AlAs
global N_AlAs
N_AlAs = l_AlAs * N_apr;
%������� ��� � ��
l_well = 36;
global N_well
N_well = l_well*N_apr;
%����� ������� �������� � ��
l_surface = 10;
%���������� ����� ��� ������������� ������� ��������
global N_surface
N_surface = l_surface*N_apr;
%����� ���� ������� � ��
global l_lattice
l_lattice = 2*(N_period/2*(l_GaAs + l_AlAs) + l_well/2 + l_surface);
%����� ����� �� �������
global N
N = l_lattice*N_apr;
%��������� ��������� ��� ����� ������� ��� ���� ������������
%������ ���� ������������ ����� GaAs � AlAs � ��
global delta_U
delta_U = 1;
%����� ������ � ��
global U_left
U_left = 3;
%������� ��������� ������ ���� �������. �� ����� ������������� � 12 �������
%� ���� GaAs
%������ ����
global N_start_delta
N_start_delta = 11*(l_GaAs + l_AlAs)*N_apr + l_surface*N_apr + 10;
global U_elstatic
%U_elstatic = zeros(1, N);
nd = zeros(1, N);
for i = N_start_delta + 1:N_start_delta + 3
    nd(i) = 4.65*10^15/(3*1/N);
end
x = linspace(0, 1, N);
%������� ���� ��������
% [psi, E] = reshenie('gamma');
% psi1_gamma = psi(1, :);
% psi2_gamma = psi(2, :);
% n1_gamma = psi1_gamma.^2*4.4*10^15;
% n2_gamma = psi2_gamma.^2*4*10^15;
% n = nd - n1_gamma - n2_gamma;
% n = symm_array(n);
% U_elstatic = puasson(n);
% ������� ��������� ��� � ����������
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
%     %������ ���������� ��� ����� ����������
%     [psi, E] = reshenie();
%     E_analize_gamma(2*i-1) = E(1);
%     E_analize_gamma(2*i) = E(2);    
%     psi1_gamma = psi(1, :);
%     psi2_gamma = psi(2, :);
%     %������ ���������� ��� � ����������
%     delta_U = -0.3;
%     [psi, E] = reshenie();
%     E_analize_X(i) = E(1);   
%     psi1_X = psi(1, :);    
%     %������ ������������� ������������ � ��������
%     n1_gamma = psi1_gamma.^2*4.4*10^15;
%     n2_gamma = psi2_gamma.^2*4*10^15;
%     n1_X = psi1_X.2*0.9;
%     %������ ������������� ������������
%     n = nd - n1;
%     n = symm_array(n);
%     U_elstatic = puasson(n);
% end
plot(x, psi);



