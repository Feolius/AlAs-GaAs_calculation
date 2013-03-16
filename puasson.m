%Решает уравнение пуассона методом рунге кутта
function fi=puasson(n)
global l_lattice
global N
%Заряд электрона
e = 1.6 * 10^-19;
%Электрическая постоянная
eps_0 = 8.85 * 10^-12;
%Диэлектрическая проницаемость
eps = 12;
%коэффициент в правой части уравнения
koeff = e/(eps*eps_0)*(l_lattice*10^-9)^2;
%Шаг сетки
h = 1/N;
%Решаем сначала до дельта слоя
%Массив с удвоенным колличеством элементов для рассчета
n_practical = zeros(1, 2*N - 1);
for i=1:N-1
    n_practical(2*i - 1) = n(i);
    n_practical(2*i) = (n(i+1) + n(i))/2;
end
n_practical(2*N - 1) = n(N);
k = zeros(1, 4);
q = zeros(1, 4);
fi1 = zeros(1, N);
fi2 = zeros(1, N);
for i=1:N-1
    k(1) = fi2(i);
    q(1) = koeff*n_practical(2*i - 1);
    k(2) = fi2(i) + 0.5*h*q(1);
    q(2) = koeff*n_practical(2*i);
    k(3) = fi2(i) + 0.5*h*q(2);
    q(3) = koeff*n_practical(2*i);
    k(4) = fi2(i) + h*q(3);
    q(4) = koeff*n_practical(2*i + 1);
    fi1(i + 1) = fi1(i) + 1/6*h*(k(1) + 2*k(2) + 2*k(3) + k(4));
    fi2(i + 1) = fi2(i) + 1/6*h*(q(1) + 2*q(2) + 2*q(3) + q(4));
end
fi1 = fi1/(l_lattice*10^-9);
fi = symm_array(fi1);
