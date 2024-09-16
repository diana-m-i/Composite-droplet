clc;
close all;
clear all;
%n-Dodecane, Н-Додекан C12H26
%Начальная температура углеводорода и окружающего каплю воздуха, К
Theta_Oil_0 = 343; 
Theta_Gas = 700;  
%Скорость обтекающего потока воздуха, м/c
U_Gas = 2;
%Безразмерная температура кипения воды, К
Water_Boil = (373 - Theta_Gas)/Theta_Gas;

%Теплофизические свойства углеводорода
%Зависимость плотности углеводорода от температуры, кг/м^3
Func_Rho_Oil = @(T) 744.11 - 0.771*(T - 300); %N-dodecane
Rho_Oil_dim = Func_Rho_Oil(Theta_Oil_0); 

%Зависимость теплоемкости углеводорода от температуры, Дж/(кг*К)
Func_C_Oil = @(T)(2.18 + 0.0041*(T - 300))*1e3;  %N-dodecane

%Зависимость теплопроводности углеводорода от температуры, Вт/(м*К)
Func_Lambda_Oil = @(T) (0.1405 - 0.00022*(T - 300)); %N-dodecane

C_Oil_dim = Func_C_Oil(Theta_Oil_0);
Lambda_Oil_dim = Func_Lambda_Oil(Theta_Oil_0);
%Коэффициент температуропроводности, м^2/c
Kappa_Oil_dim = Lambda_Oil_dim/(Rho_Oil_dim*C_Oil_dim); 
Theta_Oil_Boil = (489 - Theta_Gas)/Theta_Gas;

%Параметры капли
%Диаметр всей капли, м 
D_Oil_dim = 200*1e-6; 
%Объем капли, м^3
V_Oil_dim = 4*pi*(D_Oil_dim/2)^3/3; 
%Диаметр капли воды, м
D_Water_dim = 0.15^(1/3)*D_Oil_dim; 
%Объемное содержание воды
Content_W = (D_Water_dim/D_Oil_dim)^3; 

%Число Пекле - При малых значениях преобладает молекулярная
%теплопроводность, а при больших — конвективный перенос теплоты.
Peclet = D_Oil_dim * U_Gas/Kappa_Oil_dim;

%Теплофизические свойства воды
%Плотность, кг/м^3
Rho_Water_dim = 1e3; 
%Теплоемкость, Дж/(К*кг)
C_Water_dim = 4.2*1e3;
%Теплопроводность, Вт/(К*м)
Lambda_Water_dim = 0.66; 
%Коэффициент температуропроводности, м^2/c;
Kappa_Water_dim = Lambda_Water_dim/(C_Water_dim*Rho_Water_dim);

%Теплофизические свойства воздуза при температуре 489 К
%Теплопроводность, Вт/(К*м)
Lambda_Gas_dim = 0.0548; 
%Плотность
Rho_Gas_dim = 0.49;
%Динамическая вязкость, Па*с
Mu_Gas_dim = 34.6*1e-6;
%Число Рейнольдса
Re = D_Oil_dim*U_Gas*Rho_Gas_dim/Mu_Gas_dim;
%Число Прандтля для воздуха
Pr = 0.7;

%Формула Ранца-Маршалла, число Нуссельта
Nu = 2 + 0.6 * Re^(1/2) * Pr^(1/3);

%Коэффициент теплоотдачи с поверхности капли, Вт/(м^2*К)
Alpha_Oil_dim = Lambda_Gas_dim*Nu/D_Oil_dim;
%Доля радиационного теплового потока от конвективного
J_conv = Alpha_Oil_dim * (Theta_Gas - Theta_Oil_0);
J_rad = 5.67 * 1e-8 * (Theta_Gas^4 - Theta_Oil_0^4); 
Relation = J_rad/J_conv; 

%Переход к безразмерному виду
Kappa_Water = Kappa_Water_dim/Kappa_Oil_dim;
Kappa_Oil = Kappa_Oil_dim/Kappa_Oil_dim;
R_Water = (0.5*D_Water_dim)/(0.5*D_Oil_dim);
Alpha_Oil = Alpha_Oil_dim*0.5*D_Oil_dim/Lambda_Oil_dim;
Lambda_Water = Lambda_Water_dim/Lambda_Oil_dim;
Lambda_Oil = 1;

%Вычисление собственных значений оператора Лапласа
a = 0; b = 120; n = 180;

w0 = [];
w = linspace(a, b, n);
P = zeros(1,length(w));
%Шаг по пространству
h = 0.0001;
r_water = 0:h:R_Water-h;
r_oil = R_Water+h:h:1;
%Пространственная сетка в воде и углеводороде
r = [r_water R_Water r_oil];
%Вычисляем значения характеристической функции Psi
for i=1:length(w)
    P(i) = Psi(w(i), Alpha_Oil, Kappa_Water, Lambda_Water, R_Water);
end

%Ищем нули функции Psi на отрезке [a,b]
iter = 1;
for i=2:n
     if (P(i-1)*P(i)<0)                       %((P(i) > 0 && P(i-1) < 0) || (P(i) < 0 && P(i-1) > 0))
        w0(iter) = fzero(@Psi, w(i), [], Alpha_Oil, Kappa_Water, Lambda_Water, R_Water);
        iter = iter + 1;
    end
end

%График характеристической функции и ее корней
plot(w,P,'b' ,'LineWidth', 1.5)
hold on
plot(w0,0,'Marker','o','MarkerEdgeColor','r','LineWidth', 1.5);
% hold on
yline(0,'k--', 'LineWidth', 1)
xlabel('\omega','FontSize', 18, 'FontName', 'Times New Roman')
ylabel('\Psi(\omega)','FontSize', 18, 'FontName', 'Times New Roman')
set(gca,'FontSize', 14, 'FontName', 'Times New Roman')
xticks([0 20 40 60 80 100])
%'Rotation',0

%Собственные функции оператора Лапласа
eigU = zeros(length(w0),length(r)); 
C1 = zeros(1,length(w0));
C2 = C1;
for i=1:length(w0)
    [C1(i), C2(i)] = reqcoeff(w0(i), R_Water, Kappa_Water, Lambda_Water);
end

for i=1:length(w0)
    for j = 1:length(r)
        eigU(i,j) = eigenfunction(w0(i), r(j), C1(i), C2(i), R_Water, Kappa_Water); 
    end
end

%График собственных функций Xn, n = 3, 5, 14
figure;
plot(r, eigU(3,:), 'LineWidth', 1.5)
hold on
plot(r, eigU(5,:), 'LineWidth', 1.5)
plot(r, eigU(14,:), 'LineWidth', 1.5)
% plot([R_Water, R_Water],[-1, 1],'k--', 'LineWidth', 1.5)
xline(R_Water,'k--', 'LineWidth', 1.5);
plot([0, 1],[0, 0],'k--', 'LineWidth', 1.5)
xlabel('r^*','FontSize', 18, 'FontName', 'Times New Roman')
ylabel('X_n(r^*)','FontSize', 18, 'FontName', 'Times New Roman')
ylim([min(eigU,[],'all'), max(eigU,[],'all')]);
set(gca,'FontSize', 14, 'FontName', 'Times New Roman')
legend(['n = ' num2str(3)],['n = ' num2str(5)],...
    ['n = ' num2str(14)], 'FontSize',18, 'FontName', 'Times New Roman')

%Норма собственных функций Xn
%Веса при интегрировании различны для воды и углеводорода
weight1 = Lambda_Water/Kappa_Water;
weight2 = Lambda_Oil/Kappa_Oil;
function_U_norm = @(r,w,C1,C2,c, Rw, kappa_w) c*(eigenfunction(w, r, C1, C2, Rw, kappa_w).*r).^2;
eigUnorm = zeros(1,length(w0));

for i=1:length(w0)
     eigUnorm(i) = integral(@(r) function_U_norm(r, w0(i), C1(i), C2(i), weight1, R_Water, Kappa_Water),0, R_Water)+...
         integral(@(r) function_U_norm(r, w0(i), C1(i), C2(i), weight2, R_Water, Kappa_Water),R_Water+h, 1);
end

%График нормы собственных функций
figure;
bar(sqrt(eigUnorm)*10^3, 'b')
xlabel('n','FontSize', 18,'FontName', 'Times New Roman')
ylabel('||X_n(r^*)||*10^3','FontSize', 18,'FontName', 'Times New Roman')
set(gca,'FontSize', 14, 'FontName', 'Times New Roman')
%Коэффициенты разложения A_n
An = zeros(length(w0),1);
%Равномерное начальное распределение температуры (в общем случае может
%зависеть от координаты r)
theta_0 = @(r) (Theta_Oil_0 - Theta_Gas)/Theta_Gas;

function_An = @(r, w, C1, C2, c, Rw, kappa_w) c*...
    eigenfunction(w, r, C1, C2, Rw, kappa_w).*theta_0(r).*r.^2;

for i=1:length(w0)
    %[C1, C2] = reqcoeff(w0(i), R_Water, Kappa_Water, Lambda_Water);
     An(i) = (integral(@(r) function_An(r, w0(i), C1(i), C2(i),...
         weight1, R_Water, Kappa_Water),0,R_Water)+...
         integral(@(r) function_An(r, w0(i), C1(i), C2(i),...
         weight2, R_Water, Kappa_Water),R_Water+h, 1))/eigUnorm(i);
end
%График 4
figure;
bar(An, 'b')
xlabel('n','FontSize', 18, 'FontName', 'Times New Roman')
ylabel('A_n','FontSize', 18, 'FontName', 'Times New Roman')
set(gca,'FontSize', 14, 'FontName', 'Times New Roman')

%Общее решение задачи 
%Сетка по времени
t = 0:0.0005:1;
theta_boil = 0;
left = Water_Boil - 0.0002;
right = Water_Boil + 0.0002;

for j = 1:length(t)   
        for n=1:length(An)  
            eig_R_Water = eigenfunction(w0(n),R_Water, C1(n), C2(n), R_Water, Kappa_Water);
            theta_boil = theta_boil + An(n)*eig_R_Water*exp(-w0(n)^2*t(j));
        end
        if (theta_boil >= left) && (theta_boil <= right)
            %Время начала кипения воды, исходя из аналитического решения
            t_boil = t(j)*((0.5*D_Oil_dim)^2)/Kappa_Oil_dim;
            tb1 = t(j);
            t_boil_j1 = j;
            break;
        end
        theta_boil = 0;
end

theta1 = zeros(length(t),length(r));
for i = 1:length(r)
    for j = 1:length(t)   
        for n=1:length(An)  
            theta1(j,i) = theta1(j,i) + An(n)*eigU(n,i)*exp(-w0(n)^2*t(j));
        end
    end
end

%Графики динамики прогрева, аналитическое решение
r1 = 0:0.03:1;
figure
plot(r, theta1(1,:),'LineWidth', 2, 'Color', [0.6350 0.0780 0.1840])
hold on 
plot(r, theta1(200,:),'LineWidth', 2, 'Color', [0.4660 0.6740 0.1880])
plot(r, theta1(300,:),'LineWidth', 2,'Color', [0.4940 0.1840 0.5560])
plot(r, theta1(t_boil_j1,:), 'LineWidth', 2,'Color', [0.8500 0.3250 0.0980])
plot(r1, theta_0(r1),'r*','MarkerSize',5)
xline(R_Water,'k--', 'LineWidth', 1.5) 
yline(Water_Boil,'k--', 'LineWidth', 1.5)
xlabel('r^*','FontSize', 18, 'FontName', 'Times New Roman')
ylabel('\theta^*(r^*,t^*)','FontSize', 18, 'FontName', 'Times New Roman')
set(gca,'FontSize', 14, 'FontName', 'Times New Roman')
legend(['t^* = ' num2str(t(1))],['t^* = ' num2str(t(200))], ...
    ['t^* = ' num2str(t(300))],['t^* = ' num2str(t(t_boil_j1))],...
    'FontSize',18, 'FontName', 'Times New Roman')

%Проверка температуры на границе капли
if theta1(t_boil_j1, end) > Theta_Oil_Boil
    print('Кипение углеводорода, некорректное решение')
end

%Численное решение
%Дискретизируем отрезки [0,R_Water],[R_Water, R_Oil]
h = 0.001;
r_water = 0:h:R_Water-h;
r_oil = R_Water+h:h:1;
r2 = [r_water R_Water r_oil];
%Количество узлов сетки в воде
N_W = length(r_water) + 1; 
%Общее количество узлов
N_O = N_W + length(r_oil);
%Коэффициент излучения
eps0 = 0.5;
sigma_sb = 5.67*1e-8;
%Критерий, необходимый для учета теплообмена за счет излучения на границе
%капли. Если Q_rad = 0, на границе происходит только конвективный теплообмен
% Q_rad = eps0*sigma_sb*Theta_Gas^3*(D_Oil_dim/2)/Lambda_Oil_dim;
Q_rad = 0;
theta20 = zeros(1,N_O);
theta20(:) = (Theta_Oil_0 - Theta_Gas)/Theta_Gas;
options=odeset('Events', @(t2,theta2) BoilingEvents(t2,theta2,Water_Boil,N_W));

[t2, theta2,te,ye,ie] = ode15s(@(t2,theta2) temp_eval(t2, theta2,...
    Kappa_Water,Lambda_Water,Alpha_Oil,r2,N_W,N_O,h,Q_rad),...
    [0 1], theta20, options);

%Время начала кипения воды, полученное из численного решения
t_boil_numerical = te*(0.5*D_Oil_dim)^2/Kappa_Oil_dim;
%Графики динамики прогрева
%Оценка корректности проводится из сопоставления аналитического решения theta_1 в
%момент начала кипения воды t_boil_j1 и численного ye(t=te)
figure
plot(r, theta1(1,:),'LineWidth', 1.5)
hold on 
plot(r, theta1(200,:),'LineWidth', 1.5)
plot(r, theta1(t_boil_j1,:), 'LineWidth', 1.5)
plot(r2, ye,'.','MarkerSize',25, 'MarkerEdgeColor','b','MarkerIndices',1:50:length(r2))
xline(R_Water,'k--', 'LineWidth', 1.5) 
yline(Water_Boil,'k--', 'LineWidth', 1.5)
xlabel('r^*','FontSize', 18, 'FontName', 'Times New Roman')
ylabel('\theta^*(r^*,t^*)','FontSize', 18, 'FontName', 'Times New Roman')
set(gca,'FontSize', 14, 'FontName', 'Times New Roman')

disp(['Время начала кипения воды, полученное из аналитического решения, с: ', num2str(t_boil)])
disp(['Время начала кипения воды, полученное из численного решения, с: ', num2str(t_boil_numerical)])
%Остановить интегрирование при достижении температуры кипения воды
%Проверяем точку на границе раздела сред
function [value,isterminal,direction] = BoilingEvents(t, y, Theta_WB, N_W)
  value = y(N_W) - Theta_WB; % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end

function dydt = temp_eval(t,y,Kappa_Water,Lambda_Water,Alpha_Oil,r,N_W,N_O,h,Q_rad)
dydt = zeros(N_O,1);
%ОДУ для температуры в центре капли
dydt(1) = 2*Kappa_Water*(y(2)-y(1))/h^2;
%ОДУ для температуры в воде
for i = 2:N_W-1
    dydt(i) = Kappa_Water*((y(i+1)-y(i))*r(i)/(1-r(i)/r(i+1)) - ...
        (y(i)-y(i-1))*r(i-1)/(1-r(i-1)/r(i)))/(h*r(i)^2); 
end
%ОДУ для температуры в углеводороде
for i = N_W+1:N_O-1
    dydt(i) = ((y(i+1)-y(i))*r(i)/(1-r(i)/r(i+1)) - ...
        (y(i)-y(i-1))*r(i-1)/(1-r(i-1)/r(i)))/(h*r(i)^2);    
end
%ОДУ для температуры на границе раздела сред
dydt(N_W) = (r(N_W)*(y(N_W+1) - y(N_W))/(1-r(N_W)/r(N_W+1)) -...
    Lambda_Water*r(N_W-1)*(y(N_W) - y(N_W-1))/(1-r(N_W-1)/r(N_W)))/(h*r(N_W)^2);
%Тепловое излучение
dydt(N_O) = -2*(Alpha_Oil*y(N_O) + Q_rad*((y(N_O)+1)^4 - 1))/h - ...
     2*r(N_O-1)*(y(N_O) - y(N_O-1))/(r(N_O)^2*h*(1-r(N_O-1)/r(N_O)));

end

%Сферические функции Бесселя
function [sphbes] = sbessel(order, w, r, kappa, flag)
x = w*r/sqrt(kappa);
if flag == 1
    %flag = 1, сферические функции Бесселя 1-го рода
    sphbes = sqrt(pi./(2*x)).*besselj(order+0.5, x);
elseif flag == 2
    %flag = 2, сферические функции Бесселя 2-го рода(Неймана)
    sphbes = sqrt(pi./(2*x)).*bessely(order+0.5, x); %flag = 2, Neyman 
end
end

%Вычисление значений характеристической функции
function [P] = Psi(w, alpha, kappa_w, lambda_w, Rw)
[C1,C2] =  reqcoeff(w, Rw, kappa_w, lambda_w);
kappa_o = 1;
P = -w*(C1*sbessel(1,w,1,kappa_o,1)+C2*sbessel(1,w,1,kappa_o,2)) +...
    alpha*(C1*sbessel(0,w,1,kappa_o,1)+C2*sbessel(0,w,1,kappa_o,2));
end

%Функция вычисления коэффициентов C1, C2
function [C1, C2] = reqcoeff(w, Rw, kappa_w, lambda_w)
 kappa_o = 1;
 Aw = 1;
 nu1 =lambda_w*Aw*(sbessel(0,w,Rw,kappa_o,1) - sbessel(0,w,Rw,kappa_o,2))*...
     sbessel(1,w,Rw,kappa_w,1)/sqrt(kappa_w) -...
 Aw*(sbessel(1,w,Rw,kappa_o,1) - sbessel(1,w,Rw,kappa_o,2))*sbessel(0,w,Rw,kappa_w,1);

 nu2 = (sbessel(0,w,Rw,kappa_o,1)*sbessel(1,w,Rw,kappa_o,2) - ...
     sbessel(0,w,Rw,kappa_o,2)*sbessel(1,w,Rw,kappa_o,1));
 
 nu = nu1/nu2;
 C1 = (Aw*sbessel(0,w,Rw,kappa_w,1) - nu*sbessel(0,w,Rw,kappa_o,2))/...
   (sbessel(0,w,Rw,kappa_o,1) - sbessel(0,w,Rw,kappa_o,2)); 
 C2 = (nu*sbessel(0,w,Rw,kappa_o,1) - Aw*sbessel(0,w,Rw,kappa_w,1))/...
   (sbessel(0,w,Rw,kappa_o,1) - sbessel(0,w,Rw,kappa_o,2)) ;
end

%Собственные функции в объеме капли
function [eigenU] = eigenfunction(w, r, C1, C2, Rw, kappa_w)
Aw = 1;
kappa_o = 1;
if r <= Rw
    eigenU = Aw*sbessel(0, w, r, kappa_w, 1);
else
    eigenU = C1*sbessel(0,w,r,kappa_o,1) + C2*sbessel(0,w,r,kappa_o,2);
end
end