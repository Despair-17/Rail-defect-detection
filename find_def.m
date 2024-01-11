clear all;
close all;
%% загрузка данных акселерометров, временных меток и показаний одометра
% load('ADXL354_X.mat')
% load('ADXL1001_Z.mat')
load('ADXL1002_Z.mat')
% load('ADXL1003_Z.mat')
load('timeStamps.mat')
load('odo.mat')
%% Создание отдельных массивов показаний акселерометров для каждого колеса
% ADXL354_1=ADXL354_X{1};
% ADXL354_2=ADXL354_X{2};
% ADXL354_3=ADXL354_X{3};
% ADXL354_4=ADXL354_X{4};
% ADXL1001_1=ADXL1001_Z{1};
% ADXL1001_2=ADXL1001_Z{2};
% ADXL1001_3=ADXL1001_Z{3};
% ADXL1001_4=ADXL1001_Z{4};
ADXL1002_1=ADXL1002_Z{1};
ADXL1002_2=ADXL1002_Z{2};
ADXL1002_3=ADXL1002_Z{3};
ADXL1002_4=ADXL1002_Z{4};
% ADXL1003_1=ADXL1003_Z{1};
% ADXL1003_2=ADXL1003_Z{2};
% ADXL1003_3=ADXL1003_Z{3};
% ADXL1003_4=ADXL1003_Z{4};
%% Создание отдельных массивов скорости, координат и временных меток для каждого колеса  
Vel_1=odo{1,2};
Vel_2=odo{2,2};
Vel_3=odo{3,2};
Vel_4=odo{4,2};
Coord_1=odo{1,1};
Coord_2=odo{2,1};
Coord_3=odo{3,1};
Coord_4=odo{4,1};
timeStamps_1=timeStamps{1};
timeStamps_2=timeStamps{2};
timeStamps_3=timeStamps{3};
timeStamps_4=timeStamps{4};
% График изменения скорости 
figure
plot(1:length(Vel_3),Vel_3,'-r',1:length(Vel_4),Vel_4,'-b')
%% График необработанных сигналов ускорений 3 и 4 колеса
figure
plot(1:length(ADXL1002_3),ADXL1002_3,'-r',1:length(ADXL1002_4),ADXL1002_4,'-b')
    title('Зависимость ускорений акселерометров от отсчетов')
    xlabel('Отсчеты')
    ylabel('Ускорения, м/с^2')
    legend('Впереди идущее','Следом идущее')
    grid on
    hold on
figure
plot(Coord_3,ADXL1002_3,'-r',Coord_4,ADXL1002_4,'-b')
    title('Зависимость ускорений акселерометров от координат')
    xlabel('Путевая координата, м')
    ylabel('Ускорение, м/с^2')
    legend('Переднее колесо','Заднее колесо')  
    grid on
figure
plot(timeStamps_3,ADXL1002_3,'-r',timeStamps_4,ADXL1002_4,'-b')
    title('Зависимость ускорений акселерометров от временных меток')
    xlabel('Временных меток, с')
    ylabel('Ускорение, м/с^2')
    grid on
%% Нормирование сигнала от скорости поезда
aan_1=(ADXL1002_1)./((Vel_1).^2);
aan_2=(ADXL1002_2)./((Vel_2).^2);
aan_3=(ADXL1002_3)./((Vel_3).^2);
aan_4=(ADXL1002_4)./((Vel_4).^2);
% Приведение ускорений к нулевой линии, вычитанием средних уровней
aan_1=aan_1-mean(aan_1);
aan_2=aan_2-mean(aan_2);
aan_3=aan_3-mean(aan_3);
aan_4=aan_4-mean(aan_4);
% График нормированных ускорений 
% figure
% plot(1:length(aan_3),aan_3,'-r',1:length(aan_4),aan_4,'-b')
    % title('Зависимость нормированных ускорений акселерометров от отсчетов')
    % xlabel('Отсчеты')
    % ylabel('Нормированные ускорения, м^-1')
    % legend('Переднее колесо','Заднее колесо')
    % grid on
    % hold on
figure
plot(Coord_3,aan_3,'-r',Coord_4,aan_4,'-b')
    title('Зависимость нормированных ускорений акселерометров от координат')
    xlabel('Путевая координата, м')
    ylabel('Нормированные ускорения, м^-1')
    legend('Переднее колесо','Заднее колесо') 
grid on
% figure
% plot(timeStamps_3,aan_3,'-r',timeStamps_4,aan_4,'-b')
    % title('Зависимость нормированных ускорений акселерометров от временных меток')
    % xlabel('Временные метки, с')
    % ylabel('Нормированные ускорения, м^-1')
    % grid on
%% Автоматический поиск дефектов
% Расстояние между колесными парами
S=2.5;
% Максимальная длинна импульсного дефекта
Lmax=0.25;
% Количество отсчетов
No=length(aan_3);
% Частота съема данных
f=No/abs((timeStamps_4(end)-timeStamps_3(1)));
% Разница во времени начала съема данных 
Delay_S1=f*abs((timeStamps_1(1)-timeStamps_2(1)));
Delay_S2=f*abs((timeStamps_3(1)-timeStamps_4(1)));
% Абсолютные величины ускорений 
aan_1F=abs(aan_1);
aan_3F=abs(aan_3);
% Огибающая показаний 
[aan_1_up1,~]=envelope(aan_1F,200,'peak');
aan_1_up1(aan_1_up1<0)=0;
[aan_3_up1,~]=envelope(aan_3F,200,'peak');
aan_3_up1(aan_3_up1<0)=0;
% Минмальный порог наличия дефекта
a1=zeros(length(aan_1_up1), 1);
a1(:,1)=1.5;
a2=zeros(length(aan_3_up1), 1);
a2(:,1)=1.5;
figure
plot(1:length(aan_3_up1), aan_3_up1, 1:length(a2), a2, 'r')
    xlabel('Отсчеты')
    ylabel('Нормированные ускорения, м^-1')
    legend('Огибающая модуля сигнала','Пороговая линия') 
    grid on
    hold on
% Цикл автоматического поиска точкек предплогаемого начала дефекта, 
% по показаниям акселерометра на первом колесе 
i=500;
m=0;
while i<length(aan_1_up1)-22000
    if aan_1_up1(i)<1.5
        i=i+1;
        continue
    else
        m=m+1;
        n1(m,1)=i-100;
        i=i+3500;
    end
end
clear i m aan_1_up1 aan_1F

i=500;
m=0;
while i<length(aan_3_up1)-22000
    if aan_3_up1(i)<1.5
        i=i+1;
        continue
    else
        m=m+1;
        n2(m,1)=i-100;
        i=i+3500;
    end
end
clear i m aan_3_up1 aan_3F
% Расчет точек показаний аксселерометра на следом идущем колесе
m=0;
for i=1:1:length(n1)
    Delay_P1G=S*f/Vel_2(n1(i));
    Delay_P1=S*f/mean(Vel_2(n1(i):n1(i)+round(Delay_P1G)));
    clear Delay_P2G
    m=m+1;
    n1_1(m,1)=n1(i);
    n2_2(m,1)=round(n1(i)+Delay_P1+Delay_S1);
end
clear i m Delay_P1 n

m=0;
for i=1:1:length(n2)
    Delay_P2G=S*f/Vel_4(n2(i));
    Delay_P2=S*f/mean(Vel_4(n2(i):n2(i)+round(Delay_P2G)));
    clear Delay_P2G
    m=m+1;
    n3_3(m,1)=n2(i);
    n4_4(m,1)=round(n2(i)+Delay_P2+Delay_S2);
end
clear i m Delay_P2 n
%% Цикл проверки предпологаемого дефекта
% Первый рельс
result1=[0,0,0];
c=0;
for a=1:1:length(n1_1)
    n1=n1_1(a);
    n2=n2_2(a);
    % Расчет длинны дефекта
    long_1D=round(Lmax*f/Vel_1(n1));
    % Задержка для второго колеса 
    Delay_P1G=S*f/Vel_1(n1);
    Delay_t1=S/mean(Vel_1(n1:n1+round(Delay_P1G)));
    % Ускорения акселерометров на участке с дефектом
    aanD_1=aan_1(n1:n1+long_1D);
    aanD_2=aan_2(n2:n2+long_1D);
    figure
    plot(1:length(aanD_1),aanD_1, '-r', 1:length(aanD_2),aanD_2, '-b')
        title('Нормированные ускорения от отсчетов на участке с дефектом')
        xlabel('Отсчеты с частотой f, Гц')
        ylabel('Нормированные ускорения, м^-1')
        legend('Переднее колесо','Заднее колесо') 
        grid on
    % Синхронизация показаний акселерометров при присутствии корреляции
    i=0;
    s=400;
    b=s*0.75;
    while i<s
        y1=circshift(aanD_1,-round(b)+i);
        y2=aanD_2(1:end);
        Ri=corr(y1,y2);
        R1(i+1,:)=Ri;
        clear Ri y1 y2
        i=i+1;
    end
    clear i
    % Расчет сдвижки сигнала по положению максимума рассчитанной корреляционной функции 
    [g,k1]=max(R1);
    u1=(k1)-round(b);
    if g<0.3
        figure
        plot(-s/2:s/2-1 ,R1,'-b')
            title('Коэффициент корреляции от числа сдвижек')
            xlabel('Сдвижки')
            ylabel('Коэффициент корреляции')
            grid on
    end
    clear k1
    % Участок с дефектом с учетом вычисленной сдвижки сигналов
    y1=aanD_1(1:end);
    y2=aan_2(n2+u1:n2+long_1D+u1);
    clear u1 
    if g<0.3
        figure
        plot(1:length(y1), y1 ,'-r', 1:length(y2), y2, '-b')
            title(['Совмещенные нормированные ускорения от отсчетов ' ...
                'на участке с дефектом'])
            xlabel('Отсчеты с частотой f, Гц')
            ylabel('Ускорение, м/с^-1')
            legend('Переднее колесо','Заднее колесо') 
            grid on
    end
    % Определение длинны дефекта и глубинны дефекта
    if g>0.3
        [Rr1,ff1]=xcorr(y1,y2, 1*s, 'normalized');
        figure
        plot(ff1, Rr1)
            title('Коэффициент корреляции от числа сдвижек')
            xlabel('Сдвижки')
            ylabel('Коэффициент корреляции')
            grid on
        [~,M1]=max(Rr1);
        [~,M_Left1]=min(Rr1(1:M1-1));
        [~,M_Right1]=min(Rr1(M1+1:end));
        M_Right1=M1+M_Right1;
        % Длинна дефекта (расстояние между минимумами корреляционной функции)
        Long_def_f1=M_Right1-M_Left1;
        Long_def_m1=(Vel_1(n1)*Long_def_f1)/f;
        clear M_Right1 M_Left1 M1 Rr1 ff1 Long_def_1
        % Глубина дефекта
        A_max1=max(y1);
        A_max2=max(y2);
        if A_max1>=A_max2
            Amax_1=A_max1;
        else
            Amax_1=A_max2;
        end
        High_def_f1=(Amax_1*Long_def_m1^2)/(4*pi^2);
        clear A_max1 A_max2 Amax_1
        % Координаты найденного дефекта
        % Округленные значения высоты и длины дефекта 
        c=c+1;
        result1(c,1)=Coord_1(n1);
        result1(c,2)=Long_def_m1*1000;
        result1(c,3)=High_def_f1*1000;
        clear g Long_def_m1 High_def_f1
    end
end
clear long_1D Delay_S1 Delay_P1G y1 y2 R1 a n1 n2
delete('result1.txt');
save result1.txt result1 -ascii;
%% Второй рельс
result2=[0,0,0];
c=0;
for a=1:1:length(n3_3)
    n3=n3_3(a);
    n4=n4_4(a);
    % Расчет длинны дефекта
    long_2D=round(Lmax*f/Vel_3(n3));
    % Задержка для второго колеса 
    Delay_P2G=S*f/Vel_3(n3);
    Delay_t2=S/mean(Vel_3(n3:n3+round(Delay_P2G)));
    % Ускорения акселерометров на участке с дефектом
    aanD_3=aan_3(n3:n3+long_2D);
    aanD_4=aan_4(n4:n4+long_2D);
    figure
    plot(1:length(aanD_3),aanD_3, '-r', 1:length(aanD_4),aanD_4, '-b')
        title('Нормированные ускорения от отсчетов на участке с дефектом')
        xlabel('Отсчеты с частотой f, Гц')
        ylabel('Нормированные ускорения, м^-1')
        legend('Переднее колесо','Заднее колесо') 
        grid on
    % Синхронизация показаний акселерометров при присутствии корреляции
    i=0;
    s=400;
    b=s*0.75;
    while i<s
        y3=circshift(aanD_3,-round(b)+i);
        y4=aanD_4(1:end);
        Ri=corr(y3,y4);
        R2(i+1,:)=Ri;
        clear Ri y3 y4
        i=i+1;
    end
    clear i
    % Расчет сдвижки сигнала по положению максимума рассчитанной 
    % корреляционной функции 
    [g,k2]=max(R2);
    u2=(k2)-round(b);
    if g<0.3
        figure
        plot(-s/2:s/2-1 ,R2,'-b')
            title('Коэффициент корреляции от числа сдвижек')
            xlabel('Сдвижки')
            ylabel('Коэффициент корреляции')
            grid on
    end
    clear k1 k2 
    % Участок с дефектом с учетом вычисленной сдвижки сигналов
    y3=aanD_3(1:end);
    y4=aan_4(n4+u2:n4+long_2D+u2);
    clear u1 u2
    if g<0.3
        figure
        plot(1:length(y3), y3 ,'-r', 1:length(y4), y4, '-b')
            title(['Совмещенные нормированные ускорения от отсчетов ' ...
                'на участке с дефектом'])
            xlabel('Отсчеты с частотой f, Гц')
            ylabel('Ускорение, м/с^-1')
            legend('Переднее колесо','Заднее колесо') 
            grid on
    end
    % Определение длинны дефекта и глубинны дефекта
    if g>0.3
        [Rr2,ff2]=xcorr(y3,y4, 1*s, 'normalized');
        figure
        plot(ff2, Rr2)
            title('Коэффициент корреляции от числа сдвижек')
            xlabel('Сдвижки')
            ylabel('Коэффициент корреляции')
            grid on
        [~,M2]=max(Rr2);
        [~,M_Left2]=min(Rr2(1:M2-1));
        [~,M_Right2]=min(Rr2(M2+1:end));
        M_Right2=M2+M_Right2;
        % Длинна дефекта (расстояние между минимумами корреляционной функции)
        Long_def_f2=M_Right2-M_Left2;
        Long_def_m2=(Vel_3(n3)*Long_def_f2)/f;
        clear M_Right2 M_Left2 M2 Rr2 ff2 Long_def_f2
        % Глубина дефекта
        A_max3=max(y3);
        A_max4=max(y4);
        if A_max3>=A_max4
            Amax_2=A_max3;
        else
            Amax_2=A_max4;
        end
        High_def_f2=(Amax_2*Long_def_m2^2)/(4*pi^2);
        clear A_max3 A_max4 Amax_2
        % Координаты найденного дефекта
        % Округленные значения высоты и длины дефекта 
        c=c+1;
        result2(c,1)=Coord_3(n3);
        result2(c,2)=Long_def_m2*1000;
        result2(c,3)=High_def_f2*1000;
        clear g Long_def_m2 High_def_f2
    end
end
clear long_2D Delay_S2 Delay_P2G y3 y4 R2 a b n3 n4
delete('result2.txt');
save result2.txt result2 -ascii;