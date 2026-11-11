
%% ACTIVIDAD FORMATIVA 2
% Simulacion y analisis de senales con la Transformada de Fourier
% Alumno: Jeronimo Bautista Torres
% Materia: Senales y Sistemas
% Profesor: Ing. Luis Osvaldo Moreno Gaytan
%
% Este programa genera tres senales, calcula su FFT, grafica magnitud
% y fase, y comprueba linealidad, desplazamiento temporal y escalamiento.

clear; close all; clc;

%% 1. Parametros de muestreo
Fs = 1000;                 % Frecuencia de muestreo (Hz)
T  = 1;                    % Duracion de observacion (s)
dt = 1/Fs;
t  = -T/2:dt:T/2-dt;       % Vector de tiempo centrado
N  = length(t);
f  = (-N/2:N/2-1)*(Fs/N);  % Vector de frecuencias centrado

%% 2. Definicion de senales elementales
% Pulso rectangular de amplitud 1 y duracion aproximada de 0.10 s
x_rect = double(abs(t) <= 0.05);

% Funcion escalon unitario: se observa dentro de una ventana finita
x_step = double(t >= 0);

% Senal senoidal de 50 Hz
f0 = 50;
x_sin = sin(2*pi*f0*t);

%% 3. Calculo de la Transformada Discreta de Fourier (FFT)
% Se normaliza por N y se centra el espectro con fftshift.
X_rect = fftshift(fft(x_rect))/N;
X_step = fftshift(fft(x_step))/N;
X_sin  = fftshift(fft(x_sin))/N;

%% 4. Graficas de las senales en el dominio del tiempo
fig1 = figure('Name','Senales en el tiempo','Color','w');
subplot(3,1,1);
plot(t,x_rect,'LineWidth',1.4); grid on;
title('Pulso rectangular'); xlabel('Tiempo (s)'); ylabel('Amplitud');
xlim([-0.25 0.25]);

subplot(3,1,2);
plot(t,x_step,'LineWidth',1.4); grid on;
title('Funcion escalon unitario'); xlabel('Tiempo (s)'); ylabel('Amplitud');
ylim([-0.1 1.1]);

subplot(3,1,3);
plot(t,x_sin,'LineWidth',1.2); grid on;
title('Senal senoidal de 50 Hz'); xlabel('Tiempo (s)'); ylabel('Amplitud');
xlim([0 0.1]);
exportgraphics(fig1,'01_senales_tiempo.png','Resolution',200);

%% 5. Magnitud y fase del espectro
fig2 = figure('Name','Espectros de Fourier','Color','w');
signals = {X_rect,X_step,X_sin};
names = {'Pulso rectangular','Funcion escalon','Senoidal de 50 Hz'};
for k = 1:3
    X = signals{k};
    subplot(3,2,2*k-1);
    plot(f,abs(X),'LineWidth',1.2); grid on;
    title(['Magnitud: ',names{k}]);
    xlabel('Frecuencia (Hz)'); ylabel('|X(f)|');
    xlim([-200 200]);

    subplot(3,2,2*k);
    phase = angle(X);
    phase(abs(X) < 1e-4) = NaN; % Evita interpretar fase donde la magnitud es casi cero
    plot(f,phase,'LineWidth',1.0); grid on;
    title(['Fase: ',names{k}]);
    xlabel('Frecuencia (Hz)'); ylabel('Fase (rad)');
    xlim([-200 200]);
end
exportgraphics(fig2,'02_espectros_magnitud_fase.png','Resolution',200);

%% 6. Propiedad de linealidad
% F{a*x1+b*x2} = a*F{x1}+b*F{x2}
a = 2;
b = 0.5;
x1 = x_rect;
x2 = x_sin;
x_comb = a*x1 + b*x2;
X_comb_direct = fftshift(fft(x_comb))/N;
X_comb_property = a*(fftshift(fft(x1))/N) + b*(fftshift(fft(x2))/N);
error_linealidad = max(abs(X_comb_direct-X_comb_property));

fig3 = figure('Name','Propiedad de linealidad','Color','w');
subplot(2,1,1);
plot(t,x_comb,'LineWidth',1.2); grid on;
title('Combinacion: 2 x pulso + 0.5 x senoidal');
xlabel('Tiempo (s)'); ylabel('Amplitud');
subplot(2,1,2);
plot(f,abs(X_comb_direct),'b','LineWidth',1.3); hold on;
plot(f,abs(X_comb_property),'r--','LineWidth',1.1); grid on;
xlim([-200 200]);
title('Comparacion de espectros: calculo directo y propiedad');
xlabel('Frecuencia (Hz)'); ylabel('Magnitud');
legend('FFT directa','Por linealidad');
exportgraphics(fig3,'03_linealidad.png','Resolution',200);

%% 7. Propiedad de desplazamiento en el tiempo
% Un desplazamiento temporal cambia la fase del espectro, no su magnitud.
tau = 0.08;                         % Desplazamiento de 80 ms
muestras = round(tau*Fs);
x_shift = circshift(x_rect,[0 muestras]);
X_shift = fftshift(fft(x_shift))/N;

fig4 = figure('Name','Desplazamiento temporal','Color','w');
subplot(2,1,1);
plot(t,x_rect,'LineWidth',1.2); hold on;
plot(t,x_shift,'--','LineWidth',1.2); grid on;
title('Pulso original y pulso desplazado');
xlabel('Tiempo (s)'); ylabel('Amplitud');
legend('Original','Desplazado');
subplot(2,1,2);
plot(f,abs(X_rect),'LineWidth',1.3); hold on;
plot(f,abs(X_shift),'--','LineWidth',1.2); grid on;
xlim([-200 200]);
title('Magnitudes espectrales: permanecen iguales');
xlabel('Frecuencia (Hz)'); ylabel('Magnitud');
legend('Original','Desplazado');
exportgraphics(fig4,'04_desplazamiento.png','Resolution',200);
error_desplazamiento = max(abs(abs(X_rect)-abs(X_shift)));

%% 8. Propiedad de escalamiento en el tiempo/frecuencia
% x(2t) es mas angosta en tiempo y su espectro se extiende en frecuencia.
x_compressed = double(abs(2*t) <= 0.05);
X_compressed = fftshift(fft(x_compressed))/N;

fig5 = figure('Name','Escalamiento','Color','w');
subplot(2,1,1);
plot(t,x_rect,'LineWidth',1.2); hold on;
plot(t,x_compressed,'--','LineWidth',1.2); grid on;
title('Escalamiento temporal: x(t) y x(2t)');
xlabel('Tiempo (s)'); ylabel('Amplitud');
legend('x(t)','x(2t)');
subplot(2,1,2);
plot(f,abs(X_rect),'LineWidth',1.2); hold on;
plot(f,abs(X_compressed),'--','LineWidth',1.2); grid on;
xlim([-300 300]);
title('Al comprimirse en tiempo, el espectro se ensancha');
xlabel('Frecuencia (Hz)'); ylabel('Magnitud');
legend('X(f)','Espectro de x(2t)');
exportgraphics(fig5,'05_escalamiento.png','Resolution',200);

%% 9. Resumen numerico de las verificaciones
fprintf('\\n===== RESULTADOS DE LA ACTIVIDAD =====\\n');
fprintf('Error maximo de linealidad: %.3e\\n',error_linealidad);
fprintf('Diferencia max. de magnitud por desplazamiento: %.3e\\n',error_desplazamiento);
fprintf('Frecuencia de la senal senoidal: %.1f Hz\\n',f0);
fprintf('Se generaron cinco archivos PNG con las graficas.\\n');
fprintf('Nota: el escalon se analiza en una ventana finita de observacion.\\n');
set(groot,'defaultFigureRenderer','painters')
