% Luan Cong Doan - Final Exam - Vibration - Question 3 - DFT
close all; clear all; clc;


%% fft code
% close all; clear all; clc; 
% Dt = 0.04;  % time step
% N = 1/Dt;   % number of time domain data points
% k = 0:1:N-1;
% n = 0:1:N-1;
% Xk = zeros(1,max(k));
% x = zeros(1,N);
% x(1:12) = 5;
% ft = fft(x);
% figure; stem(n,abs(ft)/25); grid on;


%% Problem 3 - part a
Dt = 0.04;  % time step
N = 1/Dt;   % number of time domain data points
k = 0:1:N-1;
n = 0:1:N-1;
Xk = zeros(1,N);
x = zeros(1,N); x(1:12) = 5;

for j = 1:N
    Xk(j)  = (1/N)*x(1)*exp(-(1i*2*pi*k(j)*n(1))/N);
    for ni = 2:N
        Xk(j)  = Xk(j) + (1/N)*x(ni)*exp(-(1i*2*pi*k(j)*n(ni))/N);
    end
 end
figure; stem(k,abs(Xk)); xlim([0,(N-1)/2]); grid on;
xlabel('Frequency domain'); ylabel('Magnitude');
title('Magnitude frequency domain for \Delta t = 0.04');
print('fn3_VB1_1','-dpng');

X = zeros(1,(N-1)/2 + 1); xn = zeros(1,N); t = 0:0.04:1; t(1) = [];
X(1) = Xk(1);
for i = 1:N
    xn(i) = X(1);
    for j = 1:(N-1)/2
        xn(i) = xn(i) + 2*(real(Xk(j+1))*cos((2*pi*j*n(i))/N) - imag(Xk(j+1))*sin((2*pi*j*n(i))/N));
    end
end
xn = [5,xn]; t = [0,t]; x = [5,x];
figure; stem(t,xn); grid on; hold on; plot(t,x);
xlabel('Time domain'); ylabel('Signal');
title('Inverse Discrete Fourier Transform for \Delta t = 0.04 s');
print('fn3_VB1_2','-dpng');

%% 
%% Problem 3 - part b
Dt2 = 0.06;  % time step
N2 = ceil(1/Dt2);   % number of time domain data points
k2 = 0:1:N2-1;
n2 = 0:1:N2-1;
Xk2 = zeros(1,N2);
x2 = zeros(1,N2); x2(1:N2/2) = 5;
for j = 1:N2
    Xk2(j)  = (1/N2)*x2(1)*exp(-(1i*2*pi*k2(j)*n2(1))/N2);
    for ni = 2:N2
        Xk2(j)  = Xk2(j) + (1/N2)*x2(ni)*exp(-(1i*2*pi*k2(j)*n2(ni))/N2);
    end
 end

figure; stem(k2,abs(Xk2)); xlim([0,(N2-1)/2]); grid on;
xlabel('Frequency domain'); ylabel('Magnitude');
title('Magnitude frequency domain for \Delta t = 0.06');
print('fn3_VB2_1','-dpng');

X2 = zeros(1,(N2-1)/2 + 1); xn2 = zeros(1,N2); t2 = 0:0.06:1; 
X2(1) = Xk2(1);
for i = 1:N2
    xn2(i) = X2(1);
    for j = 1:(N2-1)/2
        xn2(i) = xn2(i) + 2*(real(Xk2(j+1))*cos((2*pi*j*n2(i))/N2) - imag(Xk2(j+1))*sin((2*pi*j*n2(i))/N2));
    end
end
figure; stem(t2,xn2); grid on; hold on; plot(t2,x2);
xlabel('Time domain'); ylabel('Signal');
title('Inverse Discrete Fourier Transform for \Delta t = 0.06 s');
print('fn3_VB2_2','-dpng');
k3 = k2(1:(N2+1)/2); Xk3 = Xk2(1:(N2+1)/2);
figure; stem(k,abs(Xk)); xlim([0,(N-1)/2]); hold on; stem(k3, abs(Xk3),'r--');
xlabel('Frequency domain'); ylabel('Magnitude');
title('Magnitude frequency domain for \Delta t = 0.04s, \Delta t = 0.06s');
legend('\Delta t = 0.04','\Delta t = 0.06');
grid on; print('fn3_VB2_3','-dpng');

figure; stem(t,xn); hold on; plot(t2,x2); stem(t2,xn2,'r--');
xlabel('Time domain'); ylabel('Signal');
title('Inverse Discrete Fourier Transform for \Delta t = 0.04s, \Delta t = 0.06s');
legend('\Delta t = 0.04','\Delta t = 0.06');
grid on; print('fn3_VB2_4','-dpng');

