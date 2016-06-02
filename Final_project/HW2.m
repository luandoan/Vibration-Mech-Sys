%% Vibration Final Project
%%

%% Problem 1
% *point mass assumption*
clear all;
close all;
clc
E = 68.9*10^9;%elastic modulus in pa
d = 2.7*10^-3*10^6;%density in kg/m^3
Db = 8.73*10^-3*10^6;%density of brass
b = 2.45*10^-2;%width of beam in meter
h = 3.2*10^-3;%thickness of the beam in meter
L = 46*10^-2;%length of the beam in meter
I = 1/12*b*h^3; %moment of inertia of the beam
% Transcendental Equation
% cosh(Beta*L)*cos(Beta*L)+1+ M/m*(Beta*L)*(cos(Beta*L)*sinh(Beta*L)-sin(Beta*L)*cosh(Beta*L)) = 0
V_brass = pi*(2.3/2)^2*1.13+pi*(0.9/2)^2*(2.3-1.13);%volume of the tip mass
V_endbeam = 2.3*2.45*3.2*10^-1; %volume of the end beam
M = Db*10^-6*V_brass + d*10^-6*V_endbeam; %mass of the tip mass
m = d*10^-6*3.2*10^-1*2.45*46; %mass of the beam
%function used to find root of Beta*L
f = @(x)cosh(x)*cos(x)+1+ M/m*(x)*(cos(x)*sinh(x)-sin(x)*cosh(x));
x1 = fzero(f,[0,2]);
x2 = fzero(f,[4,6]);
x3 = fzero(f,[6,8]);
x4 = fzero(f,[10,12]);

% natural frequencies in Hz
w1 = x1^2*sqrt(E*I/(d*b*h*L^4))/(2*pi);
w2 = x2^2*sqrt(E*I/(d*b*h*L^4))/(2*pi)
w3 = x3^2*sqrt(E*I/(d*b*h*L^4))/(2*pi)
w4 = x4^2*sqrt(E*I/(d*b*h*L^4))/(2*pi)

% Let H = A/B
H1 = -(sinh(x1)+sin(x1))/(cosh(x1)+cos(x1));
H2 = -(sinh(x2)+sin(x2))/(cosh(x2)+cos(x2));
H3 = -(sinh(x3)+sin(x3))/(cosh(x3)+cos(x3));
H4 = -(sinh(x4)+sin(x4))/(cosh(x4)+cos(x4));

X = 0:0.01:0.46;
for i = 1:length(X)
    y1(i) = H1*(cosh(x1*X(i)/L)-cos(x1*X(i)/L))+sinh(x1*X(i)/L)-sin(x1*X(i)/L);
    y2(i) = H2*(cosh(x2*X(i)/L)-cos(x2*X(i)/L))+sinh(x2*X(i)/L)-sin(x2*X(i)/L);
    y3(i) = H3*(cosh(x3*X(i)/L)-cos(x3*X(i)/L))+sinh(x3*X(i)/L)-sin(x3*X(i)/L);
    y4(i) = H4*(cosh(x4*X(i)/L)-cos(x4*X(i)/L))+sinh(x4*X(i)/L)-sin(x4*X(i)/L);
end
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
plot(X/L,-y1)
title('Mode 1')
xlabel('X/L'),ylabel('mode shape')
grid on 
grid minor
subplot(2,2,2)
plot(X/L,-y2)
title('Mode 2')
xlabel('X/L'),ylabel('mode shape')
grid on 
grid minor
subplot(2,2,3)
plot(X/L,-y3)
title('Mode 3')
xlabel('X/L'),ylabel('mode shape')
grid on 
grid minor
subplot(2,2,4)
plot(X/L,-y4)
title('Mode 4')
xlabel('X/L'),ylabel('mode shape')
grid on 
grid minor
%%
% *Solid having mass and inertia assumption

Mbot = pi*(2.3/2*10^-2)^2*1.13*10^-2*Db;%mass of the bottom cylinder
Mtop = pi*(0.9/2*10^-2)^2*(2.3-1.13)*10^-2*Db; %mass of the top cylinder
I0 = 1/4*Mbot*(2.3/2*10^-2)^2 + 1/3*Mbot*(1.13*10^-2)^2 + Mtop*(2.3*10^-2)^2;
syms x
K1 = M/m*x;
K2 = I0*x^3/(m*L^2);
%function used to find root of Beta*L
q =@(x)(M/m*x*I0*x^3/(m*L^2)-1)*(cosh(x)-cos(x))-(M/m*x*I0*x^3/(m*L^2)+1)+(M/m*x+I0*x^3/(m*L^2))*(cosh(x)*sin(x))-(M/m*x-I0*x^3/(m*L^2))*(sinh(x)*cos(x));

%find four roots of Beta*L
q1 = fzero(q,[0,3]);
q2 = fzero(q,[3,6]);
q3 = fzero(q,[7,8]);
q4 = fzero(q,[8,11]);

%natural frequency in Hz
f1 = q1^2/(2*pi)*sqrt(E*I/(m*L^3))
f2 = q2^2/(2*pi)*sqrt(E*I/(m*L^3))
f3 = q3^2/(2*pi)*sqrt(E*I/(m*L^3))
f4 = q4^2/(2*pi)*sqrt(E*I/(m*L^3))

% let R = A/B
R1 = -(M/m*q1*(sinh(q1)-sin(q1))+(cosh(q1)+cos(q1)))/(M/m*q1*(cosh(q1)-cos(q1))+(sinh(q1)-sin(q1)));
R2 = -(M/m*q2*(sinh(q2)-sin(q2))+(cosh(q2)+cos(q2)))/(M/m*q2*(cosh(q2)-cos(q2))+(sinh(q2)-sin(q2)));
R3 = -(M/m*q3*(sinh(q3)-sin(q3))+(cosh(q3)+cos(q3)))/(M/m*q3*(cosh(q3)-cos(q3))+(sinh(q3)-sin(q3)));
R4 = -(M/m*q4*(sinh(q4)-sin(q4))+(cosh(q4)+cos(q4)))/(M/m*q4*(cosh(q4)-cos(q4))+(sinh(q4)-sin(q4)));

Q = 0:0.01:0.46;
for i = 1:length(Q)
    yr1(i) = R1*(cosh(q1*Q(i)/L)-cos(q1*Q(i)/L))+sinh(q1*Q(i)/L)-sin(q1*Q(i)/L);
    yr2(i) = R2*(cosh(q2*Q(i)/L)-cos(q2*Q(i)/L))+sinh(q2*Q(i)/L)-sin(q2*Q(i)/L);
    yr3(i) = R3*(cosh(q3*Q(i)/L)-cos(q3*Q(i)/L))+sinh(q3*Q(i)/L)-sin(q3*Q(i)/L);
    yr4(i) = R4*(cosh(q4*Q(i)/L)-cos(q4*Q(i)/L))+sinh(q4*Q(i)/L)-sin(q4*Q(i)/L);
end

%plot four mode shapes
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
plot(Q/L,-yr1)
title('Mode 1')
xlabel('X/L'),ylabel('mode shape')
grid on 
grid minor
subplot(2,2,2)
plot(Q/L,-yr2)
title('Mode 2')
xlabel('X/L'),ylabel('mode shape')
grid on 
grid minor
subplot(2,2,3)
plot(Q/L,-yr3)
title('Mode 3')
xlabel('X/L'),ylabel('mode shape')
grid on 
grid minor
subplot(2,2,4)
plot(Q/L,-yr4)
title('Mode 4')
xlabel('X/L'),ylabel('mode shape')
grid on 
grid minor

%% Problem 2
% *part (a)*

%Transcendental equation
%F(Beta*L) = cosh(Beta*L)*cos(Beta*L)+1 =0
%function used to find root of Beta*L
f = @(x)cosh(x)*cos(x)+1;
root1 = fzero(f,[0,2]);
root2 = fzero(f,[4,6]);
% natural frequencies in Hz
nf1 = root1^2*sqrt(E*I/(d*b*h*L^4))/(2*pi)
nf2 = root2^2*sqrt(E*I/(d*b*h*L^4))/(2*pi)

% Let P = A/B
P1 = -(sinh(root1)+sin(root1))/(cosh(root1)+cos(root1));
P2 = -(sinh(root2)+sin(root2))/(cosh(root2)+cos(root2));
X = 0:0.01:0.46;
for i = 1:length(X)
    y1(i) = P1*(cosh(root1*X(i)/L)-cos(root1*X(i)/L))+sinh(root1*X(i)/L)-sin(root1*X(i)/L);
    y2(i) = P2*(cosh(root2*X(i)/L)-cos(root2*X(i)/L))+sinh(root2*X(i)/L)-sin(root2*X(i)/L);
end
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
plot(X/L,-y1)
title('Mode 1')
xlabel('X/L'),ylabel('mode shape')
grid on 
grid minor
subplot(2,1,2)
plot(X/L,-y2)
title('Mode 2')
xlabel('X/L'),ylabel('mode shape')
grid on 
grid minor
%%
% *part (b)*

%number of elements
% num_element = 3;
% num_modes = 2;
%length of each element
% l_element = L/num_element;
% X = 0:0.01:0.46;
%number of nodes
% n = 1:num_element+1;
% node1 = 1:num_element;
% node2 = 2:num_element+1;
% mx = num_element+1;
% k = zeros(2*mx);
% m = zeros(2*mx);
cross_area = b*h;
M=((d*cross_area*L)/1260)*[312 0 54 -13*L/3 0 0;...
             0 8*L^2/9 13*L/3 -L^2/3 0 0;...
             54 13*L/3 312 0 54 -13*L/3;...
             -13*L/3 -L^2/3 0 8*L^2/9 13*L/3 -L^2/3;...
             0 0 54 13*L/3 156 -22*L/3;...
             0 0 -13*L/3 -L^2/3 -22*L/3 4*L^2/9];
K=((27*E*I)/(L^3))*[24 0 -12 2*L 0 0;...
    0 8*L^2/9 -2*L 2*L^2/9 0 0;...
    -12 -2*L 24 0 -12 2*L;...
    2*L 2*L^2/9 0 8*L^2/9 -2*L 2*L^2/9;...
    0 0 -12 -2*L 12 -2*L;...
    0 0 2*L 2*L^2/9 -2*L 4*L^2/9];
% Ke = (E*I)/(l_element^3)*[12,6*l_element,-12,6*l_element;...
%             6*l_element,4*l_element^2,-6*l_element,2*l_element^2;...
%             -12,-6*l_element,12,-6*l_element;...
%             6*l_element,2*l_element^2,-6*l_element,4*l_element^2];
% Me = (mass_per_length/420)*l_element*[156,-22*l_element,54,13*l_element;...
%                            -22*l_element,4*l_element^2,-13*l_element,-3*l_element^2;...
%                            54,-13*l_element,156,22*l_element;...
%                            13*l_element,-3*l_element^2,22*l_element,4*l_element^2];
% Kg = zeros(2*(num_element+1),2*(num_element+1));
% Mg = zeros(2*(num_element+1),2*(num_element+1));
% for i=1:num_element
%     Kg(2*i-1:2*(i+1),2*i-1:2*(i+1)) = Kg(2*i-1:2*(i+1),2*i-1:2*(i+1))+Ke;
%     Mg(2*i-1:2*(i+1),2*i-1:2*(i+1)) = Mg(2*i-1:2*(i+1),2*i-1:2*(i+1))+Me;
% end
% Kg(1:2,:) = [];
% Kg(:,1:2) = [];
% Mg(1:2,:) = [];
% Mg(:,1:2) = [];
% Q = zeros(2*num_element,1);
Akm = inv(M)*K;
[V,gama] = eig(Akm);
Wkm1 = sqrt(gama(6,6))/(2*pi)
Wkm2 = sqrt(gama(5,5))/(2*pi)
Xe = linspace(0,L,4);
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
plot(Xe/L,[0,-V(1,6),-V(3,6),-V(5,6)])
title('Mode 1')
xlabel('X/L'),ylabel('mode shape')
grid on 
grid minor
subplot(2,1,2)
plot(Xe/L,[0,-V(1,5),-V(3,5),-V(5,5)])
title('Mode 2')
xlabel('X/L'),ylabel('mode shape')
grid on 
grid minor
%% Problem 3
% *25 Hz sampling rate*
clear all
% part(a)
fs = 25; %sampling frequency
ts = 1/25;
tau = 1; %sampling period
N = fs*tau;
ta = 0.04:1/fs:1;
for i = 1:length(ta)
    xa(1:12) = 5;
    xa(13:25) = 0;
end
fr = fs/N*[0:N-1];
Y = fft(xa);

figure('units','normalized','outerposition',[0 0 1 1])
stem(fr,abs(Y)/25)
xlim([0,fs/2])
title('DFT response at 25 Hz sampling rate')
xlabel('Frequency (Hz)'),ylabel('Magnitude')
grid on
grid minor
% part(b)
iY = real(ifft(Y));
figure('units','normalized','outerposition',[0 0 1 1])
stem(ta,iY,'r','LineWidth',2)
hold on
plot(ta,xa,'b-.')
% ylim([0,5])
title('Sampling rate at 25 Hz')
legend('IDFT','Origin square wave')
xlabel('Time (seconds)'),ylabel('x(t)')
grid on
grid minor

%%
% *16.67 Hz sampling rate*

fs2 = 16.67;
ts2 = 1/16.67;
N = round(fs2*tau);
tb = 0:ts2:1;
for i=1:length(tb)
    xb(1:9) = 5;
    xb(10:17) = 0;
end
frb = 0:N-1;
Yb = fft(xb);

figure('units','normalized','outerposition',[0 0 1 1])
stem(frb,abs(Yb))
xlim([0,fs2/2])
title('DFT response at 16.67 Hz sampling rate')
xlabel('Frequency (Hz)'),ylabel('Magnitude')
grid on
grid minor

iYb = ifft(Yb,10000);
figure('units','normalized','outerposition',[0 0 1 1])
plot(1/10000:1/10000:1,iYb*(10000/N),'r','LineWidth',2)
hold on
plot(tb,xb,'b-.')
legend('IDFT','Origin square wave')
title('sampling rate at 16.67 Hz')
xlabel('Time (seconds)'),ylabel('x(t)')
grid on
grid minor
%%
%part(a)
figure('units','normalized','outerposition',[0 0 1 1])
stem(fr,abs(Y),'r')
hold on
stem(frb,abs(Yb))
xlim([0,fs/2])
legend('25 Hz sampling rate','16.67 Hz sampling rate')
xlabel('Frequency (Hz)'),ylabel('Magnitude')
grid on
grid minor
%%
N = 25;
n = 0:1:24;
k = 0:1:24;
x(1:12) = 5;
x(13:25) = 0;
for j = 1:length(k)
    X(j) = (1/N)*x(1)*exp(-(1i*2*pi*k(j)*n(1))/N);
   for i = 2:length(k)
         X(j) = X(j) + (1/N)*x(i)*exp((-1i*2*pi*k(j)*n(i))/N);
   end
end
stem(n,abs(X))
t = 0.04:0.04:1;
for i = 1:length(k)
    xn(i) = X(1);
    for j = 1:(length(k)-1)/2
         xn(i) = xn(i)+2*(real(X(j+1))*cos((2*pi*j*n(i))/N) - imag(X(j+1))*sin((2*pi*j*n(i))/N));
    end
end
figure
stem(t,xn)

% N = 12;
% Np = 6;
% Xk0 = 2.5;
% k = [1:1:12];
% for i = 1:length(k)
%     Xk(i) = 2.5*1/N*(sin(k(i)*Np*pi/N)/sin(k(i)*pi/N));
% end
% k_total =0:1:12;
% Xk_total = [Xk0,Xk];
% figure('units','normalized','outerposition',[0 0 1 1])
% stem(k_total,Xk_total)
% xlabel('Frequency (Hz)'),ylabel('Magnitude')
% grid on
% grid minor




%% Problem 4
% *Part(a)*
F = 400;
T = 1/400;
S_r = randn(401,1) + i*randn(401,1);
M_r = randn(401,1) + i*randn(401,1);
N_r = randn(401,1) + i*randn(401,1);
U = 5*S_r;
syms wr
V_r = -wr^2*U*(1/(-wr^2+i*100*wr+1.14*10^6)+1/(-wr^2+i*80*wr+2.1*10^6));
V_ravg2 = sum(V_r)/2;
V_ravg10 = sum(V_r)/10;
U_avg2 = sum(U)/2;
U_avg10 = sum(U)/10;
H_trueAvg2 = V_ravg2/U_avg2;
H_trueAvg10 = V_ravg10/U_avg10;
wr = 0:2*pi*F;
Fr = linspace(0,400,length(wr));
H_trueAvg2 = subs(H_trueAvg2);
H_trueAvg10 = subs(H_trueAvg10);
% figure('units','normalized','outerposition',[0 0 1 1])
figure; subplot(2,1,1)
plot(Fr,abs(H_trueAvg2))
title('FRF Hc with 2 averages')
xlabel('Frequency (Hz)'),ylabel('Magnitude Hc')
grid on
grid minor
subplot(2,1,2)
plot(Fr,abs(H_trueAvg10),'r')
title('FRF Hc with 10 averages')
xlabel('Frequency (Hz)'),ylabel('Magnitude Hc')
grid on
grid minor
% figure; plot(Fr,abs(H_trueAvg2)); hold on; plot(Fr,abs(H_trueAvg10),'r--'); grid on;

%%
% *part(b)*
syms wr
Guv_avg2 = sum(conj(U).*V_r)/2;
Gvv_avg2 = sum(conj(V_r).*V_r)/2;
Gnn_avg2 = sum(conj(N_r).*N_r)/2;
Guu_avg2 = sum(conj(U).*U)/2;
Gmm_avg2 = sum(conj(M_r).*M_r)/2;

Guv_avg10 = sum(conj(U).*V_r)/10;
Gvv_avg10 = sum(conj(V_r).*V_r)/10;
Gnn_avg10 = sum(conj(N_r).*N_r)/10;
Guu_avg10 = sum(conj(U).*U)/10;
Gmm_avg10 = sum(conj(M_r).*M_r)/10;

coherence_avg2 = (abs(Guv_avg2))^2/((Gvv_avg2+Gnn_avg2)*(Guu_avg2+Gmm_avg2));
coherence_avg10 = (abs(Guv_avg10))^2/((Gvv_avg10+Gnn_avg10)*(Guu_avg10+Gmm_avg10));

wr = 0:2*pi*F;
Fr = linspace(0,400,length(wr));
coherence_avg2 = subs(coherence_avg2);
coherence_avg10 = subs(coherence_avg10);
figure; %('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
plot(Fr,coherence_avg2)
title('Coherence with 2 averages')
xlabel('Frequency (Hz)'),ylabel('Coherence')
grid on
grid minor
subplot(2,1,2)
plot(Fr,coherence_avg10,'r')
title('Coherence with 10 averages')
xlabel('Frequency (Hz)'),ylabel('Coherence')
grid on
grid minor;
% figure; plot(Fr,coherence_avg2); hold on; plot(Fr,coherence_avg10,'r--'); grid on;


%%
% *part(c)*

syms wr
Gff_avg2 = 2*(Guu_avg2+Gmm_avg2);
Gxx_avg2 = 2*(Gvv_avg2+Gnn_avg2);
Gff_avg10 = 2*(Guu_avg10+Gmm_avg10);
Gxx_avg10 = 2*(Gvv_avg10+Gnn_avg10);

wr = 0:2*pi*F;
Fr = linspace(0,400,length(wr));
Gff_avg2 = subs(Gff_avg2);
Gxx_avg2 = subs(Gxx_avg2);
Gff_avg10 = subs(Gff_avg10);
Gxx_avg10 = subs(Gxx_avg10);

figure; %('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
plot([0,400],[double(Gff_avg2),double(Gff_avg2)])
title('single_sided density of Gff with 2 averages')
xlabel('Frequency (Hz)'),ylabel('Gff')
grid on
grid minor
subplot(2,2,2)
plot(Fr,Gxx_avg2)
title('single_sided density of Gxx with 2 averages')
xlabel('Frequency (Hz)'),ylabel('Gxx')
grid on
grid minor
subplot(2,2,3)
plot([0,400],[double(Gff_avg10),double(Gff_avg10)])
title('single_sided density of Gff with 10 averages')
xlabel('Frequency (Hz)'),ylabel('Gff')
grid on
grid minor
subplot(2,2,4)
plot(Fr,Gxx_avg10)
title('single_sided density of Gxx with 10 averages')
xlabel('Frequency (Hz)'),ylabel('Gxx')
grid on
grid minor
