%Luan Cong Doan - Final Exam - Vibration Question 2 

close all; clear all; clc;
%% Part 1
x3 = 0:0.01:5;
fx3 = cosh(x3).*cos(x3) + 1;
figure; plot(x3,fx3,'k'); grid on;
xlabel('\beta L'); ylabel('f(\beta L)'); 
title('Characteristic equation for Closed form solution');
print('fn2_VB1_1','-dpng');

fn3 = @(n) cosh(n)*cos(n) + 1;
s31 = fzero(fn3,2); s32 = fzero(fn3,4.5);
s3 = [s31,s32];
for i=1:length(s3)
    w3(i) = s3(i).^2*22.0499;
    frac3(i) = -(sinh(s3(i)) + sin(s3(i)))/(cosh(s3(i)) + cos(s3(i)));
end
X = 0:0.01:0.46;
for i = 1:2
    for j = 1:length(X)
    mode(i,j) = frac3(i)*(cosh(s3(i)*X(j)/0.46) - cos(s3(i)*X(j)/0.46)) + sinh(s3(i)*X(j)/0.46) - sin(s3(i)*X(j)/0.46);
    end;
    figure; plot(X,-mode(i,:)); grid on;
    xlabel('x'); ylabel(['Y_',i]); title(['Mode Shape 2-',num2str(i)]);
    print(['fn2_VB1_mode',num2str(i)],'-dpng');
end

%% Finite element solution for cantilever beam vibration
% close all; clear all; clc;

p = 2.7*10^3;                   % density of material;
A = 3.2*10^(-3)*2.54*10^(-2);   % crosssectional area of beam
l = 0.46/3;                     % length of each element
wi=0:1:600;

km = (p*A*l)/420;
kk = (68.9*10^9*6.69*10^(-11))/l^3; 
M = [312, 0, 54, -13*l, 0,0; 0, 8*l^2, 13*l, -3*l^2, 0,0; 54, 13*l, 312, 0, 54, -13*l; -13*l, -3*l^2, 0, 8*l^2, 13*l, -3*l^2; 0,0,54, 13*l, 156, -22*l; 0,0,-13*l, -3*l^2,-22*l, 4*l^2];
K = [24,0,-12,6*l,0,0; 0,8*l^2, -6*l, 2*l^2, 0,0; -12, -6*l, 24, 0, -12, 6*l; 6*l, 2*l^2, 0, 8*l^2, -6*l, 2*l^2; 0,0,-12, -6*l , 12, -6*l; 0,0, 6*l, 2*l^2, -6*l, 4*l^2];

dw = @(w) det(-km*M*w.^2 + kk*K);
for i=1:length(wi)
    dwi(i) = det(-km*M*wi(i).^2 + kk*K);
end;
figure; plot(dwi); grid on; xlabel('\omega'); ylabel('det[-\omega^2M - K]');
title('Characteristic equation for Finite Element Method');
print('fn2_VB2_1','-dpng');
w1 = fzero(dw,100);
w2 = fzero(dw,500);

MK1 = -km*M*w1.^2 + kk*K;
MK2 = -km*M*w2.^2 + kk*K;
[V1,D1] = eig(MK1); 
[V2,D2] = eig(MK2);
figure; plot([0 0.46/3 2*0.46/3 0.46], [0 -V1(1,1) -V1(3,1) -V1(5,1)]); grid on;
xlabel('x'); ylabel('Y'); title('Mode Shape 1 - Finite Element Method');
print('fn2_VB2_3','-dpng');
figure; plot([0 0.46/3 2*0.46/3 0.46], [0 V2(1,2) V2(3,2) V2(5,2)]); grid on;
xlabel('x'); ylabel('Y'); title('Mode Shape 2 - Finite Element Method');
print('fn2_VB2_4','-dpng');



