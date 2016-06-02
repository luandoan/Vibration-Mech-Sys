% Luan Cong Doan - Vibration - Final Project
close all; clear all; clc;

%% Problem 1 - part 1
x1_1 = 0:0.01:5;
x1_2 = 0:0.01:10.5;
m = 14.69/97.4;
% m = 52.4/97.4;

fx1 = 1 + cos(x1_1).*cosh(x1_1) + m.*x1_1.*cos(x1_1).*sinh(x1_1) - m.*x1_1.*cosh(x1_1).*sin(x1_1);
fx2 = 1 + cos(x1_2).*cosh(x1_2) + m.*x1_2.*cos(x1_2).*sinh(x1_2) - m.*x1_2.*cosh(x1_2).*sin(x1_2);
fn1 = @(n) 1 + cos(n).*cosh(n) + m.*n.*cos(n).*sinh(n) - m.*n.*cosh(n).*sin(n);
s1 = fzero(fn1,1.5);
s2 = fzero(fn1,4);
s3 = fzero(fn1,7);
s4 = fzero(fn1,10.2);
figure; plot(x1_1,fx1); grid on;
xlabel('\beta L'); ylabel('f(\beta L)'); title('Charateristic equation');
% print('fn1_VB1_1','-dpng');
figure; plot(x1_2,fx2); grid on;
xlabel('\beta L'); ylabel('f(\beta L)'); title('Charateristic equation');
% print('fn1_VB1_2','-dpng');

s = [s1,s2,s3,s4];

for i = 1:4
    w(i) = s(i).^2*0.7;
    frac(i) = -(sinh(s(i)) + sin(s(i)))/(cosh(s(i)) + cos(s(i)));
end

X = 0:0.01:0.46;
for i = 1:4
    for j = 1:length(X)
    mode(i,j) = (frac(i)*(cosh(s(i)*X(j)/0.46) - cos(s(i)*X(j)/0.46)) + sinh(s(i)*X(j)/0.46) - sin(s(i)*X(j)/0.46));
    end;
    figure; plot(X,-mode(i,:)); grid on;
    xlabel('x'); ylabel(['Y_',i]); title(['Mode Shape 1-',num2str(i)]);
    print(['fn1_VB1_mode',num2str(i)],'-dpng');
end


%% Problem 1 - part 2
% close all; clear all; clc;
x2_1 = 0:0.1:3.5;
x2_2 = 0:0.1:9;
for i =1:length(x2_1)
    k1a(i) = 0.1508*x2_1(i); 
    k2a(i) = 0.06*x2_1(i)^3;
    f2a(i) = (k1a(i)*k2a(i)-1)*(cosh(x2_1(i))*cos(x2_1(i))) - (k1a(i)*k2a(i)+1) + (k1a(i)+k2a(i))*(cosh(x2_1(i))*sin(x2_1(i))) - (k1a(i)-k2a(i))*(sinh(x2_1(i))*sin(x2_1(i)));
end
for i =1:length(x2_2)
    k1b(i) = 0.1508*x2_2(i); 
    k2b(i) = 0.06*x2_2(i)^3;
    f2b(i) = (k1b(i)*k2b(i)-1)*(cosh(x2_2(i))*cos(x2_2(i))) - (k1b(i)*k2b(i)+1) + (k1b(i)+k2b(i))*(cosh(x2_2(i))*sin(x2_2(i))) - (k1b(i)-k2b(i))*(sinh(x2_2(i))*sin(x2_2(i)));
end
figure; plot(x2_1,f2a); grid on;
xlabel('\beta L'); ylabel('f(\beta L)'); title('Charateristic equation 2');
% print('fn1_VB2_1','-dpng');
figure; plot(x2_2,f2b); grid on;
xlabel('\beta L'); ylabel('f(\beta L)'); title('Charateristic equation 2');
% print('fn1_VB2_2','-dpng');

f2x = @(x)(0.1508*x*0.06*x^3-1)*(cosh(x)*cos(x)) - (0.1508*x*0.06*x^3+1) + (0.1508*x+0.06*x^3)*(cosh(x)*sin(x)) - (0.1508*x-0.06*x^3)*(sinh(x)*sin(x));
s21 = fzero(f2x,1.5);
s22 = fzero(f2x,3.2);
s23 = fzero(f2x,6);
s24 = fzero(f2x,9);

s2 = [s21,s22,s23,s24];
for i = 1:4
    w2(i) = s2(i).^2*22.05;
    frac2(i) = -(0.1508*s2(i)*(sinh(s2(i)) - sin(s2(i))) + (cosh(s2(i)) + cos(s2(i))))/(0.1508*s2(i)*(cosh(s2(i)) + cos(s2(i))) + (sinh(s2(i))-sin(s2(i))));
end

X = 0:0.01:0.46;
for i = 1:4
    for j = 1:length(X)
    mode2(i,j) = frac2(i)*(cosh(s2(i)*X(j)/0.46) - cos(s2(i)*X(j)/0.46)) + sinh(s2(i)*X(j)/0.46) - sin(s2(i)*X(j)/0.46);
    end;
    figure; plot(X,-mode2(i,:)); grid on;
    xlabel('x'); ylabel(['Y_',i]); title(['Mode Shape 2-',num2str(i)]);
    print(['fn1_VB2_mode',num2str(i)],'-dpng');
end


% %% ## Problem 2 ## %%
% close all; clear all; clc;
% %% Part 1
% x3 = 0:0.01:5;
% fx3 = cosh(x3).*cos(x3) + 1;
% figure; plot(x3,fx3,'k'); grid on;
% xlabel('\beta L'); ylabel('f(\beta L)'); title('Characteristic equation');
% print('fn2_VB1_1','-dpng');
% 
% fn3 = @(n) cosh(n)*cos(n) + 1;
% s31 = fzero(fn3,2); s32 = fzero(fn3,4.5);
% s3 = [s31,s32];
% for i=1:length(s3)
%     w3(i) = s3(i).^2*0.7;
%     frac3(i) = -(sinh(s3(i)) + sin(s3(i)))/(cosh(s3(i)) + cos(s3(i)));
% end
% X = 0:0.01:0.46;
% for i = 1:2
%     for j = 1:length(X)
%     mode(i,j) = frac3(i)*(cosh(s3(i)*X(j)/0.46) - cos(s3(i)*X(j)/0.46)) + sinh(s3(i)*X(j)/0.46) - sin(s3(i)*X(j)/0.46);
%     end;
%     figure; plot(X,mode(i,:)); grid on;
%     xlabel('x'); ylabel(['Y_',i]); title(['Mode Shape 2-',num2str(i)]);
%     print(['fn2_VB1_mode',num2str(i)],'-dpng');
% end
% 

