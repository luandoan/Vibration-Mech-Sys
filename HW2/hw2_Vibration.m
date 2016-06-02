% Luan Cong Doan - Vibration HW2
close all; clear all; clc;
%%
% x = 0:0.1:10.5;
% % m = 14.69/97.4;
% f = 1 + cos(x).*cosh(x) + m.*x.*cos(x).*sinh(x) - m.*x.*cosh(x).*sin(x);
% 
% figure; plot(x,f); grid on;
% xlabel('\beta L'); ylabel('f(\beta L)'); title('Charateristic equation');
% % print('hw2_VB1','-dpng');
%%
%x = 0:0.01:7.5;
 m = 14.69/97.4;
% m = 52.4/97.4;

%fx = 1 + cos(x).*cosh(x) + m.*x.*cos(x).*sinh(x) - m.*x.*cosh(x).*sin(x);
fx = @(x) 1 + cos(x).*cosh(x) + m.*x.*cos(x).*sinh(x) - m.*x.*cosh(x).*sin(x);
s1 = fzero(fx,1.5);
s2 = fzero(fx,4);
s3 = fzero(fx,7);
s4 = fzero(fx,10.2);
% figure; plot(x,fx); grid on;
% xlabel('\beta L'); ylabel('f(\beta L)'); title('Charateristic equation');
% print('hw2_VB1','-dpng');

s = [s1,s2,s3,s4];

for i = 1:4
    w(i) = s(i).^2*0.7;
    frac(i) = -(sinh(s(i)) + sin(s(i)))/(cosh(s(i)) + cos(s(i)));
end

%%
X = 0:0.01:2;
% X = 0:0.01:0.46;
for i = 1:4
    for j = 1:length(X)
    mode(i,j) = (frac(i)*(cosh(s(i)*X(j)/0.46) - cos(s(i)*X(j)/0.46)) + sinh(s(i)*X(j)/0.46) - sin(s(i)*X(j)/0.46));
    end;
    figure; plot(X,mode(i,:)); grid on;
    xlabel('x'); ylabel(['Y_',i]); title(['Mode Shape ',i]);
   
end
% print('hw2_VB_mode','-dpng');

%%
close all; clear all; clc;
x2 = 0:0.1:9;
for i =1:length(x2)
    k1(i) = 0.1508*x2(i); 
    k2(i) = 0.06*x2(i)^3;
    f2(i) = (k1(i)*k2(i)-1)*(cosh(x2(i))*cos(x2(i))) - (k1(i)*k2(i)+1) + (k1(i)+k2(i))*(cosh(x2(i))*sin(x2(i))) - (k1(i)-k2(i))*(sinh(x2(i))*sin(x2(i)));
end

figure; plot(x2,f2); grid on;
xlabel('\beta L'); ylabel('f(\beta L)'); title('Charateristic equation 2');
% print('hw2_VB3','-dpng');
%%
f2x = @(x)(0.1508*x*0.06*x^3-1)*(cosh(x)*cos(x)) - (0.1508*x*0.06*x^3+1) + (0.1508*x+0.06*x^3)*(cosh(x)*sin(x)) - (0.1508*x-0.06*x^3)*(sinh(x)*sin(x));
s21 = fzero(f2x,1.5);
s22 = fzero(f2x,3.2);
s23 = fzero(f2x,6);
s24 = fzero(f2x,9);

s2 = [s21,s22,s23,s24];

%%
for i = 1:4
    w2(i) = s2(i).^2*22.05;
    frac2(i) = -(0.1508*s2(i)*(sinh(s2(i)) - sin(s2(i))) + (cosh(s2(i)) + cos(s2(i))))/(0.1508*s2(i)*(cosh(s2(i)) + cos(s2(i))) + (sinh(s2(i))-sin(s2(i))));
end

%%
X = 0:0.01:0.46;
for i = 1:4
    for j = 1:length(X)
    mode2(i,j) = frac2(i)*(cosh(s2(i)*X(j)/0.46) - cos(s2(i)*X(j)/0.46)) + sinh(s2(i)*X(j)/0.46) - sin(s2(i)*X(j)/0.46);
    end;
    figure; plot(X,mode2(i,:)); grid on;
    xlabel('x'); ylabel(['Y_',i]); title(['Mode Shape 2 ',i]);
   
end
% print('hw2_VB2_mode','-dpng');