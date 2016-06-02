% Luan Cong Doan - Final Exam - Vibration - Question 4
close all; clear all; clc;

F = 400; T = 1/400; w=0:1:400; wx = w*2*pi;
k1 = 2; k2 = 10;
% for 2 average - the magnitude and coherence
for i=1:k1
    s2(i,:) = randn(401,1) + 1i*randn(401,1);
    m2(i,:) = randn(401,1) + 1i*randn(401,1);
    n2(i,:) = randn(401,1) + 1i*randn(401,1);
end
U2 = 5*s2; 
for j = 1:k1
    for i = 1:length(wx)
        V2(j,i) = -wx(i)^2*U2(j,i)*(1/(-wx(i)^2 + 1i*100*wx(i) + 1.14*10^6) + 1/(-wx(i)^2 + 1i*80*wx(i) + 2.1*10^6));
        H_true(i) = -wx(i)^2*(1/(-wx(i)^2 + 1i*100*wx(i) + 1.14*10^6) + 1/(-wx(i)^2 + 1i*80*wx(i) + 2.1*10^6));
    end
end

for i=1:length(wx)
        Gsx2(i) = conj(s2(1,i))*(n2(1,i)+V2(1,i)) + conj(s2(2,i))*(n2(2,i)+V2(2,i));
        Gsf2(i) = conj(s2(1,i))*(m2(1,i)+U2(1,i)) + conj(s2(2,i))*(m2(2,i)+U2(2,i));
        Hc2(i) = abs(Gsx2(i)/Gsf2(i));
        
        Guu2(i) = conj(U2(1,i))*U2(1,i) + conj(U2(2,i))*U2(2,i);
        Gvv2(i) = conj(V2(1,i))*V2(1,i) + conj(V2(2,i))*V2(2,i);
        Gmm2(i) = conj(m2(1,i))*m2(1,i) + conj(m2(2,i))*m2(2,i);
        Gnn2(i) = conj(n2(1,i))*n2(1,i) + conj(n2(2,i))*n2(2,i);
        Gff2(i) = Guu2(i) + Gmm2(i);
        Gxx2(i) = Gvv2(i) + Gnn2(i);
        Gfx2(i) = conj(U2(1,i))*V2(1,i)+ conj(U2(2,i))*V2(2,i);
        Gamma2(i) = (abs(Gfx2(i))^2)/(abs(Gxx2(i))*abs(Gff2(i)));
end
% for magnitude
figure; plot(w,Hc2); grid on; grid minor; xlabel('Frequency'); ylabel('Amplitude');
title('2 average with true value FRF'); %print('fn4_VB1_1','-dpng');
% for the cohenrence
figure; plot(w,Gamma2); grid on; grid minor; xlabel('Frequency'); ylabel('\gamma^2_{FX}');
title('2 average with coherence \gamma^2_{FX}'); %print('fn4_VB2_1','-dpng');
% for single-sided spectral density
figure; plot(w,Gff2); grid on; grid minor; xlabel('Frequency'); ylabel('G_{FF}');
title('Single-sided spectral densities G_{FF} - 2 average'); %print('fn4_VB3_1','-dpng');
figure; plot(w,Gxx2); grid on; grid minor; xlabel('Frequency'); ylabel('G_{XX}');
title('Single-sided spectral densities G_{XX} - 2 average'); %print('fn4_VB3_2','-dpng');

% for 10 average
for i=1:k2
    s10(i,:) = randn(401,1) + 1i*randn(401,1);
    m10(i,:) = randn(401,1) + 1i*randn(401,1);
    n10(i,:) = randn(401,1) + 1i*randn(401,1);
end
U10 = 5*s10; 
for j = 1:k2
    for i = 1:length(wx)
        V10(j,i) = -wx(i)^2*U10(j,i)*(1/(-wx(i)^2 + 1i*100*wx(i) + 1.14*10^6) + 1/(-wx(i)^2 + 1i*80*wx(i) + 2.1*10^6));
    end
end

for i=1:length(wx)
    Gsx10(i) = conj(s10(1,i))*(n10(1,i)+V10(1,i));
    Gsf10(i) = conj(s10(1,i))*(m10(1,i)+U10(1,i));
    Gfx10(i) = conj(n10(1,i))*V10(1,i) + conj(m10(1,i))*V10(1,i) + conj(U10(1,i))*n10(1,i) + conj(m10(1,i))*n10(1,i);
    Gff10(i) = conj(U10(1,i))*U10(1,i) + conj(m10(1,i))*U10(1,i) + conj(U10(1,i))*m10(1,i) + conj(m10(1,i))*m10(1,i);                
    Guu10(i) = conj(U10(1,i))*U10(1,i);
    Gvv10(i) = conj(V10(1,i))*V10(1,i);
    Gmm10(i) = conj(m10(1,i))*m10(1,i);
    Gnn10(i) = conj(n10(1,i))*n10(1,i);
    Gfx10(i) = conj(U10(1,i))*V10(1,i);
    
    for j = 2:k2
        Gsx10(i) = Gsx10(i) + conj(s10(j,i))*(n10(j,i)+V10(j,i)); 
        Gsf10(i) = Gsf10(i) + conj(s10(j,i))*(m10(j,i)+U10(j,i));
        Guu10(i) = Guu10(i) + conj(U10(j,i))*U10(j,i);
        Gvv10(i) = Gvv10(i) + conj(V10(j,i))*V10(j,i);
        Gmm10(i) = Gmm10(i) + conj(m10(j,i))*m10(j,i);
        Gnn10(i) = Gnn10(i) + conj(n10(j,i))*n10(j,i);
        Gfx10(i) = Gfx10(i) + conj(U10(j,i))*V10(j,i);
    end
    Hc10(i) = abs(Gsx10(i)/Gsf10(i));
    Gff10(i) = Guu10(i) + Gmm10(i);
    Gxx10(i) = Gvv10(i) + Gnn10(i);
    Gamma10(i) = (abs(Gfx10(i))^2)/abs(Gxx10(i)*Gff10(i));
end

figure; plot(w,Hc10); grid on; grid minor; xlabel('\omega'); ylabel('Amplitude');
title('10 average with true value FRF'); %print('fn4_VB1_2','-dpng');
figure; plot(w,Gamma10); grid on; grid minor; xlabel('Frequency'); ylabel('\gamma^2_{FX}');
title('10 average with coherence \gamma^2_{FX}'); %print('fn4_VB2_2','-dpng');
% for single-sided spectral density
figure; plot(w,Gff10/5); grid on; grid minor; xlabel('Frequency'); ylabel('G_{FF}');
title('Single-sided spectral densities G_{FF} - 10 average'); %print('fn4_VB3_3','-dpng');
figure; plot(w,Gxx10/5); grid on; grid minor; xlabel('Frequency'); ylabel('G_{XX}');
title('Single-sided spectral densities G_{XX} - 10 average'); %print('fn4_VB3_4','-dpng');

figure; plot(w,abs(H_true),'k', 'linewidth',1); hold on; plot(w,Hc2,'r--'); plot(w,Hc10,'b-.'); grid on;
xlabel('Frequency'); ylabel('H_c'); title('The Frequency response: 2-10 average and true value');
legend('true value','2 average','10 average'); print('fn4_VB4_1','-dpng');



