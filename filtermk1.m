clear;
close all;
% FIR band pass filter design
% Window method-kaiser window
% specifications
% Maximum passband ripple - 0.12 dB
% Minimum stopband attenuation - 54 dB
% Lower passband edge - 1000 rad/s
% Upper passband edge - 1300 rad/s
% Lower stopband edge - 900 rad/s
% Upper stopband edge - 1450 rad/s
% Sampling frequency - 3400 rad/s

%define specifications
Ap = 0.12;
Aa = 54;
omega_p1 = 1000;
omega_a1 = 900;
omega_p2 = 1300;
omega_a2 = 1450;
omega_s = 3400;
T=(2*pi)/omega_s;

%find the transition bandwidth
b_t = min((omega_p1-omega_a1),(omega_a2-omega_p2));

%find the cutoff frequenceis
omega_c1=omega_p1-0.5*b_t;
omega_c2=omega_p2+0.5*b_t;

%find the delta values
delta_a = 10^-(0.05*Aa);
delta_p = ((10^(0.05*Ap))+1)/((10^(0.05*Ap))+1);
delta = min(delta_a,delta_p);

%find the new values of the stop band attenuation due to choosen delta values
Aa_new = -20*log10(delta);

%find the alpha
if (Aa_new<=21)
        alpha=0;
elseif (21<Aa_new && Aa_new<=50)
        alpha=0.5842*(Aa_new-21)^0.4+0.07886*(Aa_new-21);
else
    alpha=0.1102*(Aa_new-8.7); 
end

%find the D
if (Aa_new<=21)
    D=0.9222;
else
    D=(Aa_new-7.95)/14.36;
end

%find the order of the filter
N_cal=ceil(((omega_s*D)/b_t)+1);
N=N_cal+mod(N_cal+1,2);

%find the window function
%calculate bessel function : lo(alpha)
k=1;
term_k=1;
lo_a=1;
while (term_k>10^-6)
    term_k=((1/factorial(k))*((alpha/2)^k))^2;
    lo_a=lo_a+term_k;
    k=k+1;
end
wk=[];
for (n=0:(N-1)/2)
    beta=alpha*sqrt(1-(((2*n)/(N-1))^2));
    %calculate lo(beta)
    k=1;
    term_k=1;
    lo_b=1;
    while (term_k>10^-6)
        term_k=((1/factorial(k))*((beta/2)^k))^2;
        lo_b=lo_b+term_k;
        k=k+1;
    end
    
    wk(n+1)=(lo_b/lo_a);
end

%find the right side fourier impulse responce of the filter
hn_right=zeros(1,(N-1)/2);
for (n=0:(N-1)/2)
    if (n==0)
        hn_right(n+1)=(2*(omega_c2-omega_c1)/omega_s);
    else
        hn_right(n+1)=(sin(omega_c2*n*T)-sin(omega_c1*n*T))/(pi*n);
    end
end

%compute the impulse responce
hw_right=wk.*hn_right;
hw_left=fliplr(hw_right);
hn=[];
hn=[hw_left(1:(N-1)/2) hw_right];

n=[0:N-1];
figure
stem(n,hn,'filled')
title('Impulse response of the digital band pass filter')
ylabel('h(nT)')
xlabel('n')
n=[-(N-1)/2:(N-1)/2];
%compute the frequency responce of the digital filter
[freq_res,omega]=freqz(hn,1);
figure
plot(omega/T,20*log10(abs(freq_res)))
title('Magnitude response of the digital band pass filter')
grid on;
ylabel('Magnitude(dB)')
xlabel('Frequency(rad/s)')
axis([0 1700 -80 10])


figure
plot(omega/T,unwrap(angle(freq_res)))
title('Phase response of the digital band pass filter')
grid on;
ylabel('Phase Angle(rad)')
xlabel('Frequency(rad/s)')


figure
plot(omega/T,20*log10(abs(freq_res)))
title('Magnitude response of the digital band pass filter')
grid on;
ylabel('Magnitude(dB)')
xlabel('Frequency(rad/s)')
axis([omega_p1 omega_p2 -Ap/2 +Ap/2]);

%% 
%implimentation of the filter
%produce x(nT)
omega_1 = omega_a1*0.5;
omega_2 = (omega_p1+omega_p2)*0.5;
omega_3 = (omega_a2+omega_s*0.5)*0.5;
n=[1:800];
xn(n) = sin(omega_1*(n-1)*T)+sin(omega_2*(n-1)*T)+sin(omega_3*(n-1)*T);

%compute the output of the bandpass filter
Xw=fft(xn,1024);
Hw=fft(hn,1024);
Yw=(Xw.*Hw);
yn=ifft(Yw);


n=[1:800];
figure
subplot(3,1,1)
stem(n,xn(n))
title('Input to the digital band pass filter')
ylabel('x(nT)')
xlabel('n')
axis([0 400 -3 3])

%ideal bandpass filter output
x_ideal(n)= sin(omega_2*(n-1)*T);
figure
subplot(3,1,2)
stem(n,x_ideal(n))
title('Output of a ideal band pass filter')
ylabel('y(nT)')
xlabel('n')
axis([0 400 -1 1])

%desinged band pass filter output
figure
subplot(3,1,3)
stem(yn(n))
title('Output of the desinged band pass filter')
ylabel('y(nT)')
xlabel('n')
axis([0 400 -1 1])

%plot the frequency responces of input and output
%input frequency responce
w = (0:1/length(Xw):0.5-1/length(Xw))*omega_s;
figure
subplot(3,1,1)
plot(w,abs(Xw(1:length(Xw)/2)));
title('Frequency response of the input signal')
ylabel('Magnitude(dB)')
xlabel('Frequency(rad/s)')

%ideal frequency responce
subplot(3,1,2)
X_ideal=fft(x_ideal,1024);
plot(w,abs(X_ideal(1:length(Xw)/2)));
title('Frequency response of ideal band pass filter')
ylabel('Magnitude(dB)')
xlabel('Frequency(rad/s)')

%ideal frequency responce
subplot(3,1,3)
plot(w,abs(Yw(1:length(Xw)/2)));
title('Frequency response of designed band pass filter')
ylabel('Magnitude(dB)')
xlabel('Frequency(rad/s)')

fvtool(hn)