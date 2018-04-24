clear all;
close all;
clc;

%f=55.*(10^(6)) + (5.*(10^(6))+5.*(10^(6)))*rand(1,1)
%use for random input oof 55 to 65 MHz

% Add extra code for Bandpass filter, Low pass etc.
%Add AWGN and analyse results
f=( 60.*(10^(6)) );
%For the normal input wave of frequency 60MHz
%f=( 60.*(10^(6)) )./ ( 100.*(10^(6)) ) ;                     
% Sampled wave at 100 MHz
t1=0:(10^(-8)):0.05;

y=@(t)sin(2*pi*f*t);            
% Making the wave a sine function

d1=0.45*(10^(-6));              
% Time duration for which the wave has non zero amplitude
t=1.45*(10^(-6));               
% Time period (Total time for which a single wave exists) for the whole wave  
fs = 10^8;
x2=0:1/fs:t;

f=@(x)[y(x).*(0<=x & x<d1)+0.*(d1<=x & x<t)];       
% Division of the wave as per the desired time intervals
pfx2=f(x2);
plot(x2,pfx2);
title('Wave for a single Time Period');
figure;

intvl2=[-2*t 2*t];
z2=linspace(0,t-0.00000001);
fx2 = repmat(f(z2),1,diff(intvl2)/t);               
% Repeats the copies of the arrays for the single wave (generated) to be repeated and appear as a complete wave with multiple clock cycles

x3 = linspace(intvl2(1),intvl2(2),length(fx2));
plot(x3,fx2);
title('The Final Input wave (Continuous Time)');
figure;

stem(x3,fx2);
title('The Final Input wave (Discreet Time)');
% The final input wave or the sampled wave
figure;

%Gaussian White Noise is added to the input signal
N0 = awgn(pfx2,10,'measured');
plot(x2,N0);
title(' Plot of the wave with external AWGN ');
figure;



F = fft(fx2);
stem(1./x3,F);
title('Fourier Transform of the Sampled Signal');
figure;

stem(1./x3,abs(F));
title('Magnitude spectrum of Fourier Transform');
figure;

fc = 40.*(10^(6));

I0 = F .* cos((2.*(pi).*fc).*x3); 
plot(1./x3,I0);
title('Wave obtained on multiplying the fft with cos(fc*t), in the frequency domain');
figure;

Q0 = F .* sin((2.*(pi).*fc).*x3);
plot(1./x3,Q0);
title('Wave obtained on multiplying the fft with sin(fc*t), in the frequency domain');

%In the following steps we are designing an FIR Filter, through which we
%are going to pass the waves we obtain on getting the sine and cosine
%products of the Fourier Transformed Sampled Sine wave

d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.1,0.15,1,60);
% We are setting the Passband Frequency at 10 MHz and the Stopband
% frequency at 15 MHz and thus, the input arguments of 0.1 and 0.15 in
% lowpass object
designmethods(d);
f = design(d, 'ellip');  
%fvtool(f)   ;                   
Q = filter(f,Q0);
I = filter(f,I0);
stem(1./x3,I);
title('Plot of the fft* with sin wave obtained after filtering, in the frequency domain');
figure;

K = ((I).^2 + (Q).^2);
K1 = log(K);
stem(x3,K); 
title('Plot of the sum of sqaures of the waves obtained after filtering (I^2 + Q^2');
figure;
plot(x3,K1);
title('Plot of log of the sum of sqaures of the waves obtained after filtering log(I^2 + Q^2');

%% for the awgn noise:
clear all;
close all;
clc;

fn=(60.*(10^(6)));
fn1= 10.*(10^(6)):(10^(6)):120.*(10^(6))
n1 = length(fn1);
%t = 0:0.0000000001:0.000001
% What we have been doing so far which is wrong ( plot dekho zoom karke curved edges aani chahiye no sharp)
x1 = 0:0.00000001:1./fn;
%What we should have done ( see the pure sine wave hamare plot me pure nahi aa rahi bilkul
y = @(x)sin(2*pi*fn*x);
y1=y(x1);
plot(x1,y1);
%figure;

y1fft = abs(fft(y1));
%xf=-1*2*pi*fn:100000:2*pi*fn;
%plot (xf,y1fft);
%figure;

%plot(1./t,xdft);

d1=0.45*(10^(-6));              
% Time duration for which the wave has non zero amplitude
d2=1.45*(10^(-6)); 
f=@(x)[y(x).*(0<=x & x<d1)+0.*(d1<=x & x<d2)];
x2=0:0.00000001:d2;
pfx2=f(x2);
t=d2;


intvl2=[-2*t 2*t];
z=linspace(0,t-0.0000001);
fx2 = repmat(f(z),1,diff(intvl2)/t); 
x3 = linspace(intvl2(1),intvl2(2),length(fx2));


plot(x2,pfx2);

figure;
plot(x3,fx2);
figure;

stem(x3,fx2);

%Gaussian White Noise is added to the input signal
N0 = awgn(pfx2,10,'measured');
plot(x2,N0);
title(' Plot of the wave with external AWGN ');
figure;

F = fft(N0);
stem(1./x2,N0);
title('Fourier Transform of the Sampled Signal');
figure;

stem(x2,abs(N0));
title('Magnitude spectrum of Fourier Transform');
figure;

fc = 40.*(10^(6));

I0 = F .* cos((2.*(pi).*fc).*x2); 
plot(1./x2,I0);
title('Wave obtained on multiplying the fft with cos(fc*t), in the frequency domain');
figure;

Q0 = F .* sin((2.*(pi).*fc).*x2);
plot(1./x2,Q0);
title('Wave obtained on multiplying the fft with sin(fc*t), in the frequency domain');

%In the following steps we are designing an FIR Filter, through which we
%are going to pass the waves we obtain on getting the sine and cosine
%products of the Fourier Transformed Sampled Sine wave

d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.1,0.15,1,60);
% We are setting the Passband Frequency at 10 MHz and the Stopband
% frequency at 15 MHz and thus, the input arguments of 0.1 and 0.15 in
% lowpass object
designmethods(d);
f = design(d, 'ellip');  
%fvtool(f)   ;                   
Q = filter(f,Q0);
I = filter(f,I0);
stem(1./x2,I);
title('Plot of the fft* with sin wave obtained after filtering, in the frequency domain');
figure;

K = ((I).^2 + (Q).^2);
K1 = log(K);
stem(x2,K); 
title('Plot of the sum of sqaures of the waves obtained after filtering (I^2 + Q^2');
figure;
plot(x2,K1);
title('Plot of log of the sum of sqaures of the waves obtained after filtering log(I^2 + Q^2');


%%

f=55.*(10^(6)) + (5.*(10^(6))+5.*(10^(6)))*rand(1,1)
%use for random input oof 55 to 65 MHz



%f=( 60.*(10^(6)) );
%For the normal input wave of frequency 60MHz
%f=( 60.*(10^(6)) )./ ( 100.*(10^(6)) ) ;                     
% Sampled wave at 100 MHz
t1=0:(10^(-8)):0.05;

y=@(t)sin(2*pi*f*t);            
% Making the wave a sine function

d1=0.45*(10^(-6));              
% Time duration for which the wave has non zero amplitude
t=1.45*(10^(-6));               
% Time period (Total time for which a single wave exists) for the whole wave  
fs = 10^8;
x2=0:1/fs:t;

f=@(x)[y(x).*(0<=x & x<d1)+0.*(d1<=x & x<t)];       
% Division of the wave as per the desired time intervals
pfx2=f(x2);
plot(x2,pfx2);
title('Wave for a single Time Period');
figure;

intvl2=[-2*t 2*t];
z2=linspace(0,t-0.00000001);
fx2 = repmat(f(z2),1,diff(intvl2)/t);               
% Repeats the copies of the arrays for the single wave (generated) to be repeated and appear as a complete wave with multiple clock cycles

x3 = linspace(intvl2(1),intvl2(2),length(fx2));
plot(x3,fx2);
title('The Final Input wave (Continuous Time)');
figure;

stem(x3,fx2);
title('The Final Input wave (Discreet Time)');
% The final input wave or the sampled wave
figure;

%Gaussian White Noise is added to the input signal
N0 = awgn(pfx2,10,'measured');
plot(x2,N0);
title(' Plot of the wave with external AWGN ');
figure;



F = fft(fx2);
stem(1./x3,F);
title('Fourier Transform of the Sampled Signal');
figure;

stem(x3,abs(F));
title('Magnitude spectrum of Fourier Transform');
figure;

fc = 40.*(10^(6));

I0 = F .* cos((2.*(pi).*fc).*x3); 
plot(1./x3,I0);
title('Wave obtained on multiplying the fft with cos(fc*t), in the frequency domain');
figure;

Q0 = F .* sin((2.*(pi).*fc).*x3);
plot(1./x3,Q0);
title('Wave obtained on multiplying the fft with sin(fc*t), in the frequency domain');

%In the following steps we are designing an FIR Filter, through which we
%are going to pass the waves we obtain on getting the sine and cosine
%products of the Fourier Transformed Sampled Sine wave
0
d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.1,0.15,1,60);
% We are setting the Passband Frequency at 10 MHz and the Stopband
% frequency at 15 MHz and thus, the input arguments of 0.1 and 0.15 in
% lowpass object
designmethods(d);
f = design(d, 'ellip');  
%fvtool(f)   ;                   
Q = filter(f,Q0);
I = filter(f,I0);
stem(1./x3,I);
title('Plot of the fft* with sin wave obtained after filtering, in the frequency domain');
figure;

K = ((I).^2 + (Q).^2);
K1 = log(K);
stem(x3,K); 
title('Plot of the sum of sqaures of the waves obtained after filtering (I^2 + Q^2');
figure;
plot(x3,K1);
title('Plot of log of the sum of sqaures of the waves obtained after filtering log(I^2 + Q^2');

