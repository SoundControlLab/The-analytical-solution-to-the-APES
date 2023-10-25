primary_audio = 'speech2.wav';
ambient_audio = 'W2.wav';

[p_truee, fs] = audioread(primary_audio);
[a_truee, fs] = audioread(ambient_audio);
framelen = 2^12;
framenum = 40;
a_true0 = a_truee(1:framelen*framenum,1);
p_true0 = p_truee(1:framelen*framenum,1);
k = 2;
PPR_range=[0.1:0.05:0.9];

for ii=1:framenum
    a_true0_frame=a_true0(1+framelen*(ii-1):framelen*ii);
    fre_a_true0=fft(a_true0_frame);
    phase_a_true0=angle(fre_a_true0);
    phase_a_true1=zeros(1,length(fre_a_true0));
    phase_a_true1(length(fre_a_true0)/2+1:length(fre_a_true0))=rand(length(fre_a_true0)/2,1)*2*pi;
    phase_a_true1(2:length(fre_a_true0)/2)=-fliplr(phase_a_true1(length(fre_a_true0)/2+2:length(fre_a_true0)));
    phase_a_true1(1)=0;
    fre_a_true1=abs(fre_a_true0).*exp(1i*phase_a_true1.');
    a_true1_frame=ifft(fre_a_true1);
    a_true1(1+framelen*(ii-1):framelen*ii)=a_true1_frame;
end
a_true1= a_true1.';
a_true0=a_truee(1:framelen*framenum,1);
p_true0=p_truee(1:framelen*framenum,1);
p_true1=k*p_true0;

pw_pri_l0 = p_true0'*p_true0;
pw_pri_r0 = p_true1'*p_true1;
pw_amb_l0 = a_true0'*a_true0;
pw_amb_r0 = a_true1'*a_true1;

for i=1:length(PPR_range)
    PPR = PPR_range(i)

    p_true = [p_true0 p_true1]*sqrt(PPR/(1-PPR))/(sqrt(pw_pri_l0+pw_pri_r0));
    a_true = [a_true0 a_true1]/sqrt(pw_amb_l0+pw_amb_r0);
    x=p_true+a_true;

    pw_pri_l0 = p_true0'*p_true0;
    pw_pri_r0 = p_true1'*p_true1;
    pw_amb_l0 = a_true0'*a_true0;
    pw_amb_r0 = a_true1'*a_true1;

    p_true = [p_true0 p_true1]*sqrt(PPR/(1-PPR))/(sqrt(pw_pri_l0+pw_pri_r0));
    a_true = [a_true0 a_true1]/sqrt(pw_amb_l0+pw_amb_r0);
    x=p_true+a_true;
    [err_AS_P,err_AS_A]=AS(x,k,p_true,a_true,framenum,framelen,PPR);
end