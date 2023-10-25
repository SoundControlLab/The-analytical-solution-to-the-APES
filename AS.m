
function [err_AS_P,err_AS_A]=AS(x,k,p_true,a_true,framenum,framelen,PPR)
%% analytic solution

for i=1:framenum
    x0=x(1+framelen*(i-1):framelen*i,1);
    x1=x(1+framelen*(i-1):framelen*i,2);
    x0_frame=fft(x0,framelen);
    x1_frame=fft(x1,framelen);
    x0_re=x0_frame(framelen/2+1:framelen);
    x1_re=x1_frame(framelen/2+1:framelen);
    thta=angle(x1_re-k*x0_re);
    rotation_theta=exp(1i*(-pi/2-thta));
    p_=k.*x0_re.*rotation_theta;
    q_=x1_re.*rotation_theta;
    x_p=real(p_);
    y_p=imag(p_);
    y_q=imag(q_);
    if k>1
        d_=sqrt(x_p.^2+((k*k*y_q-y_p)/(k*k-1)).^2);
        r_=abs(x1_re-k*x0_re)*k/(k*k-1);
        p_min=abs(d_-r_);
        theta_p = exp(1i*angle(x_p+1i*((k*k*y_q-y_p)/(k*k-1))));
        for f=1:framelen/2
            if d_(f)>r_(f)
                PP1_est(f) = p_min(f)*theta_p(f)/rotation_theta(f);
            else PP1_est(f) = p_min(f)*-1*theta_p(f)/rotation_theta(f);
            end
        end
    else
        p_min=abs(y_p+y_q)/2;
        for f=1:framelen/2
            if (y_p(f)+y_q(f))>0
                theta_p = exp(1i*pi/2);
                PP1_est(f) = p_min(f)*theta_p/rotation_theta(f);
            else
                theta_p = exp(1i*3*pi/2);
                PP1_est(f) = p_min(f)*theta_p/rotation_theta(f);
            end
        end
    end
    p1_fre_AS=[0 conj(fliplr(PP1_est(2:framelen/2))) PP1_est];
    p0_fre_AS=p1_fre_AS/k;
    a1_fre_AS=x1_frame.'-p1_fre_AS;
    a0_fre_AS=x0_frame.'-p0_fre_AS;
    A1_est_AS=ifft(a1_fre_AS,framelen);
    A0_est_AS=ifft(a0_fre_AS,framelen);
    P1_est_AS=ifft(p1_fre_AS,framelen);
    p_est_AS(1+framelen*(i-1):framelen*i,:)=[1/k*P1_est_AS.' P1_est_AS.'];
    a_est_AS(1+framelen*(i-1):framelen*i,:)=[A0_est_AS.' A1_est_AS.'];
end

err_AS_pri1=(p_true(1:framenum*framelen,1)-real(p_est_AS(:,1))).'*(p_true(1:framenum*framelen,1)-real(p_est_AS(:,1)))/(p_true(1:framenum*framelen,1).'*p_true(1:framenum*framelen,1));
err_AS_pri2=(p_true(1:framenum*framelen,2)-real(p_est_AS(:,2))).'*(p_true(1:framenum*framelen,2)-real(p_est_AS(:,2)))/(p_true(1:framenum*framelen,2).'*p_true(1:framenum*framelen,2));
err_AS_amb1=(a_true(1:framenum*framelen,1)-real(a_est_AS(:,1))).'*(a_true(1:framenum*framelen,1)-real(a_est_AS(:,1)))/(a_true(1:framenum*framelen,1).'*a_true(1:framenum*framelen,1));
err_AS_amb2=(a_true(1:framenum*framelen,2)-real(a_est_AS(:,2))).'*(a_true(1:framenum*framelen,2)-real(a_est_AS(:,2)))/(a_true(1:framenum*framelen,2).'*a_true(1:framenum*framelen,2));

err_AS_P=10*log10(0.5*err_AS_pri1+0.5*err_AS_pri2);
err_AS_A=10*log10(0.5*err_AS_amb1+0.5*err_AS_amb2);

end