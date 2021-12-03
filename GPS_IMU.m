clc
clear all
close all

%% DATA Reading
fileID = fopen('data_1.txt','r');
DATA = fscanf(fileID,'%f',[11 Inf]); %
N = size(DATA);
Nsamples = N(2); %length of DATA

%% INITIALIZING
g = 9.8;
unit_transform_acc = 2048;
unit_transform_gyro = (pi/(180*16.4));
unit_transform_mag=0.6;
n=0;
ch=0;
N_Q = 1;
N_R = 100;
N_P = 1;
v(1)=0;
fr=0;
df=0;
Bias_AccX_sum=0;
Bias_AccY_sum=0;
Bias_AccZ_sum=0;

%% RAW DATA TO SI DATA
for k = 1:Nsamples
    if DATA(1,k)==1
        DATA_SI(1,k)=1;
        DATA_SI(2,k)=DATA(10,k);
        DATA_SI(3,k)=DATA(11,k);
        DATA_SI(11,k)=1000*(ch);
        ch=ch+1;
    else
        %Acc LSB -> m/s^2
        DATA_SI(1,k)=0;
        DATA_SI(2,k)= (g/unit_transform_acc)*DATA(2,k);
        DATA_SI(3,k)= (g/unit_transform_acc)*DATA(3,k);
        DATA_SI(4,k)= (g/unit_transform_acc)*DATA(4,k);
        %Gyro LSB -> deg/s -> rad/s
        DATA_SI(5,k)= (unit_transform_gyro)*DATA(5,k);
        DATA_SI(6,k)= (unit_transform_gyro)*DATA(6,k);
        DATA_SI(7,k)= (unit_transform_gyro)*DATA(7,k);
        %Mag LSB -> uT
        DATA_SI(8,k)= unit_transform_mag*DATA(8,k);
        DATA_SI(9,k)= unit_transform_mag*DATA(9,k);
        DATA_SI(10,k)= unit_transform_mag*DATA(10,k);
        %Time ms -> s
        DATA_SI(11,k)=DATA(11,k)+(1000*(ch-1));
        let(k-ch)=DATA(11,k)+(1000*(ch-1));
    end
end

%% Accel Compensation
for k = 1:40
    if DATA_SI(1,k)==0
        Bias_AccX_sum=Bias_AccX_sum+DATA_SI(2,k);
        Bias_AccY_sum=Bias_AccY_sum+DATA_SI(3,k);
        df=df+1;
    end
end

Bias_Acc = [Bias_AccX_sum/df Bias_AccY_sum/df];

for k = 1:Nsamples
    if DATA_SI(1,k)==0
        DATA_SI(2,k)=DATA_SI(2,k)-Bias_Acc(1);
        DATA_SI(3,k)=DATA_SI(3,k)-Bias_Acc(2);
    end
end



%% Magnetometer Compensation

ch=0;
for k = 1:Nsamples
    
    if DATA_SI(1,k)==0
        Y(k-ch,:) = [DATA_SI(8,k)^2+DATA_SI(9,k)^2+DATA_SI(10,k)^2];
        X(k-ch,:) = [DATA_SI(8,k) DATA_SI(9,k) DATA_SI(10,k) 1];
        n=n+1;
    end
    
    if DATA_SI(1,k)==1
        ch=ch+1;
    end
end

N_X = [size(X)];
Xsamples = N_X(2);
Bias_Mag = 0.5*((X'*X)\eye(Xsamples))*X'*Y;

for k = 1:Nsamples
    if DATA_SI(1,k)==0
        DATA_SI(8,k) = DATA_SI(8,k) - Bias_Mag(1);
        DATA_SI(9,k) = DATA_SI(9,k) - Bias_Mag(2);
        DATA_SI(10,k) = DATA_SI(10,k) - Bias_Mag(3);
    end
end

%% Set Reference Magnetic vector (NORMALIZATION)
[~,ref_mag]=max(DATA_SI(8,:));
M=sqrt(DATA_SI(8,ref_mag)^2+DATA_SI(9,ref_mag)^2+DATA_SI(10,ref_mag)^2);
B=[DATA_SI(8,ref_mag)/M DATA_SI(9,ref_mag)/M DATA_SI(10,ref_mag)/M];

%% EKF Algorithm
ch=0;
for k = 1:Nsamples
    if DATA(1,k)==1
        
        lat(k)=DATA_SI(2,k);
        long(k)=DATA_SI(3,k);
        
        origin=[lat(k),long(k),0];
        
        [x(k),y(k),~] = latlon2local(DATA_SI(2,k),DATA_SI(3,k),0,origin);
        
        lat_gps(ch+1)=lat(k);
        long_gps(ch+1)=long(k);
        
        ch=ch+1;
    else
        %Assignment
        ax=DATA_SI(2,k);
        ay=DATA_SI(3,k);
        az=DATA_SI(4,k);
        p=DATA_SI(5,k);
        q=DATA_SI(6,k);
        r=DATA_SI(7,k);
        mx=DATA_SI(8,k);
        my=DATA_SI(9,k);
        mz=DATA_SI(10,k);
        dt=(DATA_SI(11,k+2-ch)-DATA_SI(11,k+1-ch))/1000;
        %Normalization
        G=sqrt(ax^2+ay^2+az^2);
        M=sqrt(mx^2+my^2+mz^2);
        ax=ax/G; ay=ay/G; az=az/G; mx=mx/M; my=my/M; mz=mz/M;
        %MAIN EKF FUNCTION
        [q0, q1, q2, q3] = EKF(p, q, r, B, mx, my, mz, ax, ay, az, dt, N_Q, N_R, N_P);
        %Conversion to Euler angle
        phi   =  atan2( 2*(q2*q3 + q0*q1), 1 - 2*(q1^2 + q2^2) );
        theta = -asin(  2*(q1*q3 - q0*q2) );
        psi   =  atan2( 2*(q1*q2 + q0*q3), 1 - 2*(q2^2 + q3^2) );
        EulerSaved(k-ch,:) = [ phi theta psi ];
        
        v=(-ax)*dt+v;
        dx=0.5*(-ax)*dt^2+v*dt;
        x(k)=x(k-1)+dx*cos(psi-0.25*pi);
        y(k)=y(k-1)+dx*sin(psi-0.25*pi);
        
        [lat(k),long(k),~] = local2latlon(x(k),y(k),0,origin);
        lat_imu(k-ch)=lat(k);
        long_imu(k-ch)=long(k);
        
    end
    
end


%% Plot
PhiSaved   = EulerSaved(:, 1) * 180/pi;
ThetaSaved = EulerSaved(:, 2) * 180/pi;
PsiSaved   = EulerSaved(:, 3) * 180/pi;
x = [0 0];
y = [-1000,1000];

figure()
P1=plot(let, PhiSaved, 'r');
hold on
P2=plot(let, ThetaSaved, 'b');
P3=plot(let, PsiSaved, 'g');
refline([0 0])
title('Euler Angle (degree)')
Timeline_1 = line('XData',x,'YData',y);
TimeValue_1= xlabel('');
legend([P1 P2 P3],{'Phi', 'Theta', 'Psi'},'Location','northwest','AutoUpdate','off');
axis([0 DATA_SI(11,Nsamples) -300 300])

figure()
gx = geoaxes;
PLOT=geoplot(gx,lat,long,'-.r*');
geobasemap(gx,'satellite');

figure()
gx = geoaxes;
hold on
PLOT_LINE=geoplot(gx,lat_imu,long_imu,'r*');
hold on
PLOT_DOT=geoplot(gx,lat_gps,long_gps,'g*');
hold on
geobasemap(gx,'satellite');

%% Dynamic Plot
% gx = geoaxes;
% PLOT_LINE=geoplot(gx,lat_imu(1),long_imu(1),'r*');
% hold on
% PLOT_DOT=geoplot(gx,lat_gps(1),long_gps(1),'g*');
% geobasemap(gx,'satellite');
% Realtime = title(' ');
% 
% ch=1;
% pause(5);
% for k=2:Nsamples
%     if DATA_SI(1,k)==1
%         ch=ch+1;
%         set(PLOT_LINE,'XData',lat_imu(1:k-ch));
%         set(PLOT_LINE,'YData',long_imu(1:k-ch));
%         set(PLOT_DOT,'XData',lat_gps(1:ch));
%         set(PLOT_DOT,'YData',long_gps(1:ch));
%     else
%         set(PLOT_DOT,'XData',lat_gps(1:ch));
%         set(PLOT_DOT,'YData',long_gps(1:ch));
%         set(PLOT_LINE,'XData',lat_imu(1:k-ch));
%         set(PLOT_LINE,'YData',long_imu(1:k-ch));
%     end
%     drawnow
% end