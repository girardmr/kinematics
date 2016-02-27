close all
clear all
clc

%% Inputs

a = 0.8; %link 2
b = 1.23; %link 3
c = 1.55; %link 4
d = 2.4; %link 1	
theta_2_input = 45;
omega_2_input = -15;
alpha_2_input = -10;
configuration = 0; %0 for open, 1 for crossed

%% Grashof Condition

links = [a b c d];
S = min(links);
L = max(links);
sum_SL = S+L;
sum_PQ = sum(links)-sum_SL;

if sum_SL < sum_PQ
    disp('Grashof Linkage');
elseif sum_SL > sum_PQ 
    disp('Non-Grashof Linkage');
    arg1 = ((a^2+d^2-b^2-c^2)/(2*a*d))+((b*c)/(a*d));
    arg2 = ((a^2+d^2-b^2-c^2)/(2*a*d))-((b*c)/(a*d));
    if arg1 < 1 && arg1 > -1;
        theta_toggle = acos(arg1);
    else
        theta_toggle = acos(arg2);
    end
    theta_toggle_d = theta_toggle*(180/pi);
    theta_2_input = -theta_toggle_d;
    theta2_nonGrashof = [-theta_toggle_d:theta_toggle_d];
elseif sum_SL == sum_PQ
    disp('Special Case Grashof');
end

%% 
delta_theta_2 = 1; %increment of crank angle in degrees
delta_theta_2 = delta_theta_2*(pi/180);
theta_2_final = theta_2_input-1;
total_num_positions = theta_2_input+(360-theta_2_input);   %number of angles in degrees
%total_num_positions_r = total_num_positions*(pi/180);  %radians
if sum_SL > sum_PQ
    total_num_positions = length(theta2_nonGrashof);
end

delta = 0;                          %angle between link 3 and coupler point?
delta_3=delta;					    %Coupler point angle in degrees. 
delta_3=delta_3*pi/180;				%Coupler point angle converted to radians.
p1=b;
p=p1;							    %Distance from point A to the coupler point...equal to length of link 3 for fourbar mechanism?

%%

theta1 = zeros(1,total_num_positions);  %degrees
theta2 = zeros(1,total_num_positions);
theta2(1) = theta_2_input;
if sum_SL > sum_PQ
    theta2 = theta2_nonGrashof;
end
theta3 = zeros(1,total_num_positions);
theta4 = zeros(1,total_num_positions);

omega1 = zeros(1,total_num_positions);
%constant omega
omega2 = zeros(1,total_num_positions);
omega2(1) = omega_2_input;
omega2 = ones(1,total_num_positions)*omega_2_input;
omega3 = zeros(1,total_num_positions);
omega4 = zeros(1,total_num_positions);

alpha1 = zeros(1,total_num_positions);
alpha2 = zeros(1,total_num_positions);
alpha2(1) = alpha_2_input;
alpha3 = zeros(1,total_num_positions);
alpha4 = zeros(1,total_num_positions);

transmission_angle=zeros(1,total_num_positions);

%% Checks if linkage is beyond limits or lengths are not permissible

K1=d/a;									%Eq. (4.8a)
K4=d/b;
K5=(c^2-d^2-a^2-b^2)/(2*a*b);           %Eq. (4.11a)
								    	%Make sure to have parentheses in the denominator.
for jj = 2:length(theta2)
    theta2(jj) = theta2(jj-1)+1;
end
theta2 = theta2*(pi/180);

for ii = 1:length(theta2)
    
    D=cos(theta2(ii))-K1+K4*cos(theta2(ii))+K5;	%Preceding Eq. (4.13)
    E=-2*sin(theta2(ii));
    F=K1+(K4-1)*cos(theta2(ii))+K5;
%     if (E^2-4.*D*F<0) 									% No solution for theta_3 in Eq. (4.13).
%         return;
%     end
    
    %% Position Analysis
    
    if configuration == 0; %open
        theta3(ii) = 2*atan((-E-sqrt(E^2-4.*D*F))/(2*D)); 
    elseif configuration == 1; %crossed
        theta3(ii) = 2*atan((-E+sqrt(E^2-4.*D*F))/(2*D)); 
    end
    K2 = d/c;
    K3 = (a^2-b^2+c^2+d^2)/(2*a*c);
    A = cos(theta2(ii))-K1-K2*cos(theta2(ii))+K3;
    B = -2*sin(theta2(ii));
    C = K1-(K2+1)*cos(theta2(ii))+K3;
    theta4(ii) = 2*atan((-B+sqrt(B^2-4.*A*C))/(2*A));
    
    transmission_angle(ii) = abs(theta3(ii)-theta4(ii));
    if (transmission_angle(ii)>pi)
        transmission_angle(ii)=2*pi-transmission_angle(ii);
    end
    if (transmission_angle(ii)>pi/2)
        transmission_angle(ii)=pi-transmission_angle(ii);
    end
    transmission_angle(ii)=transmission_angle(ii)*180/pi;
    if (ii == 1)
        transmission_max = transmission_angle(1);
        transmission_min = transmission_angle(1);  
    end
    transmission_max = max(transmission_max, transmission_angle(ii));
    transmission_min = min(transmission_min, transmission_angle(ii));
    
    %% Velocity Analysis
    
    omega3(ii)=a*omega2(ii)*(sin(theta4(ii)-theta2(ii)))/(b*(sin(theta3(ii)...
        -theta4(ii))));														%Eq. (6.18a)
    omega4(ii)=a*omega2(ii)*(sin(theta2(ii)-theta3(ii)))/(c*(sin(theta4(ii)...
        -theta3(ii))));	
    VA(ii) = a*omega2(ii)*(-sin(theta2(ii))+j*cos(theta2(ii)));
    VB(ii) = c*omega4(ii)*(-sin(theta4(ii))+j*cos(theta4(ii)));
    %VA = velocity point A w/r/t origin..aka point that connects links 2 and 3
    %VB = velocity point B w/r/t origin..aka point that connects links 3 and 4...aka
    %coupler point
    %VA+VBA-VB = 0
    VBA(ii) = VB(ii)-VA(ii);
    
    %% Acceleration Analysis
    
    A = c*sin(theta4(ii));
    B = b*sin(theta3(ii));
    C = a*alpha2(ii)*sin(theta2(ii))+a*(omega2(ii)^2)*cos(theta2(ii))...
        +b*omega3(ii)^2*cos(theta3(ii))-c*omega4(ii)^2*cos(theta4(ii));
    D = c*cos(theta4(ii));
    E = b*cos(theta3(ii));
    F = a*alpha2(ii)*cos(theta2(ii))-a*omega2(ii)^2*sin(theta2(ii))...
        -b*omega3(ii)^2*sin(theta3(ii))+c*omega4(ii)^2*sin(theta4(ii));
    alpha3(ii) = (C*D-A*F)/(A*E-B*D);
    alpha4(ii) = (C*E-B*F)/(A*E-B*D);
    %acceleration of point A...aka point that connects links 2 and 3
    AA = a*alpha2(ii)*(-sin(theta2(ii))+j*cos(theta2(ii)))...
        -a*omega2(ii)^2*(cos(theta2(ii))+j*sin(theta2(ii)));
    %acceleration of point B...aka point that connects links 3 and 4...aka
        %coupler point
    AB = c*alpha4(ii)*(-sin(theta4(ii))+j*cos(theta4(ii)))...
        -c*omega4(ii)^2*(cos(theta4(ii))+j*sin(theta4(ii)));
    %A_A+A_BA-A_B = 0
    ABA = AB-AA;
    
end
    
    
%% Animation 
    
XC = [];YC=[];
I13X=[];I13Y=[];I24X=[];I24Y=[];
    
for k = 1:ii-1
    
% Calculates the coordinates of joints and the coupler point based on vector equation with nomenclature in Figure 4-7.

xo2=0;                              %Coordinates of point O2
yo2=0;
xo4=xo2+d*cos(theta1(1));			%Coordinates of point O4
yo4=yo2+d*sin(theta1(1));	
	%theta(1,k) is 0 in all HW problems.
xa=xo2+a*cos(theta2(k));			%Coordinates of point A
ya=yo2+a*sin(theta2(k));
xb=xa+b*cos(theta3(k));             %Coordinates of point B
yb=ya+b*sin(theta3(k));	
	%theta(3,k) is calculated in the position analysis section.
%Alternative codes
%xb=xo4+c*cos(theta(4,k));			%Coordinates of point B
%yb=yo4+c*sin(theta(4,k));	
xc=xa+p*cos(theta3(k)+delta_3);	%Coordinates of coupler point 
yc=ya+p*sin(theta3(k)+delta_3);

% Creates the lines representing binary links and sides of the coupler link connected to the coupler point.
x_link2=[xo2 xa];			%Line of link 2
y_link2=[yo2 ya];
x_link3=[xa xb];			%Line of link 3
y_link3=[ya yb];
x_link4=[xb xo4];			%Line of link 4
y_link4=[yb yo4];
x_link1=[xo4 xo2];      	%Line of link 1
y_link1=[yo4 yo2];
x_c1=[xa xc];				%Line of coupler_link side 1
y_c1=[ya yc];
x_c2=[xb xc];				%Line of coupler_link side 2
y_c2=[yb yc];
XC=[XC,xc];
YC=[YC,yc];

% Calculate the coordinate of instant centers
I13x=[xb*(yb-yo4)*(xa-xo2)-xa*(ya-yo2)*(xb-xo4)+(ya-yb)*(xb-xo4)*(xa-xo2)]/[(yb-yo4)*(xa-xo2)-(ya-yo2)*(xb-xo4)];
I13y=[yb*(xb-xo4)*(ya-yo2)-ya*(xa-xo2)*(yb-yo4)+(xa-xb)*(yb-yo4)*(ya-yo2)]/[(xb-xo4)*(ya-yo2)-(xa-xo2)*(yb-yo4)];
I24x=[xb*(yb-ya)*(xo4-xo2)-xo4*(yo4-yo2)*(xb-xa)+(yo4-yb)*(xb-xa)*(xo4-xo2)]/[(yb-ya)*(xo4-xo2)-(yo4-yo2)*(xb-xa)];
I24y=[yb*(xb-xa)*(yo4-yo2)-yo4*(xo4-xo2)*(yb-ya)+(xo4-xb)*(yb-ya)*(yo4-yo2)]/[(xb-xa)*(yo4-yo2)-(xo4-xo2)*(yb-ya)];
I13X=[I13X,I13x];
I13Y=[I13Y,I13y];
I24X=[I24X,I24x];
I24Y=[I24Y,I24y];

plot(xo2,yo2,'-.ko',x_link2,y_link2,xa,ya,'-.ko',x_link3,y_link3,xb,yb,'-.ko',x_link4,y_link4,xo4,yo4,'-.ko',x_link1,y_link1,x_c1,y_c1,x_c2,y_c2,...
    I13x,I13y,'o', I24x,I24y,'o', xc,yc,'o', XC,YC);
hold on;
xmin=(xo2-a);
xmax=(xo4+c);
ymin=max(yo2-a, yo4-c)-p1*cos(delta*pi/180);
ymax=max(yo2+a, yo4+c)+p1*cos(delta*pi/180);
axis([xmin  xmax ymin ymax]*1.5)	
%axis equal 
%axis([xo2-a-2 xo4+c+2 yo2-a-2 yo2-a-2+(xo4+c-xo2+a+4)]*2);
plot([xmin  xmax ]*1.4,[0 0], '-.k', [0 0],[ymin ymax]*1.2,'-.k')
hold off;

% %plot coupler curve
% plot(xo2,yo2,'-.ko',x_link2,y_link2,xa,ya,'-.ko',x_link3,y_link3,xb,yb,'-.ko',x_link4,y_link4,xo4,yo4,'-.ko',x_link1,y_link1,x_c1,y_c1,x_c2,y_c2, xc,yc,'o', XC,YC)
% axis([xo2-a-2 xo4+c+2 yo2-a-2 yo2-a-2+(xo4+c-xo2+a+4)]*2)	
% %plot fixed centrodes
% plot(xo2,yo2,'-.ko',x_link2,y_link2,xa,ya,'-.ko',x_link3,y_link3,xb,yb,'-.ko',x_link4,y_link4,xo4,yo4,'-.ko',x_link1,y_link1, I13x,I13y,'o', I13X,I13Y)
% axis([xo2-a-2 xo4+c+2 yo2-a-2 yo2-a-2+(xo4+c-xo2+a+4)]*2)

drawnow;
end

theta2 = theta2.*(180/pi);
theta3 = theta3.*(180/pi);
theta4 = theta4.*(180/pi);




