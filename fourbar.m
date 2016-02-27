close all
clear all
clc

%inputs

a=0.8; %link 2
b=1.23; %link 3
c=1.55; %link 4
d=2.4; %link 1	
theta_2_input = 45;
omega_2_input = -15;
alpha_2_input = -10;
    
% ---------------------------------------------------------------------------
% Grashof Condition
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
    %theta_2_input = -theta_toggle_d;
elseif sum_SL == sum_PQ
    disp('Special Case Grashof');
end

%-------------------------------------------------------------------------------
   
delta_theta_2=2;									%Increment of crank angle in degrees. 
delta_theta_2=delta_theta_2*pi/180;			%Increment of crank angle converted to radians.
n=2*pi/delta_theta_2;							%Total number of animation intervals in a 
														% complete revolution.
	%Decreases delta_theta to get a smoother animation. 
   %Increases delta_theta to 360 degrees (so that n=1) for a quick calculation at 
   %	one position with no animation. This will be useful for checking HW results.    
% ----------------------------------------------------------------------------
theta=zeros(4,n+1);					%Initializes matrices.  
omega=zeros(4,n+1);					 
alpha=zeros(4,n+1);					
transmission_angle=zeros(1,n+1);
theta1=0;
theta(1,1)=theta1;
% ----------------------------------------------------------------------------
%The first index of the matrix, from 1 to 4, denotes the link number, and the 
%	corresponding elements are theta_1, omega_2, etc.
%The second index of the matrix, from 1 to n+1, denotes the iteration number in a 
%	motion simulation. 
%For example, theta(2,1) denotes the initial joint angle of theta_2.
%The maximum number of iteration is n+1 so that the link returns to the original 
%	position. 
% ----------------------------------------------------------------------------
delta = 0;                          %angle between link 3 and coupler point?
delta_3=delta;					    %Coupler point angle in degrees. 
delta_3=delta_3*pi/180;				%Coupler point angle converted to radians.
p1=b;
p=p1;							    %Distance from point A to the coupler point...equal to length of link 3 for fourbar mechanism?
% ----------------------------------------------------------------------------
theta(2,1)=theta_2_input;			%Initial joint angle theta_2 in degrees 
	%Modify this angle to get a new inital linkage position
theta(2,1)=theta(2,1)*pi/180; 	%Initial joint angle theta_2 converted to radians
omega(2,1)= omega_2_input;		       		%Initial angular velocity in rad/s
%Note that data given in HW is already in rad/s. No need for conversion.
alpha(2,1)= alpha_2_input;                %Initial angular acceleration in rad/s^2
% ----------------------------------------------------------------------------
wflag=1; 					%Initializes the flag of linkage existance. 
i=1;							%Initializes the counter of while loop for motion simulation.
% ----------------------------------------------------------------------------
%Checks if the linkage is beyond its limits or lengths are not permissible for a 4bar     
% ----------------------------------------------------------------------------
K1=d/a;									%Eq. (4.8a)
K4=d/b;
K5=(c^2-d^2-a^2-b^2)/(2*a*b); 	%Eq. (4.11a)
											%Make sure to have parentheses in the denominator.
D=cos(theta(2,i))-K1+K4*cos(theta(2,i))+K5;	%Preceding Eq. (4.13)
E=-2*sin(theta(2,i));
F=K1+(K4-1)*cos(theta(2,i))+K5;
if (E^2-4.*D*F<0) 									% No solution for theta_3 in Eq. (4.13).
   wflag=0;
   	% If wflag=0, The linkage is beyond its limits, or the lengths 
      %		are not permissible for a 4bar. 
      % For example, if a^2+d^2+2a*d*cos(theta_2) > (b+c)^2,
      %		the fourbar_export linkage can not be constructed.
      % At the motion limit, the fourbar_export becomes a triagle, as shown in 
      %		Figure 4.16.  If the crank is over the limit, the law of cosine will not 
      %		be valid, and the fourbar_export can not be connected.
end 
% ----------------------------------------------------------------------------
%Perform position, velocity, and acceleration analyses during motion (while loop). 
%If wflag=0, there is no need to set the linkage in motion in the following while loop.
% ----------------------------------------------------------------------------
while wflag==1     
   	%Note that "==" is a relational operator to return logical true or false.  
 		%If linkage is permissible (wflag==1), the while loop proceeds.       
% ----------------------------------------------------------------------------
%Position Analysis
% ----------------------------------------------------------------------------
theta(3,i)=2*atan((-E+sqrt(E^2-4.*D*F))/(2*D)); 				%Eq. (4.13)
	% Open configuration 
	% For crossed configuration, change the sign before sqrt. 
    % - = open; + = closed
K2=d/c;
K3=(a^2-b^2+c^2+d^2)/(2*a*c);
A=cos(theta(2,i))-K1-K2*cos(theta(2,i))+K3;
B=-2*sin(theta(2,i));
C=K1-(K2+1)*cos(theta(2,i))+K3;
theta(4,i)=2*atan((-B+sqrt(B^2-4.*A*C))/(2*A));					%Eq. (4.10)  
% ----------------------------------------------------------------------------
transmission_angle(i)=abs(theta(3,i)-theta(4,i));

if (transmission_angle(i)>pi)
transmission_angle(i)=2*pi-transmission_angle(i);
end
if (transmission_angle(i)>pi/2)
transmission_angle(i)=pi-transmission_angle(i);
end
transmission_angle(i)=transmission_angle(i)*180/pi;
%
if (i == 1)
    transmission_max = transmission_angle(1);
    transmission_min = transmission_angle(1);  
end
transmission_max = max(transmission_max, transmission_angle(i));
transmission_min = min(transmission_min, transmission_angle(i));

% set(handles.max_button,'String',transmission_max)
% set(handles.min_button,'String',transmission_min)

% ----------------------------------------------------------------------------
%Velocity Analysis
% ----------------------------------------------------------------------------
omega(3,i)=a*omega(2,i)*(sin(theta(4,i)-theta(2,i)))/(b*(sin(theta(3,i)...
  -theta(4,i))));														%Eq. (6.18a)
omega(4,i)=a*omega(2,i)*(sin(theta(2,i)-theta(3,i)))/(c*(sin(theta(4,i)...
  -theta(3,i))));														%Eq. (6.18b)
																% "..." means continuation
VA =a*omega(2,i)*(-sin(theta(2,i))+j*cos(theta(2,i)));		%Eq. (6.19a)
VB =c*omega(4,i)*(-sin(theta(4,i))+j*cos(theta(4,i)));		%Eq. (6.19c)
	%	Note that VA and VB are complex numbers.
%VA = velocity point A w/r/t origin..aka point that connects links 2 and 3
%VB = velocity point B w/r/t origin..aka point that connects links 3 and 4...aka
%coupler point
%VA+VBA-VB = 0
VBA = VB-VA;
% ----------------------------------------------------------------------------
%Acceleration Analysis
% ----------------------------------------------------------------------------
A=c*sin(theta(4,i));														%Eq. (7.12c)
B=b*sin(theta(3,i));												
C=a*alpha(2,i)*sin(theta(2,i))+a*(omega(2,i)^2)*cos(theta(2,i))...
	+b*omega(3,i)^2*cos(theta(3,i))-c*omega(4,i)^2*cos(theta(4,i));
D=c*cos(theta(4,i));
E=b*cos(theta(3,i));
F=a*alpha(2,i)*cos(theta(2,i))-a*omega(2,i)^2*sin(theta(2,i))...
	-b*omega(3,i)^2*sin(theta(3,i))+c*omega(4,i)^2*sin(theta(4,i));
alpha(3,i)= (C*D-A*F)/(A*E-B*D);										%Eq. (7.12a)
alpha(4,i)= (C*E-B*F)/(A*E-B*D);	%Eq. (7.12b)
%acceleration of point A...aka point that connects links 2 and 3
AA=a*alpha(2,1)*(-sin(theta(2,i))+j*cos(theta(2,i)))...
	-a*omega(2,i)^2*(cos(theta(2,i))+j*sin(theta(2,i)));		%Eq. (7.13a)
%acceleration of point B...aka point that connects links 3 and 4...aka
%coupler point
AB=c*alpha(4,i)*(-sin(theta(4,i))+j*cos(theta(4,i)))...
	-c*omega(4,i)^2*(cos(theta(4,i))+j*sin(theta(4,i)));		%Eq. (7.13c)
%A_A+A_BA-A_B = 0
ABA = AB-AA;
% ----------------------------------------------------------------------------
i=i+1;		%Increses the counter to continue the while loop.
%Checks if the linkage completes a cycle. 
% ----------------------------------------------------------------------------
if (i>n+1)				%The linkage completes a cycle. 
     wflag=0;			%Signals the flag to end the while loop.
else
% ----------------------------------------------------------------------------
% Calculate theta_2, omega_2, and alpha_2 for the next position. 
% ----------------------------------------------------------------------------
omega(2,i)=omega(2,i-1); 												%Constant omega
theta(2,i)=theta(2,i-1)+delta_theta_2;								%Increases theta_2
%Alternative codes
%theta(2,i)=theta(2,1)+omega(2,1)*(delta_time*i)				%Constant omega	
%delta_theta_2=omega(2,i-1)*delta_time								%Constant omega
% Note that delta_time is not used in this program. However, delta_time is useful
%	if the results need to be compared with those from Working Model files.
% ----------------------------------------------------------------------------
%alpha(2,i)=alpha(2,i-1); 												%Constant alpha
%omega(2,i)=omega(2,i-1)+alpha(2,i-1); 							%Constant alpha
%theta(2,i)=theta(2,1)+.5*alpha(2,1)*(delta_time*i)^2			%Constant alpha
%To toggle between these two options, constant velocity or constant acceleration, 
% comment out the constant omega motion and un-comment constant angular alpha motion 
% ----------------------------------------------------------------------------
%Checks if the linkage is beyond its limits with the new theta_2.
% ----------------------------------------------------------------------------
D=cos(theta(2,i))-K1+K4*cos(theta(2,i))+K5;
E=-2*sin(theta(2,i));
F=K1+(K4-1)*cos(theta(2,i))+K5;
% if theta_2_input == theta_toggle
% 
if theta(2,1) == theta_toggle_d;
    wflag = 0;
end
if(E^2-4.*D*F<0)
    wflag = 0;
end
    
% 	if  (E^2-4.*D*F<0)				%The linkage reaches its moton limits.
%   		 wflag=0;						%Signal the flag to end the while loop.
% 	end 									%End if
% ----------------------------------------------------------------------------
end 		%End if
end		%End while 

% ----------------------------------------------------------------------------
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%The following section is for animation.
% ----------------------------------------------------------------------------
XC=[];YC=[];
I13X=[];I13Y=[];I24X=[];I24Y=[];
%output=menu('Simulation Options','4bar simulation','4bar simulation with coupler curve','4bar simulation with instant centers','4bar simulation with fixed centrode');
for k=1:i-1						%"for loop" to animate linkage motion. 
	%Animation stops before the linkage reaches its moton limits or a complete cycle.   
% ----------------------------------------------------------------------------
%Calculates the coordinates of joints and the coupler point
%	based on vector equation with nomenclature in Figure 4-7.
% ----------------------------------------------------------------------------
xo2=0;									%Coordinates of point O2
yo2=0;
xo4=xo2+d*cos(theta(1,1));			%Coordinates of point O4
yo4=yo2+d*sin(theta(1,1));	
	%theta(1,k) is 0 in all HW problems.
xa=xo2+a*cos(theta(2,k));			%Coordinates of point A
ya=yo2+a*sin(theta(2,k));
xb=xa+b*cos(theta(3,k));			%Coordinates of point B
yb=ya+b*sin(theta(3,k));	
	%theta(3,k) is calculated in the position analysis section.
%Alternative codes
%xb=xo4+c*cos(theta(4,k));			%Coordinates of point B
%yb=yo4+c*sin(theta(4,k));	
xc=xa+p*cos(theta(3,k)+delta_3);	%Coordinates of coupler point 
yc=ya+p*sin(theta(3,k)+delta_3);	   
% ----------------------------------------------------------------------------
%Creates the lines representing binary links and sides of the coupler link 
% connected to the coupler point.
% ----------------------------------------------------------------------------
x_link2=[xo2 xa];			%Line of link 2
y_link2=[yo2 ya];
x_link3=[xa xb];			%Line of link 3
y_link3=[ya yb];
x_link4=[xb xo4];			%Line of link 4
y_link4=[yb yo4];
x_link1=[xo4 xo2];		%Line of link 1
y_link1=[yo4 yo2];
x_c1=[xa xc];				%Line of coupler_link side 1
y_c1=[ya yc];
x_c2=[xb xc];				%Line of coupler_link side 2
y_c2=[yb yc];
XC=[XC,xc];
YC=[YC,yc];
% ----------------------------------------------------------------------------
% Calculate the coordinate of instant centers
% ----------------------------------------------------------------------------
I13x=[xb*(yb-yo4)*(xa-xo2)-xa*(ya-yo2)*(xb-xo4)+(ya-yb)*(xb-xo4)*(xa-xo2)]/[(yb-yo4)*(xa-xo2)-(ya-yo2)*(xb-xo4)];
I13y=[yb*(xb-xo4)*(ya-yo2)-ya*(xa-xo2)*(yb-yo4)+(xa-xb)*(yb-yo4)*(ya-yo2)]/[(xb-xo4)*(ya-yo2)-(xa-xo2)*(yb-yo4)];
I24x=[xb*(yb-ya)*(xo4-xo2)-xo4*(yo4-yo2)*(xb-xa)+(yo4-yb)*(xb-xa)*(xo4-xo2)]/[(yb-ya)*(xo4-xo2)-(yo4-yo2)*(xb-xa)];
I24y=[yb*(xb-xa)*(yo4-yo2)-yo4*(xo4-xo2)*(yb-ya)+(xo4-xb)*(yb-ya)*(yo4-yo2)]/[(xb-xa)*(yo4-yo2)-(xo4-xo2)*(yb-ya)];
I13X=[I13X,I13x];
I13Y=[I13Y,I13y];
I24X=[I24X,I24x];
I24Y=[I24Y,I24y];

% val = get(h,'Value');
% switch val
% case 2 
    % --------------------------------------------------------------------
% a1= str2double(get(handles.edit_a, 'String'));
% b1= str2double(get(handles.edit_b, 'String'));
% c1= str2double(get(handles.edit_c, 'String'));
% d1= str2double(get(handles.edit_d, 'String'));
% %theta1= str2double(get(handles.theta_1_button, 'String'));
% theta2= str2double(get(handles.edit_theta2, 'String'));
% delta= str2double(get(handles.edit_delta3, 'String'));
% p1= str2double(get(handles.edit_p, 'String'));

plot(xo2,yo2,'-.ko',x_link2,y_link2,xa,ya,'-.ko',x_link3,y_link3,xb,yb,'-.ko',x_link4,y_link4,xo4,yo4,'-.ko',x_link1,y_link1,x_c1,y_c1,x_c2,y_c2)
hold on;
xmin=(xo2-a);
xmax=(xo4+c);
ymin=max(yo2-a, yo4-c)-p1*cos(delta*pi/180);
ymax=max(yo2+a, yo4+c)+p1*cos(delta*pi/180);
% axis([xmin  xmax ymin ymax]*1.5)	
axis equal      
plot([xmin  xmax ]*1.4,[0 0], '-.k', [0 0],[ymin ymax]*1.2,'-.k')
hold off;
% case 3 % instant centers
%     
%  a1= str2double(get(handles.edit_a, 'String'));
% b1= str2double(get(handles.edit_b, 'String'));
% c1= str2double(get(handles.edit_c, 'String'));
% d1= str2double(get(handles.edit_d, 'String'));
% %theta1= str2double(get(handles.theta_1_button, 'String'));
% theta2= str2double(get(handles.edit_theta2, 'String'));
% delta= str2double(get(handles.edit_delta3, 'String'));
% p1= str2double(get(handles.edit_p, 'String'));
% 
% plot(xo2,yo2,'-.ko',x_link2,y_link2,xa,ya,'-.ko',x_link3,y_link3,xb,yb,'-.ko',x_link4,y_link4,xo4,yo4,'-.ko',x_link1,y_link1, I13x,I13y,'o', I24x,I24y,'o')
%       axis([xo2-a-2 xo4+c+2 yo2-a-2 yo2-a-2+(xo4+c-xo2+a+4)]*2)	
% 
% 
% 
% case 4 % coupler curve
% a1= str2double(get(handles.edit_a, 'String'));
% b1= str2double(get(handles.edit_b, 'String'));
% c1= str2double(get(handles.edit_c, 'String'));
% d1= str2double(get(handles.edit_d, 'String'));
% %theta1= str2double(get(handles.theta_1_button, 'String'));
% theta2= str2double(get(handles.edit_theta2, 'String'));
% delta= str2double(get(handles.edit_delta3, 'String'));
% p1= str2double(get(handles.edit_p, 'String'));
% 
% plot(xo2,yo2,'-.ko',x_link2,y_link2,xa,ya,'-.ko',x_link3,y_link3,xb,yb,'-.ko',x_link4,y_link4,xo4,yo4,'-.ko',x_link1,y_link1,x_c1,y_c1,x_c2,y_c2, xc,yc,'o', XC,YC)
%       axis([xo2-a-2 xo4+c+2 yo2-a-2 yo2-a-2+(xo4+c-xo2+a+4)]*2)	
% 
% case 5 % fixed centrode
% a1= str2double(get(handles.edit_a, 'String'));
% b1= str2double(get(handles.edit_b, 'String'));
% c1= str2double(get(handles.edit_c, 'String'));
% d1= str2double(get(handles.edit_d, 'String'));
% %theta1= str2double(get(handles.theta_1_button, 'String'));
% theta2= str2double(get(handles.edit_theta2, 'String'));
% delta= str2double(get(handles.edit_delta3, 'String'));
% p1= str2double(get(handles.edit_p, 'String'));
% 
% plot(xo2,yo2,'-.ko',x_link2,y_link2,xa,ya,'-.ko',x_link3,y_link3,xb,yb,'-.ko',x_link4,y_link4,xo4,yo4,'-.ko',x_link1,y_link1, I13x,I13y,'o', I13X,I13Y)
%       axis([xo2-a-2 xo4+c+2 yo2-a-2 yo2-a-2+(xo4+c-xo2+a+4)]*2)	
% 
% 

% end

% guidata(h,handles)

% axis([xo2-a-2 xo4+c+2 yo2-a-2 yo2-a-2+(xo4+c-xo2+a+4)]*2)	
   %axis([xmin xmax ymin ymax])
	% Automatically adjusts the range of axes to show the whole linkage in the screen 
	%   when a new linkage is created.
   % This is accomplised by making the width (xmax-xmin) equal to the height (ymax-ymin)  
   %  for proper aspect ratio.
% axis off 					%Turns off all axis lines, tick marks, and labels.
drawnow; 					%Completes pending drawing events and updates the Figure window
end 							%end "for loop"
% ----------------------------------------------------------------------------
theta=theta.*180/pi;			%Converts radians to degrees for output.




% % --- Executes on button press in pushbutton_instant_centers.
% function pushbutton_instant_centers_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton_instant_centers (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% a1= str2double(get(handles.edit_a, 'String'));
% b1= str2double(get(handles.edit_b, 'String'));
% c1= str2double(get(handles.edit_c, 'String'));
% d1= str2double(get(handles.edit_d, 'String'));
% % theta1= str2double(get(handles.theta_1_button, 'String'));
% theta1= 0;
% theta2= str2double(get(handles.edit_theta2, 'String'));
% delta= str2double(get(handles.edit_delta3, 'String'));
% p1= str2double(get(handles.edit_p, 'String'));
