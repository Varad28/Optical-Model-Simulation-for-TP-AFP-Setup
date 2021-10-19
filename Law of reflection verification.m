%Varad Kulkarni
%13.12.2020
%Verifying Vector Treatment of Law of Reflection of Light


clc
clear all
close all

%Calculating reflected ray vector 

theta=input('Enter Angle of Incidence with Surface: ');
fp=input('Enter Co-ordinates of focus point:(in []): ');

x=(fp(2)\1)*(fp(2)*fp(1)-(fp(2)^2/tand(theta)));
%fprintf('\nx is: %d; ',x);
i=[x-fp(1),-fp(2)];
n=[0,1];

r=i-2*dot(i,n)*n; %Formula to be verified  
g=sprintf('%d ', r);
fprintf('\nThe reflected vector r is: [ %s]\n',g);

a=abs(dot(i,n)); 
b=abs(dot(r,n)); %Calculating absolute values of dot products; 'absolute' to always measure accute angles

mod_i=norm(i); 
mod_n=norm(n); 
mod_r=norm(r); %Calculating magnitudes of vectors

theta_i=acosd((mod_i*mod_n)\a); 
theta_r=acosd((mod_r*mod_n)\b); %Calculating angle b/w vectors using dot product equation: a.b=|a||b|cos(theta)
fprintf('Theta i is: %d',theta_i);
fprintf('\nTheta r is: %d\n',theta_r);

%Ray Diagram

ui=mod_i.\i; ur=mod_r.\r; %Calculating unit vectors

figure('name','Ray Diagram');

v_i=quiver(x-ui(1),-ui(2),ui(1),ui(2),0,'r'); %q=quiver(X,Y,U,V,S,C) - X,Y-coordinates of tail, U,V-components, S-Scaling factor, C- colour
hold on
set (v_i, "maxheadsize", 0.033);

v_r=quiver(x,0,ur(1),ur(2),0,'g');
set (v_r, "maxheadsize", 0.033);

v_n=quiver(x,0,n(1),n(2),0.4,'--k'); % -- is line style
set (v_n, "maxheadsize", 0.033); 

legend('i','r','n','Fontsize',24,'location','northeastoutside')
t=title('Ray Diagram of Input and Calculated  Vectors','color','b'); set(t,"fontsize",24); set(t,'fontname','Times');

%text(0.06,0.5,['\angle(i,n)= ',num2str(theta_i),'\deg'],'Fontsize', 24) %num2str - to print whichever value of theta is calculated
%text(-0.3,0.5,['\angle(r,n)= ',num2str(theta_r),'\deg'],'Fontsize', 24)
%grid on
hold off
