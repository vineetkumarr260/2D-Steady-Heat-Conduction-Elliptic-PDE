%% Initialization %%

clear all;
close all;
clc;

material = ('Silver');
density = 10490;  % kg/m^3
conductivity = 429;   % W/(mK)
heat_capacity = 235;  % J/(kgK)

% material = ('Copper');
% density = 8960;  % kg/m^3
% conductivity = 401;   % W/(mK)
% heat_capacity = 384.4;  % J/(kgK)

% material = ('Aluminium');
% density = 2700;  %kg/m^3
% conductivity = 237;   %W/(mK)
% heat_capacity = 904;  %J/(kgK)

% material = ('Gold');
% density = 19300;  kg/m^3
% conductivity = 318;   W/(mK)
% heat_capacity = 129.1;  J/(kgK)

alpha = conductivity/(density*heat_capacity);

L_x = 1;
L_y = 1;
N_x = 20;
N_y = 20;
t_end = 1000;

dx = L_x/N_x;
dy = L_y/N_y;
dt_ideal = 1/(2*alpha*(1/dx^2+1/dy^2));
fprintf('Critical time-step size : %.2f\n\n',dt_ideal);

dt = 3.58;  % CHANGE TIME-STEP SIZE HERE %

% Checking stability condition
if(dt>1/(2*alpha*(1/dx^2+1/dy^2)))
    txt=sprintf('Stability condition is not satisfied. Please chose a value for dt smaller than %.2f\n',dt_ideal);
    errordlg(txt,'Stability Error');
    return;
end

% Initial condition and boundary conditions
T_initial = 120;
T_top = 70;    % RIGHT
T_bottom = 150; % LEFT
T_left = 50;    % BOTTOM
T_right = 100;  % TOP

% Applying initial and boundary conditions to the grid
T = zeros(N_x+2,N_y+2,75000);
T(:,1,:) = T_bottom;
T(1,:,:) = T_left;
T(:,N_y+1,:) = T_top;
T(N_x+1,:,:) = T_right;
T(:,:,1)=T_initial;

% Arrays to store and update steady state solution at the end of each iteration
T1=zeros(N_x+2,N_y+2);  % For Jacobi Iterative Method
T2=zeros(N_x+2,N_y+2);  % For Jacobi Iterative Method
T3=zeros(N_x+2,N_y+2);  % For Gauss-Seidel Method
T4=zeros(N_x+2,N_y+2);  % For Gauss-Seidel Method
T5=zeros(N_x+2,N_y+2);  % For Successive Over Relaxation Method
T6=zeros(N_x+2,N_y+2);  % For Successive Over Relaxation Method

T1(:,1) = T_bottom;
T1(1,:) = T_left;
T1(:,N_y+1) = T_top;
T1(N_x+1,:) = T_right;

T2(:,1) = T_bottom;
T2(1,:) = T_left;
T2(:,N_y+1) = T_top;
T2(N_x+1,:) = T_right;

T3(:,1) = T_bottom;
T3(1,:) = T_left;
T3(:,N_y+1) = T_top;
T3(N_x+1,:) = T_right;

T4(:,1) = T_bottom;
T4(1,:) = T_left;
T4(:,N_y+1) = T_top;
T4(N_x+1,:) = T_right;

T5(:,1) = T_bottom;
T5(1,:) = T_left;
T5(:,N_y+1) = T_top;
T5(N_x+1,:) = T_right;

T6(:,1) = T_bottom;
T6(1,:) = T_left;
T6(:,N_y+1) = T_top;
T6(N_x+1,:) = T_right;

tolerance_sst = 0.005;  % tolerance value for steady state solution
tolerance = 0.5;    % tolerance value for finite difference methods i.e Euler and Range-Kutta Method

%% Steady Solution %%

% Jacobi Iterative Method
JIM = 1;    % Iteration counter
error_max = 100;    % Initial error value (over-predicted)
ejim(1) = 0;    % Array to store error values of successive iterations
while (error_max >= tolerance_sst)
    for i = 2:N_x
        for j = 2:N_y
            T2(i,j) = (T1(i+1,j)+T1(i-1,j)+T1(i,j-1)+T1(i,j+1))/4;
        end
    end
    JIM = JIM+1;
    error_max = max(max(T2-T1));
    ejim(JIM) = error_max;
    T1 = T2;
end

% Gauss-Siedel Method
GSM = 1;
error_max = 100;
egsm(1) = 0;
while(error_max >= tolerance)
    for i = 2:N_x
        for j = 2:N_y
            T4(i,j) = (T4(i+1,j)+T4(i-1,j)+T4(i,j-1)+T4(i,j+1))/4;
        end
    end
    GSM = GSM+1;
    error_max = max(max(T4-T3));
    egsm(GSM) = error_max;
    T3 = T4;
end

% Successive Over Relaxation Method
SOR = 1;
lambda = 2-2*pi/N_x;    % Optimum relaxation factor for a square grid
error_max = 100;
esor(1) = 0;
while(error_max >= tolerance)
    for i = 2:N_x
        for j = 2:N_y
            T6(i,j) = lambda*0.25*(T6(i+1,j)+T6(i-1,j)+T6(i,j-1)+T6(i,j+1))+(1-lambda)*T6(i,j);
        end
    end
    SOR = SOR+1;
    error_max = max(max(T6-T5));
    esor(SOR) = error_max;
    T5 = T6;
end

%-------------- Steady State Solution Convergence Plots ------------------%
figure('Name','Steady State Solution Convergence Plots','NumberTitle','off','units','normalized','outerposition',[0 0 1 1])

subplot(2,2,1)
plot(2:JIM,ejim(2:JIM),'r');
title('Jacobi Iterative Method');
xlabel('Number of Iterations \rightarrow');
ylabel('Absolute Error \rightarrow');

subplot(2,2,2)
plot(2:GSM,egsm(2:GSM),'g');
title('Gauss-Seidel Method');
xlabel('Number of Iterations \rightarrow');
ylabel('Absolute Error \rightarrow');

subplot(2,2,3)
plot(2:SOR,esor(2:SOR),'b');
title('Successive Over Relaxation Method');
xlabel('Number of Iterations \rightarrow');
ylabel('Absolute Error \rightarrow');

subplot(2,2,4)
plot(2:JIM,ejim(2:JIM),'r');
hold on;
plot(2:GSM,egsm(2:GSM),'g');
plot(2:SOR,esor(2:SOR),'b');
legend('Jacobi Iterative','Gauss-Seidel','SOR');
hold off;

%% Finite-Difference Methods %%

dlgTitle    = 'Select Method';
dlgQuestion = 'Choose the numerical method utilizing finite-difference scheme to solve the associated PDE?';
choice = questdlg(dlgQuestion,dlgTitle,'Euler','2nd Order Runge-Kutta','Euler');

message = msgbox('Your computer is now iterating to solve the PDE, Please wait... '); % Busy message

% Euler Method
k = 1;
if strcmp(choice, 'Euler')
    error_max(k) = 100;
    error_min(k) = 100;
    while(error_max(k) >= tolerance || error_min(k) >= tolerance)
        for i = 2:N_x
            for j = 2:N_y
                T(i,j,k+1) = T(i,j,k)+dt*alpha*(((T(i-1,j,k)-2*T(i,j,k)+T(i+1,j,k))/dx^2)+((T(i,j-1,k)-2*T(i,j,k)+T(i,j+1,k))/dy^2));
            end
        end
        k=k+1;
        error_max(k) = abs(max(max(T(:,:,k)-T1)));
        error_min(k) = abs(min(min(T(:,:,k)-T1)));
        if round(error_min(k),5)==round(error_min(k-1),5) && error_min(k)~=0
            errordlg('The solution is not converging. Please choose a larger tolerance value.','Tolerence Error');
            return;
        end
        if round(error_max(k),5)==round(error_max(k-1),5) && error_max(k)~=0
            errordlg('The solution is not converging. Please choose a larger tolerance value.','Tolerence Error');
            return;
        end
    end
else
     % Range-Kutta Method    
     error_max(k) = 100;
     error_min(k) = 100; 
     while error_max(k) >= tolerance || error_min(k) >= tolerance
         for i = 2:N_x
             for j = 2:N_y
                 k1 = alpha*(((T(i-1,j,k)-2*T(i,j,k)+T(i+1,j,k))/dx^2)+((T(i,j-1,k)-2*T(i,j,k)+T(i,j+1,k))/dy^2));
                 Tk = T(:,:,k)+k1*dt;
                 k2 = alpha*(((Tk(i-1,j)-2*Tk(i,j)+Tk(i+1,j))/dx^2)+((Tk(i,j-1)-2*Tk(i,j)+Tk(i,j+1))/dy^2));
                 T(i,j,k+1) = T(i,j,k)+(dt/2)*(k1+k2);
            end
         end
         k = k+1;
         error_max(k) = abs(max(max(T(:,:,k)-T1))); 
         error_min(k) = abs(min(min(T(:,:,k)-T1))); 
         if round(error_max(k),5)==round(error_max(k-1),5) && error_max(k)~= 0
             errordlg('The solution is not converging. Please choose a larger tolerance value.','Tolerence Error');
             close(message)
             return
         end
         if round(error_min(k),5)==round(error_min(k-1),5) && error_min(k)~= 0
             errordlg('The solution is not converging, Please choose a larger tolerance value.','Tolerence Error');
             close(message)
             return
         end
     end
end

close(message);

% Free the unallocated array space %
T=T(:,:,1:k);

sol_time = k*dt;

T_max = max([T_top,T_bottom,T_right,T_left,T_initial]);
T_min = min([T_top,T_bottom,T_right,T_left,T_initial]);

%% ------------------------------ SOLUTION ----------------------------- %%
fprintf('Solution for plate of dimensions %im x %im made of material %s \n',L_x,L_y,material);
fprintf('Thermal conductivity of %s : %.1f W/(mK) \n',material,conductivity);
fprintf('Plate was subjected to Dirichlet boundary conditions : \n');
fprintf('T(x,y,0) = %i, T(x,0,t) = %i, T(x,L,t) = %i, T(0,y,t) = %i, T(W,y,t) = %i \n',T_initial,T_left,T_right,T_bottom,T_top);
fprintf('Time taken to reach steady state temperature is %.2f seconds under tolerance of %.2f \n', sol_time,tolerance);
fprintf('Jacobi Iterative Method took %i iterations to converge.\n',JIM);
fprintf('Gauss-Seidel Method took %i iterations to converge.\n',GSM);
fprintf('Successive Over Relaxation Method took %i iterations to converge.\n\n',SOR);
%-------------------------------------------------------------------------%

%% Generating the plate %%
x=zeros(1,N_x+2);
y=zeros(1,N_y+2);
for i = 1:N_x+2
x(i) =(i-1)*dx;
end
for i = 1:N_y+2
y(i) =(i-1)*dy;
end

%% -------------------------- Plotting Section ------------------------- %%
figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,2,1)
surf(x,y,T1);
title('Temperature Distribution At Steady State');
xlabel('Length');
xlim([0 L_x+dx]);
ylabel('Width');
ylim([0 L_y+dy]);
zlabel('Temperature');
zlim([T_min T_max]);
colorbar;
caxis([T_min T_max]);
hold on;
plot3(L_x/4,L_y/4,T_max,'ko','markerfacecolor','r');
plot3(L_x/4,L_y/4,T_min,'ko','markerfacecolor','r');
plot3(L_x/2,L_y/2,T_max,'ko','markerfacecolor','g');
plot3(L_x/2,L_y/2,T_min,'ko','markerfacecolor','g');
plot3(3*L_x/4,3*L_y/4,T_max,'ko','markerfacecolor','b');
plot3(3*L_x/4,3*L_y/4,T_min,'ko','markerfacecolor','b');
view(0,90);
drawnow;
hold off;

subplot(2,2,2)
val1 = num2str(T1(round(N_x/4),round(N_y/4)),'%.2f');
txt1 = strcat(' T = ',val1);
val2 = num2str(T1(round(N_x/2),round(N_y/2)),'%.2f');
txt2 = strcat(' T = ',val2);
val3 = num2str(T1(round(3*N_x/4),round(3*N_y/4)),'%.2f');
txt3 = strcat(' T = ',val3);
hold on;
xlabel('Number of Iterations');
ylabel('Temperature');
title('Steady State Temperature At Selected Points');
grid on;
xlim([0 2*k]);
ylim([T_min T_max]);
scatter(k,T1(round(N_x/4),round(N_y/4)),'ko','markerfacecolor','r');
text(k,T1(round(N_x/4),round(N_y/4)),txt1,'HorizontalAlignment','left');
scatter(k,T1(round(N_x/2),round(N_y/2)),'ko','markerfacecolor','g');
text(k,T1(round(N_x/2),round(N_y/2)),txt2,'HorizontalAlignment','left');
scatter(k,T1(round(3*N_x/4),round(3*N_y/4)),'ko','markerfacecolor','b');
text(k,T1(round(3*N_x/4),round(3*N_y/4)),txt3,'HorizontalAlignment','left');
legend('Red Point','Green Point','Blue Point');
drawnow;
hold off;

for l=1:k

subplot(2,2,3)
surf(x,y,T(:,:,l));
xlabel('Length');
xlim([0 L_x+dx]);
ylabel('Width');
ylim([0 L_y+dy]);
zlabel('Temperature');
zlim([T_min T_max]);
colorbar;
caxis([T_min T_max]);
title(sprintf('Temperature Distribution At Iteration %i And Time = %.0f seconds',l,l*dt));
hold on;
plot3(L_x/4,L_y/4,T_max,'ko','markerfacecolor','r');
plot3(L_x/4,L_y/4,T_min,'ko','markerfacecolor','r');
plot3(L_x/2,L_y/2,T_max,'ko','markerfacecolor','g');
plot3(L_x/2,L_y/2,T_min,'ko','markerfacecolor','g');
plot3(3*L_x/4,3*L_y/4,T_max,'ko','markerfacecolor','b');
plot3(3*L_x/4,3*L_y/4,T_min,'ko','markerfacecolor','b');
view(0,90);
drawnow;
hold off;

subplot(2,2,4)
hold on;
xlabel('Number of Iterations');
ylabel('Temperature');
title('Temperature Convergence At Selected Points');
grid on;
xlim([0 l+1]);
scatter(l,T(round(N_x/4),round(N_y/4),l),'r.');
scatter(l,T(round(N_x/2),round(N_y/2),l),'g.');
scatter(l,T(round(3*N_x/4),round(3*N_y/4),l),'b.');
legend('Red Point','Green Point','Blue Point');
drawnow;
hold off;

end
%-----------------------------------------------------------------------%