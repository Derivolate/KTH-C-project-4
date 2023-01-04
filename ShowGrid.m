close all;
x_data = fopen("outfileX.bin");
y_data = fopen("outfileY.bin");
u_data = fopen("u_out.bin");
udx_data = fopen("dx_out.bin");
udy_data = fopen("dy_out.bin");
ddxy_data = fopen("ddxy_out.bin");
x = fread(x_data, 'double');
y = fread(y_data, 'double');
u = fread(u_data, 'double');
udx = fread(udx_data, 'double');
udy = fread(udy_data, 'double');
ddxy = fread(ddxy_data, 'double');
discrete_functions = {u, udx, udy, ddxy};
coloraxis = {[-Inf,Inf],[-Inf,Inf],[0.75,1.25],[-Inf,Inf]};
titles= ["U Function", "du/dx", "du/dy", "Laplacian"];
ax = [-10 5 0 3];
m = 20;
n = 50;
sz = 20;
border3 = [x(1:m),y(1:m)];
border1 = [x(n*m-m+1:n*m),y(n*m-m+1:n*m)];
border0 = [x(1:m:n*m-m+1), y(1:m:n*m-m+1)];
border2 = [x(m:m:n*m), y(m:m:n*m)];

a_u = sin((x./10).^2).*cos(x./10)+y;
a_udx = 1/50.*(x.*cos(x./10).*cos(x.^2./100) - 5.*sin(x./10).*sin(x.^2./100));
a_udy = ones(size(x));
a_laplace = (-10.*x.*sin(x./10).*cos(x.^2./100)-cos(x./10).*((x.^2.+25).*sin(x.^2./100)-50.*cos(x.^2./100)))/2500;
% a_u = x.*x;
% a_udx = 2*x;
% a_udy = zeros(size(x));
% a_laplace = ones(size(x))*2;
analytical_functions = {a_u,a_udx,a_udy,a_laplace};


% ---------------------------- Plot 1 --------------------------------- 
figure
for k = 1:4 
    subplot(2,2,k);
    hold on
    title(titles(k));
    for i = m+1:m:n*m-m+1
        vline = [x(i:i+m-1),y(i:i+m-1)];
        plot(vline(:,1),vline(:,2), 'k-');
    end
    for j = 2:m
        hline = [x(j:m:j+(n-1)*m), y(j:m:j+(n-1)*m)];
        plot(hline(:,1),hline(:,2), 'k-');
    end
    plot(border0(:,1),border0(:,2), 'm-');
    plot(border1(:,1),border1(:,2), 'b-');
    plot(border2(:,1),border2(:,2), 'g-');
    plot(border3(:,1),border3(:,2), 'r-');
    scatter(x,y,sz,discrete_functions{k}, 'filled');
    axis(ax);
    caxis(coloraxis{k});
    colorbar
end

% ---------------------------- Plot ERROR --------------------------------- 


figure
for k = 1:4 
    subplot(2,2,k);
    f = discrete_functions{k}-analytical_functions{k};
    hold on
    title(titles(k));
    for i = m+1:m:n*m-m+1
        vline = [x(i:i+m-1),y(i:i+m-1)];
        plot(vline(:,1),vline(:,2), 'k-');
    end
    for j = 2:m
        hline = [x(j:m:j+(n-1)*m), y(j:m:j+(n-1)*m)];
        plot(hline(:,1),hline(:,2), 'k-');
    end
    plot(border0(:,1),border0(:,2), 'm-');
    plot(border1(:,1),border1(:,2), 'b-');
    plot(border2(:,1),border2(:,2), 'g-');
    plot(border3(:,1),border3(:,2), 'r-');
    scatter(x,y,sz,f, 'filled');
    axis(ax);
    colorbar
end



figure
subplot(2,2,1);
scatter3(x,y,u,'filled');


subplot(2,2,2);
scatter3(x,y,udx,'filled');

subplot(2,2,3);
scatter3(x,y,udy,'filled');

subplot(2,2,4);
scatter3(x,y,ddxy,'filled');





