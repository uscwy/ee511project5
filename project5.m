%EE511 Project 5
%Author: Yong Wang <yongw@usc.edu>
x=linspace(0,1);
fx = 0.6.*(1-x).^7./beta(1,8) + 0.4.*x.^8./beta(9,1);
figure;
plot(x,fx);

k=2000;
x=zeros(k,1);
v=.2;
x(1)=rand; %the intial point
for t=1:k-1
	%generate from normal(xt,1)
    xp = x(t) + v*randn;
    if xp < 0 || (xp > 1)
        fp = 0;
    else
        fp = 0.6.*(1-xp).^7./beta(1,8) + 0.4.*xp.^8./beta(9,1);
    end
    f = 0.6.*(1-x(t)).^7./beta(1,8) + 0.4.*x(t).^8./beta(9,1);
    a = fp*normpdf(x(t),xp)/f*normpdf(xp,x(t));
    if a > 1
        a = 1;
    end
    if a >= rand
        %accept
        x(t+1) = xp;
    else
        %reject
        x(t+1) = x(t);
    end
end
figure;
y=linspace(0,k,k);
plot(x,y);
xlabel('Sample Value');
ylabel('Iterations');
title(sprintf('Sample Path with var=%f',v));
figure;
histogram(x);
xlabel('X');
ylabel('Frequency');
title('Histogram of samples');
%% Schwefel Function
x=linspace(-500,500);
y=linspace(-500,500);
[X,Y] = meshgrid(x,y);
Z = 418.9829*2 - X.*sin(sqrt(abs(X))) - Y.*sin(sqrt(abs(Y)));
contour(x,y,Z);
xlabel('x');
ylabel('y');
title('Sample Path');

iter = 1000;
minima=zeros(100,1);
t0=8000;
for j=1:100 
    x=zeros(iter+1,2);
    %initial point
    x(1,:) = -500+1000*rand(1,2);
    for i=1:iter 
        %cooling schedule
        t = t0*0.99^(i-1);  %exponential
        t = t0/i; %polynomial
        t = t0/(1+80*log(i)); %logarithmic
        v = 400;
        xp = x(i,:);
        if(mod(i,2) == 0) 
            %generate new x
            xp(1) = xp(2)+v*randn;
        else
            %generate new y
            xp(2) = xp(1)+v*randn;
        end
        if(abs(xp(1))>500)
            xp(1)=sign(xp(1))*500;
        end
        if(abs(xp(2))>500)
            xp(2)=sign(xp(2))*500;
        end
        z=418.9829*2 - x(i,1)*sin(sqrt(abs(x(i,1)))) - x(i,2)*sin(sqrt(abs(x(i,2))));
        zp=418.9829*2 - xp(1)*sin(sqrt(abs(xp(1)))) - xp(2)*sin(sqrt(abs(xp(2))));

        a = exp(-(zp-z)/t);
        if  zp < z ||  rand < a
            x(i+1,:)=xp;
            m=zp;
        else
            x(i+1,:)=x(i);
            m=z;
        end
    end
    minima(j) = m;
end
hold on;
plot(x(:,1),x(:,2),'r');
plot(x(iter+1,1),x(iter+1,2),'bx');
figure;
histogram(minima);
xlabel('Minima');
ylabel('Frequency');
title(sprintf('Minimas with iteration=%d, t0=%d exponential cooling',iter,t0));
