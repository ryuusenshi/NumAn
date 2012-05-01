clc;
%%
%INITIALIZATION
%   seting up the grid dimensions, voltages of the container and the key
c_voltage = 0;
key_voltage = 5;
hor_dim = 51;
ver_dim = 51;

%%
%DEFININING THE CONTAINER
%   assumed rectangular
xv=[1 10 10 1 1];
yv=[1 1 6 6 1];
hold on;
plot(xv,yv);
axis equal;

%%
%DEFINING THE KEY
%   assumed a polygon
xv2=[5 7 7 5 5 4 3.75 4 5 5];
yv2=[4 4 5 5 4.75 4.75 4.5 4.25 4.25 4];
plot(xv,yv);
plot(xv2,yv2);


%should correct this to work even if the container is not rectangulat
%   define the grid and then use only grid points that are inside the
%   container!

%%
%DEFINING GRID POINTS AND PLOTTING THEM
% (blue - outside the key)
% (red - inside the key)

x=linspace(min(xv),max(xv),hor_dim);
dx=(max(xv)-min(xv))/(hor_dim-1);
xh=ones(hor_dim^2,1);
for i=0:hor_dim-1
    xh(i*hor_dim+1:(i+1)*hor_dim) = x(i+1)*ones(hor_dim,1);
end

y=linspace(min(yv),max(yv),ver_dim);
dy=(max(yv)-min(yv))/(ver_dim-1);
yh=ones(ver_dim^2,1);
for i=0:ver_dim-1
    yh(i*ver_dim+1:(i+1)*ver_dim) = y(1:ver_dim);
end
in = inpolygon(xh,yh,xv2,yv2);
plot(xh(in),yh(in),'.r',xh(~in),yh(~in),'.b');



%%
%SETING UP THE POTENTIAL MATRIX
% v=c_voltage - outside the key
% v=key_voltage - inside or on the key

v= c_voltage * ones(hor_dim,ver_dim);
indices = find(in==1);
indx=floor(indices/hor_dim)+1;
indy=mod(indices,ver_dim);
for i=1:length(indx)
    v(indy(i),indx(i))=key_voltage;
end


fixed = zeros(hor_dim,ver_dim);
fixed(1,:) = 1;
fixed(hor_dim,:)=1;
fixed(:,1) = 1;
fixed(:,ver_dim)=1;

for i=1:length(indx)
    v(indy(i),indx(i))=key_voltage;
    fixed(indy(i),indx(i))=1;
end


%%
% ITERATIVE CALCULATION OF POTENTIAL
%
figure;
surf(x,y,v);
pause;
%while(1)
for loopn = 1:50
for i=1:ver_dim
    mask = ones(3,3);
    mask(2,2) = 0;
    for j=1:hor_dim
        if (fixed(i,j)==0)
            if (i+1>ver_dim)
                mask(3,:) = 0;
            end
            if (i-1<1)
                mask(1,:) = 0;
            end
            if (j+1>hor_dim)
                mask(:,3) = 0;
            end
            if (j-1<1)
                mask(:,1) = 0;
            end

            [xh,yh]=find(mask==1);

            l=[];
            for k=1:length(xh)
                l=[l v(xh(k)+i-2,yh(k)+j-2)];
            end
            v(i,j)=mean(l);
        end
    end
end
%surf(x,y,v);
%pause;
end

surf(x,y,v);

%%
%CALCULATION OF ELECTRIC FIELD (E=-grad(V))
[fx,t]=gradient(v/dx);
[t,fy]=gradient(v/dy);
Ex=-fx
Ey=-fy


