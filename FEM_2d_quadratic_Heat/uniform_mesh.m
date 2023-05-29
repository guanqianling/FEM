function [p1, s1, t1]=uniform_mesh(n,left,right,draw)

if nargin<4
    draw=0;
end
if nargin<2
    left=0;
    right=1;
end
if nargin<1
    n=4;
end

h=(right-left)/n;
p=zeros(3,(n+1)^2);
t=zeros(3,2*n^2);
for i=1:n+1
    for j=1:n+1
        m=tr(i,j,n);
        p(1,m)=right-(j-1)*h;
        p(2,m)=left+(i-1)*h;
        p(3,m)=0;
        if (i-1)*(j-1)*(i-n-1)*(j-n-1)==0
            p(3,m)=1;
        end
    end
end
for i=1:n
    for j=1:n
        m=trs(i,j,n);
        k=2*n^2-trs1(n+1-i,n+1-j,n)+1;
        t(:,m)=[tr(i,j,n);tr(i+1,j,n);tr(i,j+1,n)];
        t(:,k)=[tr(i,j+1,n);tr(i+1,j,n);tr(i+1,j+1,n)];
    end
end
e=zeros(2,4*n);
for i=1:n
    e(:,i)=[tr(i,1,n);tr(i+1,1,n)];
    e(:,n+i)=[tr(n+1,i,n);tr(n+1,i+1,n)];
    e(:,2*n+i)=[tr(n+2-i,n+1,n);tr(n+1-i,n+1,n)];
    e(:,3*n+i)=[tr(1,n+2-i,n);tr(1,n+1-i,n)];
end
[s, t]=sort_edge(e,t);
for i=1:size(t,2)
    [~,order]=sort(-t(1:3,i));
    order_inv(order)=1:3;
    t(4:6,i)=t(3+order_inv,i);
end

for i=1:size(s,2)
    if s(4,i)==0
        if p(1,s(1,i))+p(1,s(2,i))==0
            s(4,i)=-2;
        elseif p(1,s(1,i))+p(1,s(2,i))==2
            s(4,i)=-3;
        else
            s(4,i)=-1;
        end
    end
end

if draw==1
    hold on
    for i=1:size(p,2)
        temp=num2str(i);
        text(p(1,i)+0.1*h,p(2,i)-0.1*h,temp,'Color','blue');
        plot(p(1,i),p(2,i),'r*');
    end
    for i=1:size(s,2)
        temp=num2str(i);
        text(sum(p(1,s(1:2,i)))/2,sum(p(2,s(1:2,i)))/2,temp,'Color','green');
    end
    for i=1:size(t,2)
        temp=num2str(i);
        text(sum(p(1,t(1:3,i)))/3,sum(p(2,t(1:3,i)))/3,temp,'Color','red');
        plot(p(1,t([1 2],i)),p(2,t([1 2],i)));
        plot(p(1,t([2 3],i)),p(2,t([2 3],i)));
        plot(p(1,t([1 3],i)),p(2,t([1 3],i)));
    end
    hold off
end

p1=[p(1:2,:) (p(1:2,s(1,:))+p(1:2,s(2,:)))/2];
s1=s;
t1=t+[zeros(3,k);(n+1)^2*ones(3,k)];

function m=tr(i,j,n)
if i+j<n+3
    m=(i+j-2)*(i+j-1)/2+j;
    if mod(i+j,2)==1
        m=(2*i+2*j-2)*(i+j-1)/2+1-m;
    end
else
    m=(n+1)^2-tr(n+2-i,n+2-j,n)+1;
end

function m=trs(i,j,n)
if i+j<n+2
    m=(i+j-2)^2+2*j-1;
    if mod(i+j,2)==1
        m=(i+j-1)^2+(i+j-2)^2+1-m;
    end
else
    m=2*n^2-((2*n+1-i-j)^2+2*(n+1-j))+1;
    if mod(i+j,2)==1
        m=4*n^2-((2*n+2-i-j)^2+(2*n+1-i-j)^2-1)-m;
    end
end

function m=trs1(i,j,n)
if i+j<n+2
    m=(i+j-2)^2+2*j-1;
    if mod(i+j,2)==0
        m=(i+j-1)^2+(i+j-2)^2+1-m;
    end
else
    m=2*n^2-((2*n+1-i-j)^2+2*(n+1-j))+1;
    if mod(i+j,2)==0
        m=4*n^2-((2*n+2-i-j)^2+(2*n+1-i-j)^2-1)-m;
    end
end