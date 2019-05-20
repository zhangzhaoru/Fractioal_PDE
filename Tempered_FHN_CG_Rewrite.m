%---------------���ڹ����ݶȷ���FHNģ�����޲�ָ�ʽ�����-----------------------------
function Tempered_FHN_CG_Rewrite()
%ʱ�������1�����Euler��ʽ
%�ռ������2��GL��ֵ�ƽ�
%���ù����ݶȷ���CG��
%ʹ��ADI��ʽ���㣬�����ٶ���Խ���

clc
clear
% clf
m1 = 256;     % x�������ʷֽڵ���
m2 = 256;     % y�������ʷֽڵ���
%��Ϊż������㣬����2�Զ��ҵ��������ڵ���е�

m1_half = m1/2;
m2_half = m2/2;
m1_quarter=m1/4;
m2_quarter=m2/4;
N = m1_half*m2_half;
M = (m1 - 1)*(m2 - 1); % ϵ������A��ά��
% r1=1.25;
% r2=0.625;


%��Բ�뾶
%��Բ�뾶


alpha1 = 1.7  ;   % x����������
alpha2 = 1.7  ;   % y����������
Dx = 2.5;       % x�����Ҷ˵�
Dy = 2.5;       % y�����϶˵�
Kx = 1e-4;    % x������ɢϵ��
Ky = 1e-4;    % y������ɢϵ��
hx = Dx/m1;   % x�������񲽳�
hy = Dy/m2;   % y�������񲽳�
tau = 0.1;   % ʱ�䲽��
t = 1000;     % ��ֹʱ��
layend = t / tau;  % ʱ���ѭ������

%-----FHNģ���в�����ѡȡ---------%
a=0.1;
epsilon = 0.01;
beta=0.5;
gamma = 1;
delta=0;
lamda1=0;
lamda2=lamda1;
% �ػ�����׿�������

Ca1 = 1.0/(2 * cos(pi * alpha1 / 2)); % x����Rieze�����׵��������ϵ��
Ca2 = 1.0/(2 * cos(pi * alpha2 / 2)); % y����Rieze�����׵��������ϵ��

r1 = (-1) * (tau* Kx * Ca1)/(hx^alpha1);
r2 = (-1) * (tau * Ky * Ca2)/(hy^alpha2);

%-----x��y�����ϵĽڵ�ֵ------%
xx1=linspace(0,Dx,m1-1);
%1*m1-1
yy1=linspace(0,Dy,m2-1);

xx=repmat(xx1,1,m1-1)';
%(m1-1)*(m1-1),1 cycle
yy=repmat(yy1,m2-1,1);
%(m2-1),(m2-1)  �и���
yy=reshape(yy,M,1);
%(m2-1)*(m2-1),1  Staggered cycle

%% initial values of u,v
u = zeros(m1-1,m2-1);
v = zeros(m1-1,m2-1);



for i = 1:m1-1
    for j = 1:m2-1
        if( ( i < m1_half ) && ( j <= m2_half ) ) % u�������ķ�֮һ����ֵΪ1
            u( i,j ) = 1.0;
        else
            u( i ,j  ) = 0.0;
        end
        if( ( i >= m1_half ) && ( i < m1 ) )    % v���ϰ�����Ϊ1
            v( i,j  ) = 0.1;
        else
            v( i,j) = 0.0;
        end
    end
end



%% initial Gruwald weights g(L)
g1 = zeros(1,m1+1);
g2 = zeros(1,m1+1);
g1(1) = 1.0;
g2(1) = 1.0;
for L = 1:m1
    g1(L+1)=-1*g1(L)*(alpha1-L+1)/L;
    g2(L+1)=-1*g2(L)*(alpha2-L+1)/L;
end


ge1 = zeros(1,m1+1);
ge2 = zeros(1,m1+1);
ge1(1) = 1.0 .* exp(-(1-1).*hx.*lamda1);
ge2(1) = 1.0 .* exp(-(1-1).*hy.*lamda2);
for L = 1:m1
    ge1(L+1)=-1*ge1(L)*(alpha1-L+1)/L* exp(-(L).*hx.*lamda1);
    ge2(L+1)=-1*ge2(L)*(alpha2-L+1)/L* exp(-(L).*hy.*lamda2);
end

%����m1+1��Ȩ��ϵ��
%
% % w1 = zeros(1,m1+1);
% % w2 = zeros(1,m1+1);
% % w1(1) = (alpha1/2)*g1(1);
% % w2(1) = (alpha2/2)*g1(1);
% % for L = 2:m1+1
% %     w1(L) = (alpha1/2)*g1(L)+(1-alpha1/2)*g1(L-1);
% %     w2(L) = (alpha1/2)*g2(L)+(1-alpha2/2)*g2(L-1);
% % end
% %
% % %����ʱ��2��ϵ��
%���Ĳ�����ӵ�Ȩ��ϵ��
%
% g1 = zeros(1,m1+1);
% g2 = zeros(1,m1+1);
% % g1(1) = 1.0;
% % g2(1) = 1.0;
%
% a = gamma(alpha1+1)/gamma(alpha1/2+1)^2;
% b = gamma(alpha2+1)/gamma(alpha2/2+1)^2;
% for L = 1:m1
%     g1(L+1)=(1-(alpha1+1)/(alpha1/2.0+L))*g1(L);
%     g2(L+1)=(1-(alpha2+1)/(alpha2/2.0+L))*g2(L);
% end
% %����m1+1��Ȩ��ϵ��



%%
% Coefficient Matix A

A=zeros(m1-1,m1-1);
B=zeros(m2-1,m2-1);

for i=1:m1-1
    for j=1:m2-1

        if(j==i)
            A(i,j)=1+2.*(-1).*r1.*g1(1+1);
        elseif (j==i+1 || j==i-1)
            A(i,j)= (-1).*r1.*(g1(0+1)+g1(2+1));
        elseif (j < i-1)
            A(i,j) = (-1).*r1.*g1(i-j+1+1);
        else
            A(i,j) = (-1).*r1.*g1(j-i+1+1);
        end


        if(j==i)
            B(i,j)=1+2.*(-1).*r2.*g2(1+1);
        elseif (j==i+1 || j==i-1)
            B(i,j)=(-1).*r2.*(g2(0+1)+g2(2+1));
        elseif (j<i-1)
            B(i,j) = (-1).*r2.*g2(i-j+1+1);
        else
            B(i,j) = (-1).*r2.*g2(j-i+1+1);
        end


    end
end

fprintf('Matrix A and B is over\n');      % ����A�������


% %����Toeplitz�����һ��
% t1 =  ./..;l  zeros(1,m1-1);
% t2 = zeros(1,m1-1);
% 
% t1(1) = 1+2.*(-1).*r1.*ge1(1+1);
% t1(2) = (-1).*r1.*(ge1(0+1)+ge1(2+1));
% t1(3:m1-1) = (-1).*r1.*ge1((3:m1-1)+1);
% %һ��Toeplitz����ϵ���ƽ���X����
% 
% t2(1) = 1+2.*(-1).*r2.*ge2(1+1);
% t2(2) = (-1).*r2.*(ge2(0+1)+ge2(2+1));
% t2(3:m1-1) = (-1).*r2.*ge2((3:m1-1)+1);
% %һ��Toeplitz����ϵ���ƽ���Y����


%% �����ʽ��ָ�ʽ

cd('C:\Users\Randoom\Desktop');  % ���õ�ǰĿ¼��current directory
mkdir('Alpha1.7bb');
cd('C:\Users\Randoom\Desktop\Alpha1.7bb'); 
%%mkdir('Alpha_1.7_lamda_5');
%%cd('C:\Users\Randoom\Desktop\Alpha1.7aa\Alpha_1.7_lamda_5');  % ���õ�ǰĿ¼��current directory
%�½�һ���ļ���
u1=zeros(m1-1,m1-1);

tstart=tic;
tol = 1.0e-8;  % the tolerance for the PCG method
it_max = 100; % the maximum number of iterations for the PCG method
for lay = 1:layend   %ʱ�䲽��
    % comute the right-hand side b
    
    for w1=1:m1-1
        b = u(:,w1) + tau * ( u(:,w1) .* ( 1 - u(:,w1) ) .* ( u(:,w1) - a ) - v(:,w1) );
        % ��ָ�ʽ�Ҷ���  ��w1��
        u1(:,w1) = bicg(A,b,tol,it_max);     % �����ݶȷ����Au=b
%        u1(:,w1) = Toeplitz_main(t1,b,tol,it_max);
%         fprintf('X������㣬����ɵڣ�%f ',w1);
%         fprintf('�м���\n');
    end
    for  w2=1:m2-1
        temp=u1(w2,:).';
        %��w2�У���Ϊ������
        u(w2,:) = bicg(B,temp,tol,it_max).';
%       u(w2,:) =Toeplitz_main(t2,temp,tol,it_max).';
%         fprintf('Y������㣬����ɵڣ�%f ',w2);
%         fprintf('�м���\n');
    end
    
    v = v + epsilon*tau*( beta * u -gamma * v -delta );  % ��ⳣ΢�ַ���
    %  ���v
    if(mod(lay,1000)==0)
        fprintf('ʱ���Ϊ:%f\n',lay);
    end
    
    %     %¼����⼰��ر���
    if( mod(lay,100) ==0 ) % ÿ��100��洢һ��
        uu=[xx,yy,reshape(u',(m1-1)*(m1-1),1)];
        [row,col]=size(uu);
        %         str=['data_',num2str(lay),'.plt'];
        folderName{i} = ['Lay', num2str(lay),'.txt'];
        fid=fopen(folderName{i},'w+');
        fprintf(fid,'TITLE="contour of the solution"\n');
        fprintf(fid,'VARIABLES="x","y","u"\n');
        fprintf(fid,'ZONE T="BOX",I= 256,J= 256,F=POINT\n');
        for i=1:row
            for j=1:col
                if j==col
                    fprintf(fid,'%g\n',uu(i,j));
                else
                    fprintf(fid,'%g\t',uu(i,j));
                end
            end
        end
        fclose(fid);
        
        fprintf('ʱ���Ϊ��%f \n',lay);
        
    end
end

tend=toc(tstart)
% 1000ʱ�䲽����ʱ��

end



%%



function u=Toeplitz_main(t,b,tol,it_max)

% The main program for solving Toeplitz systems T\bu = {\bf b}.


n=size(b,1);


% t = kern(n,fchoice); % t is the first column of the Toeplitz matrix T
% b = ones(n,1);     % the right-hand side vector b
pchoice = 1;
[gev,ev] = genevs(t,pchoice); % compute the eigenvalues; see A.3
ig = zeros(n,1); % the initial guess
u = pcg(gev,ev,b,ig,tol,it_max); % call the PCG method

end




%%
function [gev,ev] = genevs(t,pchoice)
% genevs computes the eigenvalues of the circulant matrix in which T
% is embedded and the eigenvalues of circulant preconditioner C.
% t: the first column of the Toeplitz matrix T;
% pchoice: choice of preconditioner;
% gev: the eigenvalues of the circulant matrix in which
% the Toeplitz matrix T is embedded;
%Toeplitz�������ɵ�ѭ���������������
% ev: the eigenvalues of the circulant preconditioner C.

n = length(t);
t1 = conj(t(n:-1:2)); % last column of T
%�������ֵ
gev = real(fft([t 0 t1].'));
%����T���ɵ�ѭ���������������
if pchoice == 1 % T. Chan  preconditioner
    coef = 1/n:1/n:1-1/n;
    ev = fft([t(1) (1-coef).*t(2:n)+coef.*t1])';
elseif pchoice == 2 % Strang  preconditioner.
    ev = fft([t(1:n/2) 0 conj(t(n/2:-1:2))].');
elseif pchoice == 3 % R. Chan  preconditioner
    ev = fft([t(1) t(2:n)+t1].');
elseif pchoice == 4 % Modified Dirichlet kernel
    c = [t(1) t(2:n)+t1].';
    c(2) = c(2)-0.5*t1(1);
    c(n) = c(n)-0.5*t(n);
    ev = fft(c);
elseif pchoice == 5 % de la Vallee Poussin kernel
    c = zeros(1,n);
    c(1) = t(1);
    m = floor(n/2);
    c(2:m+1) = (1:m).*conj(t(2*m:-1:m+1))/m+t(2:m+1);
    c(m+2:2*m) = (m-1:-1:1).*t(m+2:2*m)/m+conj(t(m:-1:2));
    ev = fft(c);
elseif pchoice == 6 % von Hann kernel
    tp = pi/(2*n);
    coef = (cos(tp:tp:(n-1)*tp)).^2;
    c = [t(1) coef.*t(2:n)+(1-coef).*t1];
    ev = fft(c)';
elseif pchoice == 7 % Hamming kernel
    tp = pi/(2*n);
    coef = (cos(tp:tp:(n-1)*tp)).^2;
    c = [t(1) 0.46*(coef.*t(2:n)+(1-coef).*t1)+0.54*(t(2:n)+t1)];
    ev = fft(c)';
elseif pchoice == 8 % Bernstein kernel
    tp = pi/n;
    coef = 0.5*(1+exp((tp:tp:(n-1)*tp)*i));
    c = [t(1) coef.*t(2:n)+(1-coef).*t1];
    ev = fft(c)';
elseif pchoice == 9 % Generalized Jackson kernel: r = 2
    coef = convol(n,2);
    c = [t(1)*coef(1) coef(2:n).*t(2:n)+coef(n:-1:2).*t1];
    ev = fft(c)';
elseif pchoice == 10 % Generalized Jackson kernel: r = 3
    coef = convol(n,3);
    c = [t(1)*coef(1) coef(2:n).*t(2:n)+coef(n:-1:2).*t1];
    ev = fft(c)';
elseif pchoice == 11 % Generalized Jackson kernel: r = 4
    coef = convol(n,4);
    c = [t(1)*coef(1) coef(2:n).*t(2:n)+coef(n:-1:2).*t1];
    ev = fft(c)';
elseif pchoice == 12 % Superoptimal preconditioner
    h = zeros(n,1); % h: first column of the circulant part
    s = zeros(n,1); % s: first column of the skew-circulant part
    h(1) = 0.5*t(1);
    s(1) = h(1);
    t = t.';
    t1 = t1.';
    h(2:n) = 0.5*(t(2:n) + t1);
    s(2:n) = t(2:n)-h(2:n);
    ev1 = fft(h);
    coef = (1:-2/n:2/n-1)';
    c = coef.*s;
    % first column of T. Chan preconditioner
    % for the skew-circulant part
    ev2 = fft(c);
    d = (exp((0:1/n:1-1/n)*pi*i)).';
    s = s.*d;
    sev = fft(s); % eigenvalues of the skew-circulant part
    sev = sev.*conj(sev);
    s = ifft(s);
    s = conj(d).*s;
    h = coef.*s;
    ev3 = fft(h);
    ev = (abs(ev1).^2 + 2*ev1.*ev2+ev3)./ev1;
elseif pchoice == 0
    ev = ones(n,1);
end
ev = real(ev);


end


%%
function coef = convol(n,r)
% convol computes the coefficients for the generalized Jackson kernel.
% K_{m,2r}, see (4.2).
% n: the size of the Toeplitz matrix T;
% r: 2, 3, or 4;
% coef: the coefficients {b_k?(m,r):|k|<=r(m-1)} of
% the generalized Jackson kernel K_{m,2r}.
m = floor(n/r);
a = 1:-1/m:1/m;
r0 = 1;
coef = [a(m:-1:2) a];
while r0 < r
    M = (2*r0+3)*m;
    b1 = zeros(M,1);
    c = zeros(M,1);
    c(1:m) = a;
    c(M:-1:M-m+2) = a(2:m);
    b1(m:m+2*r0*(m-1)) = coef;
    tp = ifft(fft(b1).*fft(c));
    coef = real(tp(1:2*(r0+1)*(m-1)+1));
    r0 = r0+1;
end
M = r*(m-1)+1;
coef = [coef(M:-1:1)' zeros(1,n-M)]';
coef = coef';


end


%%
function u = pcg(gev,ev,b,ig,tol,it_max)
% pcg uses PCG to solve the Toeplitz system Tx=b with circulant
% preconditioner C.
% gev: the eigenvalues of the circulant matrix in which
% the Toeplitz matrix T is embedded, see A.3;
% ev: the eigenvalues of the circulant preconditioner C
% which depend on the choice of preconditioner, see A.3;
% b: the right-hand side vector b;
% ig: the initial guess;
% tol: the tolerance;
% it_max: the maximal number of iterations;
% u: the output solution.
r = b-tx(ig,gev);
u = ig;
aa = norm(r);
t1 = 1;
d = zeros(length(b),1);
iter = 1;
e(1) = aa;
% fprintf('\n at step %1.0f, residual=%e', iter-1, e(iter));
while iter < it_max & e(iter)/e(1) > tol,
    z = cinvx(r,ev);
    t1old = t1;
    t1 = z'*r;
    beta = t1/t1old;
    d = z+beta*d;
    s = tx(d,gev);
    suma = d'*s;
    tau = t1/suma;
    u = u+tau*d;
    r = r-tau*s;
    iter = iter+1;
    e(iter) = sqrt(r'*r);
    % fprintf('\n at step %1.0f, relative residual = %e',...
    % iter-1,e(iter)/e(1));
end
if (iter == it_max),
    fprintf('\n Maximum iterations reached');
end

end



%%
function y = cinvx(d,ev)
% cinvx solves the preconditioning (circulant) system Cy=d.
% d: the right-hand side vector;
% ev: the eigenvalues of the circulant matrix C;
% y: the output vector.
y = ifft(fft(d)./ev);
if norm(imag(y)) < 1.0e-14 % check if v is real
    y = real(y);
end
end


%%

function y = tx(v,gev)
% tx multiplies the Toeplitz matrix T with a vector v.
% v: the vector to be applied;
% gev: the eigenvalues of the circulant matrix
% in which T is embedded, see A.3;
% y: the result of the multiplication.
n = length(v);
y = zeros(2*n,1);
y(1:n) = v;
y = ifft(fft(y).*gev);
y = y(1:n);
if norm(imag(y)) < 1.0e-14 % check if y is real
    y = real(y);
end
end



%%




