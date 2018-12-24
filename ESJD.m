function [ answer ] = simulate(  )
global YMIN
global AMIN
global SIGMASTAR
global XSTAR
global SIGMA
global KAPPA
global THETA
global JUMPA
global JUMPB
global LAMBDA0
global LAMBDA1

%defining constants
%we are interested in simulating the SDE
%dr_t=KAPPA(THETA-r_t)dt+SIGMA*sqrt(r_t)dW_t+dJ_t
%where Lambda(x)=LAMBDA0+LAMBDA1*x
%\Delta(x)=x+U(JUMPA,JUMPB)
%
%we will evaluate E[r_T] for T=3)
JUMPA=0.0113;
JUMPB=0.0312;
LAMBDA0=0.011;
LAMBDA1=0.1;
SIGMA=0.013;
KAPPA=0.0117;
THETA=0.0422;
x0=THETA;
XSTAR=x0;
SIGMASTAR=0.25;
T=3;

NN=1e4;

rand('twister', sum(100*clock));
randn('state', sum(100*clock));
format long e

%finding min and max values for phi and A
[YMIN fval]=fminbnd(@(x) phi(x0,x),-2*sqrt(x0)/SIGMA,100);
[AMIN fval]=fminbnd(@(x) -A(x0,x),-2*sqrt(x0)/SIGMA,100);

%calculations
M=1e6;
c=zeros(1,NN);
tic;
for j=1:NN
    at=generateFinal(x0,T); 
	c(1,j)=at;
    if (~mod(j,round(NN/20)))
        disp(['Completed: ',num2str(5*round(j/round(NN/20))),'%']);
    end
end
time_elapsed=toc
mean_c=mean(c)
var_c=var(c)



function output_lambda=Ylambda(x0,y) %Lambda(Finverse(y))
output_lambda=lambda(Finverse(x0,y));

function output_lambda=mYlambda(x0,y) %-Lambda(Finverse(y))
output_lambda=-lambda(Finverse(x0,y));

function varsigma=GenerateNextExitTime(Li) %generate next exit time for standard BM
V=xing();
varsigma=Li^2*V;

function t=phi(x0,x)                       %phi function defined in (14)
t=0.5*(muYprime(x0,x)+muY(x0,x)^2);

function H=GenerateWsigma(Li) %exit value generation
U=rand;
if U<=0.5
    H=-Li;
else
    H=Li;
end

function tau=GenerateJumpCandidates(lambdabar, varsigmanext) %Poisson arrivals generation
%generate exponentially
tau(1)=0;
if lambdabar==0
    tau(2)=varsigmanext;
else
    i=1;
    while (tau(i)<varsigmanext)
        dt=exprnd(1/lambdabar);
        if (tau(i)+dt<varsigmanext)
            tau(i+1)=tau(i)+dt;
        else
            tau(i+1)=varsigmanext;
        end
        i=i+1;
    end
end


function tau=GenerateBrownianThinning(m,M,varsigma) %Poisson arrivals generation for acceptance
tau(1)=0;
i=1;
while (tau(i)<varsigma)
    dt=exprnd(1/(M-m));
    if (tau(i)+dt<varsigma)
        tau(i+1)=tau(i)+dt;
    else
        tau(i+1)=varsigma;
    end
    i=i+1;
end

function b=GenerateBBTriplets(varsigma,tau) %generate a triplet of Brownian Bridges used in meander construction
for j=1:3
    b(1,j)=0;
    b(length(tau),j)=0;
end
z=randn(3,length(tau));
for i=2:length(tau)-1
    for j=1:3
        b(i,j)=(varsigma-tau(i))*b(i-1,j)/(varsigma-tau(i-1))+z(j,i)*sqrt((varsigma-tau(i))*(tau(i)-tau(i-1))/(varsigma-tau(i-1)));
    end
end

function b=GenerateBB(varsigma,tau,Triplets) %construct meander from triplets
for i=1:length(tau)
	b(i)=sqrt(((varsigma-tau(i))/varsigma+Triplets(i,1))^2+(Triplets(i,2))^2+(Triplets(i,3))^2);
end

function t=SigmaNormalize(tau,L) %scaling meander to standard BM with exit boundary L
t=zeros(1,length(tau));
for i=1:length(tau)
    t(i)=tau(i)/L^2;
end

function I=TestI(U,tau1,varsigma,BB) %meander testing function
dif=0;
Accepted=true;
Resolved=false;
q=0;
n=0;
x=BB(length(BB)-1);
s=tau1(length(tau1))-tau1(length(tau1)-1);
while ~Resolved
    qprev=q;
    n=n+1;
    if (mod(n,2))
        m=ceil(n/2);
        q=q-(1/x)*(4*m-x)*exp(-4*m*(2*m-x)/s);
    else
        m=n/2;
        q=q+(1/x)*(4*m+x)*exp(-4*m*(2*m+x)/s);
    end
    if n>=log(4)*(s)/8*1/x+2
        if ((qprev<=q)&&(q+1<U(length(U))))
            Resolved=true;
            Accepted=false;
        else
            if ((qprev>=q)&&(q>U(length(U))-1))
                Resolved=true;
                Accepted=Accepted&&true;
            end
        end
    end
end
for i=1:length(U)-1
    t=tau1(length(tau1))-tau1(length(tau1)-i-1);
    s=tau1(length(tau1))-tau1(length(tau1)-i);
    y=BB(length(BB)-i-1);
    x=BB(length(BB)-i);
    p1=1/(1-exp(-2*x*y/(t-s)));
	p=p1;
    n=0;
    Resolved=false;
    while ~Resolved
        pprev=p;
        n=n+1;
        if (mod(n,2))
            m=ceil(n/2);
            p=p-p1*(exp(-2*(2*m-x)*(2*m-y)/(t-s))+exp(-2*(2*(m-1)+x)*(2*(m-1)+y)/(t-s)));
        else
            m=n/2;
            p=p+p1*(exp(-2*m*(4*m+2*(x-y))/(t-s))+exp(-2*m*(4*m+2*(x-y))/(t-s)));
        end
        if n>=log(3)*(t-s)/8*max(1/x,1/y)+1
            if ((pprev<=p)&&(p<U(i)))
                Resolved=true;
                Accepted=false;
            else
                if ((pprev>=p)&&(p>U(i)))
                    Resolved=true;
                    Accepted=Accepted&&true;
                end
            end
        end
    end
    if ~Accepted
        break;
    end
end
I=Accepted;

function [t ind]=MergeTimes(t1,t2) %merges two times' arrays
t=zeros(1,length(t1)+length(t2)-2);
ind=t;
i=2;
j=2;
t(1)=t1(1);
t(length(t))=t1(length(t1));
for k=2:length(t)
    if t1(i)<=t2(j)
        t(k)=t1(i);
        i=i+1;
    else
        t(k)=t2(j);        
        j=j+1;
        ind(k)=1;
    end
end  

function [I0 j]=TestAll(B,t,x0,x,M,m,hf,ind,lbar,mA,Tb,V,W) %main candidate testing function
j=0;
I0=true;
J=false;
I=rand(1,length(B));
for i=2:length(t)-1-hf
    if (ind(i)==1) %B(i) is a real jump candidate
        if (I(i)<(Ylambda(x0,x+B(i)))/lbar)
            I0=I0&&(V<exp(-m*t(i))/max(1,exp(-m*Tb)));
            if I0              
                I0=(W<exp(A(x0,x+B(i)))/mA);
            end
            j=i;    %jump at point i
            J=true; %jump happened
            break;
        end
    else
        if (ind(i)==0) %B(i) is a helper arrival candidate
            if (I(i)<((phi(x0,x+B(i))-m)/(M-m)))
                I0=false;
                break;
            end        
        end
    end
end
if (~J)&&I0 %if no jump happened, accept the sample
    I0=I0&&(V<exp(-m*t(length(t)-hf))/max(1,exp(-m*Tb)));
    if I0              
        I0=(W<exp(A(x0,x+B(length(B)-hf)))/mA);
    end
end

function X = xing()
%this function copied from Burq and Jones paper
accepted = 0;
while ~accepted
	X = gamrnd(1.088870, 0.810570);
	Y = rand*1.243707*gampdf(X, 1.088870, 0.810570);
	sqrt2piX3 = sqrt(2*pi*X^3);
    N = max([ceil(0.275*X), 3]);
	K = (1+2*[-N:N]);
	fN0 = sum((-1).^[-N:N].*K.*exp(-K.^2./(2*X)))/sqrt2piX3;
    N = N + 1;
    fN1 = fN0 + (-1)^N*((1-2*N)*exp(-(1-2*N)^2/(2*X))+(1+2*N)*exp(-(1+2*N)^2/(2*X)))/sqrt2piX3;
    while sign((Y - fN0)*(Y - fN1)) == -1
        fN0 = fN1;
        N = N + 1;
        fN1 = fN0 + (-1)^N*((1-2*N)*exp(-(1-2*N)^2/(2*X)) + (1+2*N)*exp(-(1+2*N)^2/(2*X)))/sqrt2piX3;
    end
    if Y <= fN1
        accepted = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HERE USER-DEFINABLE FUNCTIONS BEGIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THE A(x) FUNCTION DEFINED IN EQN (12)
%EDIT THE FUNCTION ACCORDING TO YOUR MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a=A(x0,x)
global SIGMA
global KAPPA
global THETA
a=(4*KAPPA*THETA-SIGMA^2)/(2*SIGMA^2)*(log(x+2*sqrt(x0)/SIGMA)-log(2*sqrt(x0)/SIGMA))-KAPPA/4*((x+2*sqrt(x0)/SIGMA)^2-(2*sqrt(x0)/SIGMA)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THE muY(x) FUNCTION DEFINED IN EQN (7)
%EDIT THE FUNCTION ACCORDING TO YOUR MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a=muY(x0,x)
global SIGMA
global KAPPA
global THETA
a=(4*KAPPA*THETA-SIGMA^2)/(2*SIGMA^2)/(x+2*sqrt(x0)/SIGMA)-KAPPA/2*(x+2*sqrt(x0)/SIGMA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THE muY(x) FUNCTION DERIVATIVE
%EDIT THE FUNCTION ACCORDING TO YOUR MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a=muYprime(x0,x)
global SIGMA
global KAPPA
global THETA
a=-(4*KAPPA*THETA-SIGMA^2)/(2*SIGMA^2)/(x+2*sqrt(x0)/SIGMA)^2-KAPPA/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THE \LAMBDA(x) INTENSITY FUNCTION
%EDIT THE FUNCTION ACCORDING TO YOUR MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_lambda=lambda(x)
global LAMBDA0
global LAMBDA1
output_lambda=LAMBDA0+LAMBDA1*x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THE \Delta_Y FUNCTION DEFINED IN (9)
%EDIT THE FUNCTION ACCORDING TO YOUR MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_c=Ycjump(x0,y)
global JUMPA
global JUMPB
R=(JUMPB-JUMPA)*rand+JUMPA;
output_c=F(x0,Finverse(x0,y)+R)-y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THE LAMPERTI TRANSFORM DEFINED IN (4)
%EDIT THE FUNCTION ACCORDING TO YOUR MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=F(x0,x)
global SIGMA
f=2*(sqrt(x)-sqrt(x0))/SIGMA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THE INVERSE LAMPERTI TRANSFORM
%EDIT THE FUNCTION ACCORDING TO YOUR MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=Finverse(x0,y)
global SIGMA
f=(y*SIGMA/2+sqrt(x0))^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THE OPTIMAL BOUNDARY CHOICE FUNCTION
%EDIT THE FUNCTION ACCORDING TO YOUR CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L=GetLi(y,ybar)
if (y<-25)
    L=-1.25/(ybar+25)*y+ybar*1.25/(ybar+25);
elseif (y<-15)
    L=0.1*y+3.75;
elseif (y<-8)
    L=2.25;    
elseif (y<-4)
    L=0.05*y+2.95;
elseif (y<2)
    L=2.75;
elseif (y<6)
    L=-0.05*y+2.85;
elseif (y<18)
    L=2.25;
else
    L=1+50/(y+22);
end

function f=generateFinal(x0,T)
%define all global constants required
global SIGMA
global YMIN
global AMIN
currentt=0; %current time
x=x0;       %starting x0
ybar=-2*sqrt(x0)/SIGMA; %only needed for AJD
while currentt<T
    restart=true;
    y=F(x0,x);
    Li=GetLi(y,ybar); %selecting new Li based on current Xt and bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THE OPTIMAL m, M, mA, h, H, \bar{\lambda} SELECTION
%EDIT THE FUNCTION ACCORDING TO THE FORM OF phi, A, lambda(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if ((y-Li-YMIN)*(y+Li-YMIN)<0)
        m=phi(x0,YMIN);
        M=phi(x0,y+Li);
        MM=phi(x0,y-Li);
        if (M<MM)
            M=MM;
        end
    else
        if (y+Li<YMIN)
            M=phi(x0,y-Li);
            m=phi(x0,y+Li);
        else
            M=phi(x0,y+Li);
            m=phi(x0,y-Li);
        end
    end

    if ((y-Li-AMIN)*(y+Li-AMIN)<0)
        mA=exp(A(x0,AMIN));
    else
        if (y+Li<AMIN)
            mA=exp(A(x0,y+Li));
        else
            mA=exp(A(x0,y-Li));
        end
    end
    lambdabar=lambda(Finverse(x0,y+Li)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    while restart
        restart=false;
        varsigma=GenerateNextExitTime(Li); %generating exit time with xing
        tau=GenerateJumpCandidates(lambdabar, min(varsigma,T-currentt));
        %(now we have T, currentt+tau, currentt+varsigma). we need the
        %lowest to play around with it
        H1=GenerateWsigma(Li); %generating W_varsigma candidate
        V=rand;
        W=rand;
        if (T-currentt<varsigma)
            %if we hit the time horizon
            vsigma=T-currentt;            
            I=false;
            while ~I %until we encounter a qualified meander, we keep on generating it
                tau0=GenerateBrownianThinning(m,M,vsigma); %generating Poisson arrivals for ARM on DIFFUSION (not JUMP)                
                [tau0 taui]=MergeTimes(tau0,tau);
                tau0(length(tau0)+1)=varsigma;
                taui(length(taui)+1)=0;
                tau1=SigmaNormalize(tau0,Li); %scaling the times to fit usual meander at hitting time 1
                
                %OLD VERSION
                BB=3;
                while max(BB)>2
                    B=GenerateBBTriplets(varsigma/Li^2,tau1);
                    BB=GenerateBB(varsigma/Li^2,tau1,B); %generating meander from triplets                    
                end
                %now BB is W_{\varsigma_1-t_1},W_{\varsigma_1-t_2} etc
                U=rand(1,length(tau1)-1);
                I=TestI(U,tau1,varsigma/Li^2,BB); %testing meander for not hitting the UPPER Li
                if I %if it qualifies, scale it back
                    if (H1<0)
                        for i=1:length(BB)
                            BB(i)=Li*(-1+BB(i));
                        end
                    else
                        for i=1:length(BB)
                            BB(i)=Li*(1-BB(i));
                        end                            
                    end
                end
            end
            [I2 j]=TestAll(BB,tau0,x0,y,M,m,1,taui,lambdabar,mA,T-currentt,V,W);
            if I2
                if j>0
                    x=Finverse(x0,y+BB(j)+Ycjump(x0,y+BB(j)));              
                    currentt=currentt+tau0(j);
                    restart=false;                    
                else
                    x=Finverse(x0,y+BB(length(BB)-1));  
                    currentt=T;
                    restart=false;
                end
            else
                restart=true;
            end
        else
                %tic;
                %HERE WE JUST GENERATE EXIT TIME AND VALUE
                tau0=GenerateBrownianThinning(m,M,varsigma);
                [tau0 taui]=MergeTimes(tau0,tau);
                if length(tau0)>2
                    tau1=SigmaNormalize(tau0,Li);
                    I=false;
                    while ~I
                        %OLD VERSION
                        BB=3;
                        while max(BB)>2
                            B=GenerateBBTriplets(varsigma/Li^2,tau1);
                            BB=GenerateBB(varsigma/Li^2,tau1,B); %generating meander from triplets                    
                        end                        
                        
                        U=rand(1,length(tau1)-1);
                        I=TestI(U,tau1,varsigma/Li^2,BB);
                        if I
                            if (H1<0)
                                for i=1:length(BB)
                                    BB(i)=Li*(-1+BB(i));
                                end
                            else
                                for i=1:length(BB)
                                    BB(i)=Li*(1-BB(i));
                                end
                            end
                        end
                    end                  
                else
                    BB=[0 H1];
                end   
                [I1 j]=TestAll(BB,tau0,x0,y,M,m,0,taui,lambdabar,mA,T-currentt,V,W); 
                if I1                  
                    if j>0
                        x=Finverse(x0,y+BB(j)+Ycjump(x0,y+BB(j)));
                        currentt=currentt+tau0(j);
                        restart=false;                    
                    else
                        x=Finverse(x0,y+H1);
                        currentt=currentt+varsigma;
                        restart=false;
                    end
                else
                    restart=true;
                end                
        end 
    end
end
f=x;
