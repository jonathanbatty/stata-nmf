function [A,X]= nmf ghals (Y,opts)

% Non−negative Matrix Factorization (NMF) with Hierarchical
% Alternating Least−Squares (ALS) algorithms
%
% INPUT
% Y : source with size of I x T
% opts : structure of optional parameters for algorithm (see defoptions)
% .tol: tolerance of stopping criteria (Frobenius distance) (1e−5)
% .J: number of components
% .algtype vector (1x2) indicates which algorithm will be used ([6 6])
% .init: vector 1 x 2 defines initialization types ([1 1])
% (see nmf initialize )
% .nrestart: number of multi initializations (1)
% .niter: number of iterations (300)
% .niniter: number of iterations for updating only A or only X (1)
% .ellnorm: ell p norm
% .Anorm and Xnorm: normalize A, X
% .lsparse , .lcorr , .lsmth: vector (1x2) indicates degrees of
% sparsness , orthogonality , smoothness for each factor A, X
% OUTPUT
% A and X with Y = A ∗ X.
%
% Copyright 2008 by A. Cichocki and A.H. Phan and R. Zdunek
% Optimized by Anh Huy Phan and Andrzej Cichocki − 01/2008
%% =======================================================================

Y(Y< 0) = eps; [I,T]= size(Y);

% Set algorithm parameters from input or by using defaults
defopts = struct('tol',1e−5,'J',I,'algtype ',[1 4],'init ',[1 1],'nrestart ',1,'niter ',300,'niniter ',1,'ellnorm ',2,'Xnorm ',0,'Anorm ',1,'alpha ',1,'beta ',1,'lsparse ',[0 0],'lcorr ',[0 0],'lsmooth ',[0 0],'tau',50,'alpha0 ',20,'A0',[],'X0',[],'verbose ',0);

if ∼exist('opts ','var'), opts = struct; end

[tol,J,algtype ,init ,nrestart ,niter ,niniter ,ellnorm ,Xnorm ,Anorm , alpha ,beta ,lspar ,lcorr ,lsmth ,tau,alpha0 ,A0,X0,verbose] = scanparam(defopts ,opts);

No_iter = 30; % number of iterations to select the best trial

for nr= 0: nrestart
pause(.0001)
if (nr == nrestart )&&( nrestart > 0) % initialize
A = A best; X = X best; No iter = niter;
else
[A,X] = nmf initialize (Y,J,init ,{A0 X0 '}); X = X';
end
cost = costfunction;
% Compute the error
Res = Y − A∗ X; Ainew = zeros(I,1); Xinew = zeros(1,T);
49
for k = 1: No iter
if (algtype (1) == 2) && (algtype (2) == 2) % alpha−HALS−1
for j = 1:J
Res = Res + [−Ainew A(:,j)] ∗ [Xinew; X(j,:)];
Resalpha = phi(Res,alpha);
55
Xinew = max(eps,phi((A(:,j)'∗Resalpha − lspar (2)) ...
./sum(A(:,j).^( alpha+1)),1/alpha));
Ainew = max(eps,phi(Resalpha∗X(j,:)fl-lspar(1),1/alpha));
if ellnorm / = 0
Ainew = Ainew/norm(Ainew ,ellnorm);
end
A(:,j) = Ainew; X(j,:) = Xinew;
end
elseif (algtype (1) == 4) && (algtype (2) == 4) % beta−HALS−1
for j = 1:J
Res = Res + [ −Ainew A(:,j)] ∗ [Xinew ;X(j,:)];
Xinew = max(eps,((A(:,j).^beta)'∗Res−lspar (2)) ...
/sum(A(:,j).^( beta +1)));
Ainew = max(eps,(Res∗(X(j,:).^beta)'−lspar (1)));
if ellnorm / = 0
Ainew = Ainew/norm(Ainew ,ellnorm);
end
A(:,j) = Ainew; X(j,:) = Xinew;
end
else
% Update for A
for t = 1: niniter % inner iterations
A=nmfcore(Y',X',A',algtype(1),lspar(1),lcorr(1),lsmth (1))';
end % t for A
if (algtype (2) / = 0) && Anorm
A = normalize(A);
end

% Update for X
for t = 1: niniter % inner iterations
X = nmfcore(Y,A,X,algtype(2),lspar(2),lcorr(2),lsmth (2));
end
if (algtype (1) / = 0) && Xnorm
X = normalize(X')';
end
end

if (nr == nrestart) && (mod(k,30)==0) % stopping condition
checkstoppingcondition
if verbose
fprintf(1,'Best trial %d, step %d, Cost value %d\n',...
nr best+1,k,cost);
end
if stop , break; end
end
end % k
if (nr < nrestart) % select best trial
cost = costfunction;
if (nr == 0) | | (cost < cost min )
A best = A; X best = X; cost min = cost; nr best = nr;
end
if verbose
fprintf(1, 'Trial %d, Cost value = %e\n',nr+1, cost);
end
end
end % nr
112
function X = nmfcore(Y,A,X,type alg , lspar ,lcorr ,lsmth)
Xsm = lsmth ∗ [X(:,2) (X(:,1:end−2)+X(:,3:end ))/2 X(:,end−1)];
switch type alg
case 1 % Fast HALS
AA = A'∗A;
normA = diag (1./( diag(AA)+ lsmth+lcorr));
AA = normA∗AA;
AAX = normA∗(A'∗Y) − AA ∗ X + Xsm;
for j = 1:J
Xn= max(eps, X(j,:) + AAX(j,:) −lspar −lcorr∗sum(X));
AAX = AAX − AA(:,j) ∗ (Xn − X(j,:)) ;
X(j,:) = Xn;
end

case {2,3} % Alpha−HALS
Res = Y−A∗X;
Ainew = zeros(size(Y,1),1);
Xinew = zeros(1,size(Y,2));
for j = 1:J
Res = Y−A∗X + A(:,j)∗X(j,:);
Resalpha = real(Res.^alpha);
Xinew = max(eps,phi((A(:,j)'∗Resalpha − lspar ...
−lcorr∗sum(X)+ Xsm(j,:))...
./(sum(A(:,j).^( alpha +1))+ lsmth),1/alpha));
X(j,:) = Xinew;
Ainew = A(:,j);
end

case {4,5} % Fast Beta−HALS
normA = diag (1./(sum(A.^( beta +1))+ lsmth + lcorr));
Abeta = real(A.^beta);
AY = normA∗Abeta '∗Y;
AA = normA ∗Abeta '∗A;
for j = 1:J
X(j,:) = max(eps, X(j,:) + AY(j,:) − AA(j,:) ∗ X − ...
lspar − lcorr∗sum(X) + Xsm(j,:));
end
150
case 6 % ALS
X = max(eps,pinv(A'∗A)∗A'∗Y);

case 7 % Regularized ALS
alpha reg = alpha0∗exp(−k/tau);
if cond(A) > 1E6 , lambdaX = 1E−6; else lambdaX = 0; end
X = max(eps,pinv(A'∗A + alpha reg + lambdaX∗I)∗A'∗Y);
end
end

function y = phi(x,alpha)
y = real(x.^alpha);
end

function A = normalize(A)
A = bsxfun(@rdivide ,A,sum(A.^ellnorm).^(1/ ellnorm));
end

function checkstoppingcondition
cost old = cost; cost = costfunction;
stop = abs(cost old−cost) ≤ tol∗cost;
end

function cost = costfunction
Yhat = A∗X+eps;
cost = sum(sum(Y.∗log(Y./Yhat + eps) − Y + Yhat));
end
end