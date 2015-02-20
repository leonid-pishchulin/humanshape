function [regParams,Bfit,ErrorStats]=absor(A,B,varargin)
%ABSOR - a tool for solving the absolute orientation problem using Horn's 
%quaternion-based method, that is, for finding the rotation, translation, and 
%optionally also the scaling, that best maps one collection of point coordinates
%to another in a least squares sense. The function works for both 2D and 3D
%coordinates, and also gives the option of weighting the coordinates non-uniformly.
%The code avoids for-loops to maximize speed.
%
%DESCRIPTION:
%
%As input data, one has
%
%  A: a 2xN or 3xN matrix whos columns are the coordinates of N source points.
%  B: a 2xN or 3xN matrix whos columns are the coordinates of N target points.
%
%The syntax
%
%     [regParams,Bfit,ErrorStats]=absor(A,B) 
%
%solves the unweighted/unscaled registration problem
%
%           min. sum_i ||R*A(:,i) + t - B(:,i)||^2 
%
%for unknown rotation matrix R and unknown translation vector t. 
%
%This is a  special case of the more general problem
%
%           min. sum_i w(i)*||s*R*A(:,i) + t - B(:,i)||^2   
%
%where s>=0 is an unknown global scale factor to be estimated along with R and t
%and w is a user-supplied length N vector of  weights. One can include either
%s or w or both in the problem formulation using the syntax,
%
%  [regParams,Bfit,ErrorStats]=absor(A,B,'param1',value1,'param2',value2,...)
%
%with parameter/value pair options
%
%  'doScale' -  Boolean flag. If TRUE, the global scale factor, s, is included.
%               Default=FALSE.
%
%  'weights' - the length N-vector of weights, w. Default, no weighting.
%
%
%
%out:
%
%
% regParams: structure output with estimated registration parameters,
%
%     regParams.R:   The estimated rotation matrix, R
%     regParams.t:   The estimated translation vector, t
%     regParams.s:   The estimated scale factor (set to 1 if doScale=false).
%     regParams.M:   Homogenous coordinate transform matrix [s*R,t;[0 0 ... 1]].
%
%     For 3D problems, the structure includes
%
%        regParams.q:   A unit quaternion [q0 qx qy qz] corresponding to R and
%                       signed to satisfy max(q)=max(abs(q))>0
%
%     For 2D problems, it includes
%
%        regParams.theta: the counter-clockwise rotation angle about the
%                         2D origin
%
%
%  Bfit: The rotation, translation, and scaling (as applicable) of A that 
%        best matches B.
%
%
% ErrorStats: structure output with error statistics. In particular,
%             defining err(i)=sqrt(w(i))*norm( Bfit(:,i)-B(:,i) ), 
%             it contains
%
%      ErrorStats.errlsq = norm(err) 
%      ErrorStats.errmax = max(err) 
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Matt Jacobson
% Copyright, Xoran Technologies, Inc.  http://www.xorantech.com


%%Input option processing and set up

 options.doScale = 0;     
 options.weights = []; 
    
for ii=1:2:length(varargin)
    param=varargin{ii};
    val=varargin{ii+1};
    if strcmpi(param,'doScale'), 
        options.doScale=val;
    elseif strcmpi(param,'weights')
        options.weights=val;
    else
       error(['Option ''' param ''' not recognized']);
    end
end

 doScale = options.doScale;     
 weights = options.weights;             



if ~isempty(which('bsxfun'))
 matmvec=@(M,v) bsxfun(@minus,M,v); %matrix-minus-vector
 mattvec=@(M,v) bsxfun(@times,M,v); %matrix-minus-vector
else
 matmvec=@matmvecHandle;
 mattvec=@mattvecHandle;
end


dimension=size(A,1);

if dimension~=size(B,1),
    error 'The number of points to be registered must be the same'
end


%%Centering/weighting of input data


if isempty(weights)
    
 sumwts=1;   
    
 lc=mean(A,2);  rc=mean(B,2);  %Centroids
 left  = matmvec(A,lc); %Center coordinates at centroids
 right = matmvec(B,rc); 
 
else
    
  sumwts=sum(weights);
  
  weights=full(weights)/sumwts;
     weights=weights(:);
  sqrtwts=sqrt(weights.');
  
  
  lc=A*weights;   rc=B*weights; %weighted centroids
  
  left  = matmvec(A,lc);  
     left  = mattvec(left,sqrtwts);
  right = matmvec(B,rc);
    right = mattvec(right,sqrtwts);  
  
end

M=left*right.';


%%Compute rotation matrix

switch dimension

    case 2
    
        
        Nxx=M(1)+M(4); Nyx=M(3)-M(2);

        N=[Nxx   Nyx;...
           Nyx   -Nxx];

        [V,D]=eig(N);

        [trash,emax]=max(real(  diag(D)  )); emax=emax(1);

        q=V(:,emax); %Gets eigenvector corresponding to maximum eigenvalue
        q=real(q);   %Get rid of imaginary part caused by numerical error


        q=q*sign(q(2)+(q(2)>=0)); %Sign ambiguity
        q=q./norm(q);

        R11=q(1)^2-q(2)^2;
        R21=prod(q)*2;

        R=[R11 -R21;R21 R11]; %map to orthogonal matrix


        
        
    case 3
        
        [Sxx,Syx,Szx,  Sxy,Syy,Szy,   Sxz,Syz,Szz]=dealr(M(:));

        N=[(Sxx+Syy+Szz)  (Syz-Szy)      (Szx-Sxz)      (Sxy-Syx);...
           (Syz-Szy)      (Sxx-Syy-Szz)  (Sxy+Syx)      (Szx+Sxz);...
           (Szx-Sxz)      (Sxy+Syx)     (-Sxx+Syy-Szz)  (Syz+Szy);...
           (Sxy-Syx)      (Szx+Sxz)      (Syz+Szy)      (-Sxx-Syy+Szz)];

        [V,D]=eig(N);

        [trash,emax]=max(real(  diag(D)  )); emax=emax(1);

        q=V(:,emax); %Gets eigenvector corresponding to maximum eigenvalue
        q=real(q);   %Get rid of imaginary part caused by numerical error

        [trash,ii]=max(abs(q)); sgn=sign(q(ii(1)));
         q=q*sgn; %Sign ambiguity

        %map to orthogonal matrix
        
        quat=q(:);
        nrm=norm(quat);
        if ~nrm
         'Quaternion distribution is 0'    
        end

        quat=quat./norm(quat);

        q0=quat(1);
        qx=quat(2); 
        qy=quat(3); 
        qz=quat(4);
        v =quat(2:4);


        Z=[q0 -qz qy;...
           qz q0 -qx;...
          -qy qx  q0 ];

        R=v*v.' + Z^2;

    otherwise
        error 'Points must be either 2D or 3D'

end

%%

if doScale
   
     summ = @(M) sum(M(:));
  
     sss=summ( right.*(R*left))/summ(left.^2);
     t=rc-R*(lc*sss);
     
     
else
    
    sss=1;
    t=rc-R*lc;
   
 
end


regParams.R=R;
regParams.t=t;
regParams.s=sss;


if dimension==2
    
    regParams.M=[sss*R,t;[0 0 1]];
    regParams.theta=atan2(q(2),q(1))*360/pi;
    
else%dimension=3
    
    regParams.M=[sss*R,t;[0 0 0 1]];
    regParams.q=q/norm(q);
    
end

if nargout>1
    
    Bfit=matmvec((sss*R)*A,-t);
    
end

if nargout>2
    
    l2norm = @(M,dim) sqrt(sum(M.^2,dim));

    err=l2norm(Bfit-B,1);
    
    if ~isempty(weights), err=err.*sqrtwts; end

    ErrorStats.errlsq=norm(err)*sqrt(sumwts); %unnormalize the weights
    ErrorStats.errmax=max(err);
    
   
end
    



function M=matmvecHandle(M,v)
%Matrix-minus-vector

    for ii=1:size(M,1)
     M(ii,:)=M(ii,:)-v(ii);
    end
 
function M=mattvecHandle(M,v)
%Matrix-times-vector

    for ii=1:size(M,1)
     M(ii,:)=M(ii,:).*v;
    end
    
function varargout=dealr(v)

  varargout=num2cell(v);     

