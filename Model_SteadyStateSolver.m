% Atri Model with Mechanics - Steady State Solver V1

 % REMEMBER: change TH(c), if applicable 

function R = Model_SteadyStateSolver(p,lm)

K1 = 4.88566;
K2 = 0.1;
K3 = 0.05;
Ve = 1e-6;
kinfty = 52e-6;
g = 0.5;



syms c h TH(c)

%--- TH(c) tuning parameters --- 

deg   = 5;

%TH(c) = 10*c^2/(1+10*c^2);
%TH(c) = 1/(1+10*(c-4)^2);

%TH(c) = 10*c/(1+10*c);

TH(c) = 0;

%-------------------------------
%{
deg = 5; % can be set after seeing number of roots produced by solC

syms rs(c) de(c) pl(c)

M = 0.2;
n = 4; % only for even n
W = 2;

rs(c) = (c^n/(M^n+c^n));
de(c) =  1-(c^n/(M^n+c^n))*heaviside(-c);
pl(c) = rs(c) - de(c-W) ;

TH(c) = pl(c);
%}
%-------------------------------

digits(16);
eqnC = ((((Ve)/(g*(kinfty*(((p)^4)/(1 + (p)^4)))))^4)*(c^6)) + ((K2^2)*(((Ve)/(g*(kinfty*(((p)^4)/(1 + (p)^4)))))^4)*(c^4)) + ((1 - K1)*(c^2)) + ((K2^2)) - ((K1*(K3^2))) == 0;
solC = vpasolve(eqnC, c)

% THE FOLLOWING CODE RETURNS ONLY ONE REAL ROOT
%% --- Root Selector ---

rt = zeros([1 deg]);
ct = 0;
pv = 0;
rn = size(solC,1);
% symbolic -> double conversion causes loss of precision & rounding errors
% use pv to refer to symbolic value corresponding to extracted root
for j=1:rn
 z = double(solC(j)); 
 % convert symbolic value to double ie explicit typecast  
 % boolean operators cant work with symbolic values
 REC = log10(abs(real(z)));
 IMC = log10(abs(imag(z)));
 
 if (isreal(z)) %deals with purely real(+ve and -ve) roots
     if(z >= 0)
         ct = ct + 1;
         rt(ct) = z;
         if (ct == 1)
             pv = j;     
         end  
     end    
     
 elseif (real(z) > 0) %deals with complex roots having +ve real part
     if((REC-IMC) > 3)
         ct = ct + 1;
         rt(ct) = real(z);
         if (ct == 1)
             pv = j;     
         end  
     end 
 end
 
end

R(1) = vpa(rt(1),16); % R(1) = steady state 'c', R(2) = steady state 'h' 
R(2) = 1/(1 + ((((Ve)/(g*(kinfty*(((p)^4)/(1 + (p)^4)))))*(R(1))^4)));
R(3) = TH(R(1))/1;
% we only use the first real root, should be enough for Atri Model
% avoid multiple root selection to reduce complexity
                     
if (ct == 0)
    Out = "No suitable real roots"
    R(1) = 0;
    R(2) = 0;
    R(3) = 0;
else
    Out = "Real root extracted"
    % Out = "Root list:"
    % solC
end

end