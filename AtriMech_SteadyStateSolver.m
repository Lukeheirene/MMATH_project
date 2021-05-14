% Atri Model - Steady State Solver V1



function R = AtriMech_SteadyStateSolver(mu,lm)

kf = 16.2;
th = 2;
k1 = 0.7;
y  = 2;
ky = 0.1;
k2 = 0.7;
b  = 0.111;
% ko = 1;

K1 = (kf*th)/k1;
T  = (y*th)/k1;
K  = ky/k1;
;

syms c h TH(c)

 

deg   = 5;



TH(c) = 0;



digits(16);
eqnC = mu*K1*(1/(1+c^2))*((b+c)/(1+c))-(T*c/(K+c))+lm*TH == 0;
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
R(2) = 1/(1+R(1)^2);
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