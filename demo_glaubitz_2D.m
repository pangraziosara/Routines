
function demo_glaubitz_2D

%--------------------------------------------------------------------------
% OBJECT
%--------------------------------------------------------------------------
% Demo used to compute integrals from samples involving scattered data.
%
% 1. Determines domain and functions to study from function defined by the
%    user.
% 2. computes integral using only evaluations at scattered data:
%
% (a) For methods based on evaluating in scattered data, using the values
%    to evaluate at points of an algebraic rule:
%    Makes experiments, evaluating the integrand at cubature nodes of an
%    algebraic rule, when only samples at scattered data are available.
%
%    Important: algebraic degree of exactness of the cubature rule is fixed.
%
% (b) By Glaubitz algorithm computes the value of the integral, with a 
%    formula based on scattered data.
%
% 3. Makes plot of the domain, scattered samples, cubature nodes.
%--------------------------------------------------------------------------
% Modified: November 6, 2024.
%--------------------------------------------------------------------------
% COPYRIGHT
%--------------------------------------------------------------------------
% Copyright (C) 2024-
%
% Authors:
% Alvise Sommariva.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------



% ...................... Problem settings  ..........................

%--------------------------------------------------------------------------
% 1:'polygon'; 2:'disk'; 3:'lune'; 4:'circular-annular-sector';
% 5:'sector'; 6:'asymmetric-circular-sector'; 7:'asymmetric-annulus';
% 8:'vertical-circular-zone'; 9:'horizontal-circular-zone';
% 10:'circular-segment'; 11:'symmetric-lens'; 12:'butterfly';
% 13:'candy'; 14:'NURBS'; 15:'union-disks';16:'asymmetric-circular-sector';
% 17:'square'; 18:'unit-square[0,1]x[0,1]'; 19:'triangle'; 20:'polygcirc';
% 21, 'unit-simplex';
%--------------------------------------------------------------------------
% Note: use deg_rule <=15 for NURBS
%--------------------------------------------------------------------------
tic
domain_type=11;

deg_rule=40; %mettere 20/30 % degree of precision of the rule; choose deg=10 for NURBS.
% grado formula di riferimento (basata su gqellblend)
card=100;    % scattered data cardinality
scat_type='halton'; % scattered data type.

f=@(x,y) 1./( (1+x.^2).*(1+y.^2)); 
%f =@(x,y) exp(-x.^2-y.^2);

method='Glaubitz-algorithm';


%nella figura punti rossi ref, neri glaubitz

% ....................... Preliminary settings  ...........................

% some additional routines are required
addpath(genpath('../EXTERNAL_ROUTINES/'));

% flag_compression: if "1" compression will be attempted, if necessary,
% otherwise the code will propose an alternative rule.
flag_compression=1;



% ...................... Making numerical tests  ..........................


% 1. Define domain
%    Note: the domain routine is in "../EXTERNAL_ROUTINES/DCUB" folder

fprintf('\n \t * defining integration domain');
domain_example=domain_str(domain_type);
domain_struct=define_domain(domain_example);
fprintf('\n \n \t \t The domain is: '); disp(domain_example);

% 1A. Define scattered data
fprintf('\n \t * defining scattered data');
[t,dbox,area_domain]=define_scattered_pointset(card,domain_struct,...
    scat_type);
ft=feval(f,t(:,1),t(:,2));

% 1B. Define cubature rule
%    Note: the domain routine is in "../EXTERNAL_ROUTINES/DCUB" folder
XYW=define_cub_rule(domain_struct,deg_rule,flag_compression);

% 1C: Compute integral from evaluation at scattered data


area_dbox=diff(dbox(1,:))*diff(dbox(2,:));
wQMC=(area_domain/area_dbox)/size(t,1);
[w,deg_rule]=glaubitz_algorithm(t,domain_struct,wQMC,'l1'); % chage approach to change formula
I=w'*ft;


%--------------------------------------------------------------------------
% 2. Display statistics.
%--------------------------------------------------------------------------
fprintf('\n \t -------------------------------------------------------');
fprintf('\n \t Domain (integration)          : '); disp(domain_struct.domain);
fprintf('\n \t Scattered data type           : '); disp(scat_type);
fprintf(   '\t Evaluation method             : '); disp(method);
fprintf('\n \t Scatt. points card. required  : %5.0f',card);
fprintf('\n \t Scatt. points card. provided  : %5.0f \n',size(t,1));

fprintf('\n \t Degree of precision GL          : %5.0f',deg_rule);
% grado di precisione formula di Glaubitz
% caso bivariato sarebbero grado 6 --> 28 punti
card_rule=size(XYW,1); % cardinality della formula di riferimento
wp = find(w>0);  %cardinality della formula di Glaubitz
fprintf('\n \t Cardinality of the formula Ref   : %5.0f',card_rule);
fprintf('\n \t Cardinality of the formula GL   : %5.0f',length(wp));
dim_poly=(deg_rule+1)*(deg_rule+2)/2;
fprintf('\n \t Dimension of the poly. space  : %5.0f',dim_poly);


fprintf('\n \t Integral approximation        : %1.15e',I);

switch domain_type
    case 17
        atol=10^(-12); rtol=10^(-12);
        IR = integral2(f,-1,1,-1,1,'AbsTol',atol,'RelTol',rtol);
    case 18
        atol=10^(-12); rtol=10^(-12);
        IR = integral2(f,0,1,0,1,'AbsTol',atol,'RelTol',rtol);
    otherwise
        ft=feval(f,XYW(:,1),XYW(:,2));
        IR=(XYW(:,3))'*ft;
end

fprintf('\n \t Integral computed by alg.rule : %1.15e',IR);
fprintf('\n \t Absolute error                : %1.1e',abs(IR-I));
fprintf('\n \t Relative error                : %1.1e',abs(IR-I)/abs(I));

fprintf('\n \n');

fprintf('\n \t ........................................................ ');

fstr = strrep(char(f),'@(x)','');
fprintf('\n \n \t The integrand is: '); disp(fstr);
fprintf('\n \t The domain is   : '); disp(domain_struct.domain);
fprintf('\n \n \t ........................................................ ');
fprintf('\n \n');



% 3. plot domain and scattered data
figure(1)
plot_2D(domain_struct,t,XYW(:,1:2));
saveas(gcf,domain_example,'epsc');




toc




function domain_example=domain_str(example)

%--------------------------------------------------------------------------
% OBJECT
%--------------------------------------------------------------------------
% Given a number "example", it provides the name of the corrsponding domain
% in the gallery,
%--------------------------------------------------------------------------

switch example

    case 1, domain_example='polygon';
    case 2, domain_example='disk';
    case 3, domain_example='lune';
    case 4, domain_example='circular-annular-sector';
    case 5, domain_example='sector';
    case 6, domain_example='asymmetric-circular-sector';
    case 7, domain_example='asymmetric-annulus';
    case 8, domain_example='vertical-circular-zone';
    case 9, domain_example='horizontal-circular-zone';
    case 10, domain_example='circular-segment';
    case 11, domain_example='symmetric-lens';
    case 12, domain_example='butterfly';
    case 13, domain_example='candy';
    case 14, domain_example='NURBS';
    case 15, domain_example='union-disks';
    case 16, domain_example='asymmetric-circular-sector';
    case 17, domain_example='square';
    case 18, domain_example='unit-square[0,1]x[0,1]';
    case 19, domain_example='triangle';
    case 20, domain_example='polygcirc';
    case 21, domain_example='unit-simplex';
end












