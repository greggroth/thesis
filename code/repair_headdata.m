function [ head_out ] = repair_headdata( head_in )
% repaid_headdata repopulates the headdata matrix
%   If any changes are made to the index column in the headdata 
% matrix, use this function to repopulate and correct the 
% parameter values before running any other functions.
%  head_in can be either 3 or 4 dimenisions


% =====================
% = Parameter Storage =
% =====================

QmSTORE = [0 0 26.1 11600 0 26.1 697 0 0 302 15575 0 500 1100 5192];
cSTORE = [1006 4600 2110 3640 3800 1300 3720 3000 4200 2300 3680 3500 3010 3150 3600];
rhoSTORE = [1.3 1057 1080 1035.5 1007 1850 1126 1076 1009 916 1035.5 1151 978.5 1100 1027.4];
kSTORE = [0.026 0.51 0.65 0.534 0.5 0.65 0.527 0.4 0.594 0.25 0.565 0.4975 0.3738 .342 .503];
wSTORE = [0 1000 3 45.2 0 1.35 40 0 0 2.8 67.1 3.8 3.3 12 23.7];

if ndims(head_in)==4
    head_in = head_in(:,:,:,1);
end

% Reassign the parameter values
head_out = cat(4,head_in, zeros(size(head_in)), QmSTORE(head_in), cSTORE(head_in), rhoSTORE(head_in), kSTORE(head_in), wSTORE(head_in));

end