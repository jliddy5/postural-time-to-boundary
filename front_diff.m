function [ xp ] = front_diff( x, dt )
%FRONT_DIFF Numerical approximation of the derivative of x via a forward
%difference which provides first order accuracy (i.e., truncation errors on
%the order of dt).
%
% ARGUMENTS
% x - Vector or matrix containing numerical data. If x is a matrix
% differentation operates along the vertical dimensions (i.e., each column
% is differentiated).
% 
% dt - Unit change in time between samples (assumed constant).
%
% RETURNS
% xp - First-order derivative of x with respect to time. If x is an m x n
% vector or matrix, xp will have a size of m-1 x n.
%
%=========================================================================%

[ r, c ] = size(x); %determine the size of x

xp = zeros(r-1, c); %create empty matrix for the derivative values

%for each column
for i = 1:c
    
    %for each row
    for j = 1:r-1
        
        xp(j,i) = (x(j+1,i) - x(j,i)) / dt; %calculate forward difference
        
    end %end row
    
end %end column

end %end front_diff

