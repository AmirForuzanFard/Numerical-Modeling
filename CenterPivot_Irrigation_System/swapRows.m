function v = swapRows(v,i,j)
% Swap rows i and j of column vector or matrix v.
% USAGE: v = swapRows(v,i,j)

temp = v(i,:);
v(i,:) = v(j,:);
v(j,:) = temp;
