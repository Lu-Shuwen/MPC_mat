function matrix = create_matrix(n)
% create the matrix needed in each component of uncertainty set
% n is the horizon of the data
if n ==1
    matrix = [1;-1];
else
    m = create_matrix(n-1);
    a = ones(2^(n-1),1);
    b = -a;
    matrix = [a,m;b,m];
end