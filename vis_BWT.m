function vis_BWT (string)

l = size(string, 2);
matrix(l, l) = 0;
% 1 = A
% 2 = G
% 3 = T
occ(1,3) = 0;

for i = 1:l
	switch string(i)
		case 'A'
			occ(1,1) = occ(1,1) + 1;
		case 'G'
			occ(1,2) = occ(1,2) + 1;
		case 'T'
			occ(1,3) = occ(1,3) + 1;
	end
	x = 1;
	for j = [i:l 1:i - 1] 	
		matrix(i, x) = string(j);
		x = x + 1;
	end
end
char(matrix)
M = char(sortrows(matrix))

[M(230, [1:l])];
[M(231, [1:l])];
bwt(l) = 0;

for i = 1:l
	bwt(i) = M(i, l);
end

bwt;
occ;

% xlswrite('matrix0.xlsx', M);

end