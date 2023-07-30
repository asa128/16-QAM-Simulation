function outVector = randbin(n,logic_0,logic_1)
%RANDBIN Generates a vector of n numbers, each pseudo-randomly assigned to either the value of logic_0 or logic_1

outVector = logic_0 + (logic_1 - logic_0).*randi([0 1], 1, n);

end