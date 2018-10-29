function [ifr] = i_free(o1,o2,opt)


if opt == 0
    ifr = o1*2.545727780163160 - o2*1.545727780163160;
elseif opt == 1
    ifr = o1*2.53125 - o2*1.53125;
elseif opt == 2
    ifr = o1*2.260604327518826 - o2*1.260604327518826;
elseif opt == 3
    ifr = o1*2.487168313616925 - o2*1.487168313616925;
end

end

