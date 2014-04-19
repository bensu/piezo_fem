function index = index_range(step,counter)
% index = index_range(step,counter)
% Yields a useful sequence
    index = (step*(counter-1)+1):(step*counter);
end