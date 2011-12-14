function [] = hmm_main(test, i)
%Girish Rao
%Usage: hmm_main(test, i)
%where test is an entire cell array and i is the index (1-4) of
%the sequence to be tested
%

  if (nargin ~= 2)              %check args
    help hmm_main
  else

    load varW.mat;
    load varR.mat;
    load hw1.mat;
    load trainrun.mat;

    sprintf('%s', 'Testing on HMMs...')
    LLw = hmm_test(test{i},sumtransW, sumpriorsW, summeansW);
    LLr = hmm_test(test{i}, sumtransR, sumpriorsR, summeansR);

    sprintf('%s %e', 'Walking LL', LLw)
    sprintf('%s %e', 'Running LL', LLr)
    
  end
end