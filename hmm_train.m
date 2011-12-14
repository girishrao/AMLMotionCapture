function [] = hmm_train()
%Girish Rao
%hmm_train() - no input parameters required
%
%hmm_train trains both the walking and running sequences
%In order to avoid having to train every time the user wants to test,
%the parameters generated are saved to a workspace .MAT file
%The trained walking params are stored in VARW.MAT and the trained running
%params are stored in VARR.MAT (these are loaded during subsequent testing)
%This program calls hmm(train,i) to execute the actual JTA EM steps
%
%NOTE: The trainr variable provided in the original hw1.mat file was
%conflicting with one of MatLab's predefined "trainr" functions. Thus, I
%loaded this variable into another .MAT file renamed as trainrun.mat 

  if (nargin ~= 0)              %check correct number of arguments
    help hmm_train
  else

    load hw1.mat;
    load trainrun.mat;
    tfile = 'varW.mat';
    tfile2 = 'varR.mat';
    
    %Train Walking
    sprintf('%s', 'Training 4 Walking HMMs...')
    [transW1 priorsW1 meansW1] = hmm(trainw,1);
    [transW2 priorsW2 meansW2] = hmm(trainw,2);
    [transW3 priorsW3 meansW3] = hmm(trainw,3);
    [transW4 priorsW4 meansW4] = hmm(trainw,4);

    sumtrans = transW1 + transW2 + transW3 + transW4;
    sumtrans = sumtrans / 4;
    sumRow = repmat(sum(sumtrans,2), 1, 2);
    sumtransW = sumtrans ./ sumRow;

    sumpriors = priorsW1 + priorsW2 + priorsW3 + priorsW4;
    sumpriors = sumpriors / 4;
    sumRow = sum(sumpriors,2);
    sumpriorsW = sumpriors ./ sumRow;

    summeans = meansW1 + meansW2 + meansW3 + meansW4;
    summeansW = summeans / 4;
    summeansW;

    %Save Parameters to Workspace
    save(tfile, 'sumtransW', 'sumpriorsW', 'summeansW');
    %Animate for Mean Pose
    animate( transpose(summeansW(:,:)) );
        
    %Train Running
    sprintf('%s', 'Training 4 Running HMMs...')
    [transR1 priorsR1 meansR1] = hmm(trainrun,1);
    [transR2 priorsR2 meansR2] = hmm(trainrun,2);
    [transR3 priorsR3 meansR3] = hmm(trainrun,3);
    [transR4 priorsR4 meansR4] = hmm(trainrun,4);

    sumtrans = transR1 + transR2 + transR3 + transR4;
    sumtrans = sumtrans / 4;
    sumRow = repmat(sum(sumtrans,2), 1, 2);
    sumtransR = sumtrans ./ sumRow;

    sumpriors = priorsR1 + priorsR2 + priorsR3 + priorsR4;
    sumpriors = sumpriors / 4;
    sumRow = sum(sumpriors,2);
    sumpriorsR = sumpriors ./ sumRow;

    summeans = meansR1 + meansR2 + meansR3 + meansR4;
    summeansR = summeans / 4;
    summeansR;

    %Save Parameters to Workspace
    save(tfile2, 'sumtransR', 'sumpriorsR', 'summeansR');
    %Animate Mean Pose
    animate( transpose(summeansR(:,:)) );
  end
  
end