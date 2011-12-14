function [a, piMatrix, eta] = hmm(train,i)
%Girish Rao
%Usage: hmm(train, i)
%where train is the entire training cell array and i is the index (1-4) of
%the sequence in the cell array to be trained on
%

states = 2;
sequences = 4;

  if (nargin ~= 2) % check correct number of arguments
    help hmm
  else
      M1 = train{i};
      [yDim, T1] = size(M1);
      
      %Init data structures
      piMatrix = (1/states) * ones(1, states);  %1*2
      a = abs(rand(states, states));            %2*2
      eta = abs(rand(states, yDim));            %2*123
      sigma = abs(rand(1, states));             %1*2
      emissions = zeros(states, T1);            %2*358

      
      %Normalize transition matrix a; sum each row (the second dimension)
      %Get log probabilities
      sums = sum(a, 2);
      a = (a ./ (repmat(sums, 1, states)));

      %Init variables
      loop = 1;
      iterationLL = zeros(1000);
      iteration = 0;
      covar = eye(yDim) * 0.25;
      zeta = log(ones(states, T1-1));
      phi = log(ones(states, T1-1));
      p = (2*pi) ^ (yDim/2);
      
      while(loop == 1)
          
          data = transpose(M1(:,1));
          for i=1:states
              mean = eta(i,:);
              phi(i,1) = log(1/( p * sqrt(det(covar)))) + (-(1/2) * (data - mean) * pinv(covar) * (data - mean)');
          end
          for j=1:(T1-1)
              data = transpose(M1(:,j+1));
              for i=1:states
                  mean = eta(i,:);
                  zeta(i,j) = log(1/( p * sqrt(det(covar)))) + (-(1/2) * (data - mean) * pinv(covar) * (data - mean)');
              end
          end


          %%%E STEP%%%
          a = log(a);
          piMatrix = log(piMatrix);

          if(iteration >= 1)
              for j=2:(T1-1)
                  phi(:,j) = log(phi(:,j));
              end
          end
          phi(:,1) = phi(:,1) + transpose(piMatrix);              
          
          psiQQ = zeros(states, states, T1-1);
          psiQQ = repmat(a, [1 1 (T1-1)]);  %states*states*T-1
          
          %Collect
          for i = 1:(T1-1)
              
              psiQQ(1,1,i) = psiQQ(1,1,i) + phi(1,i) + zeta(1,i);
              psiQQ(1,2,i) = psiQQ(1,2,i) + phi(1,i) + zeta(2,i);
              psiQQ(2,1,i) = psiQQ(2,1,i) + phi(2,i) + zeta(1,i);
              psiQQ(2,2,i) = psiQQ(2,2,i) + phi(2,i) + zeta(2,i);
                                     
              if (i < (T1-1))
                  sumCols = numericalTrickCol(psiQQ, i);
                  sumCols = sumCols - transpose(phi(:,i+1));
                  phi(:,i+1) = transpose(sumCols);
                                    
              end
          end

          
          %Distribute
          old_phi = phi;
          for i = (T1-1):-1:1
              sumCols = numericalTrickCol(psiQQ, i);
              %sumCols = sumCols - transpose( (zeta(:,i)) );
              zeta(:,i) = transpose(sumCols);
              
              sumRows = numericalTrickRow(psiQQ, i);
              %sumRows = sumRows - phi(:,i);
                 phi(:,i) = sumRows;
              
              if(i > 1)
                  sumCols = transpose(phi(:,i) - old_phi(:,i));
                  psiQQ(1,1,i-1) = psiQQ(1,1,i-1) + sumCols(1,1);
                  psiQQ(1,2,i-1) = psiQQ(1,2,i-1) + sumCols(1,1);
                  psiQQ(2,1,i-1) = psiQQ(2,1,i-1) + sumCols(1,2);
                  psiQQ(2,2,i-1) = psiQQ(2,2,i-1) + sumCols(1,2);
                  
              end
          end
          
          LL = psiQQ(1,1,T1-1) + psiQQ(1,2,T1-1) + psiQQ(2,1,T1-1) + psiQQ(2,2,T1-1);
          iteration = iteration + 1;
          iterationLL(iteration) = LL;
          %sprintf('%e %d', LL, iteration)
          if(iteration >= 2)
              if( abs( (iterationLL(iteration) - iterationLL(iteration-1)) / iterationLL(iteration-1) ) < 1e-06)
                 %sprintf('%s', '*******************')
                 loop = 0;
              end
          end

          %Normalize
          for i=1:(T1-1)
              sumRows = numericalTrickRow(psiQQ, i);
              psiQQ(1,1,i) = psiQQ(1,1,i) - sumRows(1,1);
              psiQQ(1,2,i) = psiQQ(1,2,i) - sumRows(1,1);
              psiQQ(2,1,i) = psiQQ(2,1,i) - sumRows(2,1);
              psiQQ(2,2,i) = psiQQ(2,2,i) - sumRows(2,1);
              
              sumRows2 = sum(phi(:,i));
              phi(1,i) = phi(1,i) - sumRows2;
              phi(2,i) = phi(2,i) - sumRows2;
              
              sumRows2 = sum(zeta(:,i));
              zeta(1,i) = zeta(1,i) - sumRows2;
              zeta(2,i) = zeta(2,i) - sumRows2;
          end
          
          %Exponentiate
          piMatrix = exp(piMatrix);
          psiQQ = exp(psiQQ);
          a = exp(a);

          %%%M STEP%%%%%
          %update pi
          piMatrix(1,1) = phi(1,1);
          piMatrix(1,2) = phi(2,1);
          
          %update a
          aSum = 0;
          denom = 0;
          for i=1:states
              for j=1:states
                  for t=1:(T1-1)
                      aSum = aSum + psiQQ(i,j,t);
                  end
                  a(i,j) = aSum;
                  for k=1:states
                      for m=1:(T1-1)
                          denom = denom + psiQQ(i,k,m);
                      end
                  end
                  a(i,j) = a(i,j) / denom;
                  aSum = 0;
                  denom = 0;
              end
          end
         
                          
          %update mu
          muSum = zeros(123,1);
          probSum = 0;
          for i=1:states
              for j=1:(T1-1)
                  muSum = muSum + (zeta(i,j) * M1(:,j+1));
                  probSum = probSum + zeta(i,j);
              end
              eta(i,:) = transpose(muSum) / probSum;
              muSum = zeros(123,1);
              probSum = 0;
          end

      end %end While
      a;
      piMatrix;
      eta;
      %animate(transpose(eta(:,:)));

  end %end else if
end %end hmm function


%matrix parameter is 2*2
function [sumCols] = numericalTrickCol(matrix, tp)

   sumCols = zeros(1,2);
   [maxVal, index] = max(matrix(:,1,tp));
   if(index == 1)
       ind = 2;
   elseif(index == 2)
       ind = 1;
   end                  
   sumCols(1,1) = maxVal + log( exp(matrix(ind,1,tp) - maxVal) + exp(matrix(index,1,tp) - maxVal) );

   [maxVal, index] = max(matrix(:,2,tp));
   if(index == 1)
       ind = 2;
   elseif(index == 2)
       ind = 1;
   end
   sumCols(1,2) = maxVal + log( exp(matrix(ind,2,tp) - maxVal) + exp(matrix(index,2,tp) - maxVal) );
   sumCols;
end

%matrix parameter is 2*2
function [sumRows] = numericalTrickRow(matrix, tp)

   sumRows = zeros(2,1);
   [maxVal, index] = max(matrix(1,:,tp));
   if(index == 1)
       ind = 2;
   elseif(index == 2)
       ind = 1;
   end                  
   sumRows(1,1) = maxVal + log( exp(matrix(1,ind,tp) - maxVal) + exp(matrix(1,index,tp) - maxVal) );

   [maxVal, index] = max(matrix(:,2,tp));
   if(index == 1)
       ind = 2;
   elseif(index == 2)
       ind = 1;
   end
   %sprintf('%s %d %e', '********', tp, matrix(2,:,tp))

   sumRows(2,1) = maxVal + log( exp(matrix(2,ind,tp) - maxVal) + exp(matrix(2,index,tp) - maxVal) );
   sumRows;
end
