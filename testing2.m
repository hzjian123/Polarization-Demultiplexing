%function [x_o,y_o,e_x,e_y,w_xx,w_xy,w_yx,w_yy] ...
 % = LMS(x_i,y_i,u_x,u_y,symbolCounter,w_xx,w_xy,w_yx,w_yy,...
  %trainingSymbolsX,trainingSymbolsY,train_flag,stepSize,nff)
% u_x: delayed line
% u_y: delayed line
[w_xx,w_xy,w_yx,w_yy]=deal(w_xx(:,ind),w_xy(:,ind),w_yx(:,ind),w_yy(:,ind));
[refx,refy]=deal(0);% ref: reference symbol% ref: reference symbol
sigConst=qammod(0:15,16);adaptAfterTraining=1;weightUpdatePeriod=1;sps=2;
% Preallocate output
outlen = (length(x_i) / sps);
x_o = zeros([outlen 1], 'like', x_i);
y_o = zeros([outlen 1], 'like', y_i);
e_x = zeros([outlen 1], 'like', x_i);
e_y = zeros([outlen 1], 'like', y_i);
trainingDelayBufferX=struct;trainingDelayBufferY=struct;
trainingDelayBufferX.Length=0;trainingDelayBufferY.Length=0;
activeDelay=0;
weightUpdateCounter=0;
numTrainSymbols=1;
adaptWeights=1;

for p = 1:outlen


  % fullfill this first sps slots (i.e. a symbol period)

  % count increment
  symbolCounter = symbolCounter + 1;
  weightUpdateCounter = weightUpdateCounter + 1;
  
  % Filter
  x_tmp = w_xx' * u_x + w_xy' * u_y;
  y_tmp = w_yx' * u_x + w_yy' * u_y;

  if symbolCounter <= numTrainSymbols
    if p <= length(trainingSymbolsX)
      dRefx = trainingSymbolsX(p,1);
      dRefy = trainingSymbolsY(p,1);
    else
      % use dummy symbols when running out of trainingSymbols
      if isempty(trainingSymbolsX)
        dummySymX = 1;
        dummySymY = 1;
      else
        dummySymX = trainingSymbolsX(1,1);
        dummySymY = trainingSymbolsY(1,1);
      end
      dRefx = dummySymX;
      dRefy = dummySymY;
    end

    if trainingDelayBufferX.Length > 0
      refx= trainingDelayBufferX.Buffer(trainingDelayBufferX.Pointer);
      trainingDelayBufferX.Buffer(trainingDelayBufferX.Pointer) = dRefx;
      trainingDelayBufferX.Pointer = trainingDelayBufferX.Pointer + 1;
      if trainingDelayBufferX.Pointer > trainingDelayBufferX.Length
        trainingDelayBufferX.Pointer = 1;
      end                                                            
      
      refy= trainingDelayBufferY.Buffer(trainingDelayBufferY.Pointer);
      trainingDelayBufferY.Buffer(trainingDelayBufferY.Pointer) = dRefy;
      trainingDelayBufferY.Pointer = trainingDelayBufferY.Pointer + 1;
      if trainingDelayBufferY.Pointer > trainingDelayBufferY.Length
        trainingDelayBufferY.Pointer = 1;
      end                                                            
    else
      refx = dRefx;
      refy = dRefy;
    end

    if symbolCounter <= activeDelay
      % detect radius
      [~, idxx] = min(abs(sigConst - x_tmp));
       refx= sigConst(idxx);
      [~, idxy] = min(abs(sigConst - y_tmp));
      refy = sigConst(idxy);
    end
  else
    % references come from symbol decision
    [~, idxx] = min(abs(sigConst - x_tmp));
    refx = sigConst(idxx);
    [~, idxy] = min(abs(sigConst - y_tmp));
    refy = sigConst(idxy);
    
    if ~adaptAfterTraining
      adaptWeights = false;
    end
  end

  e_tmp_x = refx - x_tmp;
  e_tmp_y = refy - y_tmp;
 
  if adaptWeights
    % Update tap weights
    if train_flag || weightUpdateCounter == weightUpdatePeriod
      weightUpdateCounter = 0;
      w_xx = w_xx + stepSize * conj(e_tmp_x) * u_x;
      w_xy = w_xy + stepSize * conj(e_tmp_x) * u_y;
      w_yx = w_yx + stepSize * conj(e_tmp_y) * u_x;
      w_yy = w_yy + stepSize * conj(e_tmp_y) * u_y;
    end
  end
  
  % Assign output
  y_o(p,1) = y_tmp;
  x_o(p,1) = x_tmp;
  e_x(p,1) = e_tmp_x;
  e_y(p,1) = e_tmp_y;
end
