function G = NRDgradientFunction(A)

global NUM_LM;
global NUM_TL_POINTS;

global LM_TL_IDX;
global NN_VERT_TL;

global TL_POINTS;
global LM_SC;
global LM_TL;

global NN_SC;
global CONF_SC_POINTS;

global TL_HEAD_IDX;
global TL_HAND_IDX;
global TL_ALL_IDX;

global NUM_TL_HEAD;
global NUM_TL_HAND;
global NUM_TL_ALL;

global IS_VALID_NN;

global W_DATA;
global W_SMOOTH_HEAD;
global W_SMOOTH_HANDS;
global W_SMOOTH_GEN;
global W_LM;

A = reshape(A,[4 4 NUM_TL_POINTS]);

Grad_LM = zeros(4,4,NUM_TL_POINTS);

for i = 1:NUM_LM
    Grad_LM(:,:,LM_TL_IDX(i,1)) = W_LM   * 2*((A(:,:,LM_TL_IDX(i,1))*LM_TL(:,:,i)) - LM_SC(:,:,i))*LM_TL(:,:,i)';
end

G_SMOOTH = zeros(4,4,NUM_TL_POINTS);

nSmooth = 0;
for i = 1:NUM_TL_HEAD
    NNvertSMPoint = NN_VERT_TL{TL_HEAD_IDX(1,i),1};
    numNNvert = size(NNvertSMPoint,2);
    for j = 1:numNNvert
        G_SMOOTH(:,:,TL_HEAD_IDX(1,i)) = G_SMOOTH(:,:,TL_HEAD_IDX(1,i)) + 2*(A(:,:,TL_HEAD_IDX(1,i)) - A(:,:,NNvertSMPoint(j)));
    end
    G_SMOOTH(:,:,TL_HEAD_IDX(1,i)) = W_SMOOTH_HEAD * G_SMOOTH(:,:,TL_HEAD_IDX(1,i));
    nSmooth = nSmooth + numNNvert;
end

for i = 1:NUM_TL_HAND
    NNvertSMPoint = NN_VERT_TL{TL_HAND_IDX(1,i),1};
    numNNvert = size(NNvertSMPoint,2);
    for j = 1:numNNvert
        G_SMOOTH(:,:,TL_HAND_IDX(1,i)) = G_SMOOTH(:,:,TL_HAND_IDX(1,i)) + 2*(A(:,:,TL_HAND_IDX(1,i)) - A(:,:,NNvertSMPoint(j)));
    end
    G_SMOOTH(:,:,TL_HAND_IDX(1,i)) = W_SMOOTH_HANDS * G_SMOOTH(:,:,TL_HAND_IDX(1,i));
    nSmooth = nSmooth + numNNvert;
end

for i = 1:NUM_TL_ALL
    NNvertSMPoint = NN_VERT_TL{TL_ALL_IDX(1,i),1};
    numNNvert = size(NNvertSMPoint,2);
    for j = 1:numNNvert
        G_SMOOTH(:,:,TL_ALL_IDX(1,i)) = G_SMOOTH(:,:,TL_ALL_IDX(1,i)) + 2*(A(:,:,TL_ALL_IDX(1,i)) - A(:,:,NNvertSMPoint(j)));
    end
    G_SMOOTH(:,:,TL_ALL_IDX(1,i)) = W_SMOOTH_GEN * G_SMOOTH(:,:,TL_ALL_IDX(1,i));
    nSmooth = nSmooth + numNNvert;
end

G_DATA = zeros(4,4,NUM_TL_POINTS);

for i=1:NUM_TL_POINTS
    G_DATA(:,:,i) = W_DATA * IS_VALID_NN(i,1) * CONF_SC_POINTS(i,1)*2*(A(:,:,i)*TL_POINTS(:,:,i) - [NN_SC(i,:) 1]')*TL_POINTS(:,:,i)';
end

G_LM = 0;
if (NUM_LM > 0)
    G_LM = (1/NUM_LM * Grad_LM);
end
G = G_LM + (1/nSmooth * G_SMOOTH) + (1/sum(IS_VALID_NN) * G_DATA);
G = reshape(G,[4*4*NUM_TL_POINTS 1]);
end
