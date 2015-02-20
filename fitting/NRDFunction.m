function E = NRDFunction(A)

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

global E_TOTAL;
global E_DATA;
global E_SMOOTH;
global E_LM;

A = reshape(A,[4 4 NUM_TL_POINTS]);
A_LM = cat(3,A(:,:,LM_TL_IDX));
A_DATA = A;

E_LM_Gen = 0;
for i = 1:NUM_LM
    E_LM_Gen = E_LM_Gen  + (norm( A_LM(:,:,i) *LM_TL(:,:,i) -  LM_SC(:,:,i)))^2;
end

E_SMOOTH_Head = 0;
E_SMOOTH_Hands = 0;
E_SMOOTH_Gen = 0;

nSmoothAdd = 0;
for i = 1:NUM_TL_HEAD;
    NNvertSMPoint = NN_VERT_TL{TL_HEAD_IDX(1,i),1};
    numNNvert = size(NNvertSMPoint,2);
    for j = 1:numNNvert
        E_SMOOTH_Head = E_SMOOTH_Head + norm(A(:,:,TL_HEAD_IDX(1,i)) - A(:,:,NNvertSMPoint(j)),'fro')^2;
    end
    nSmoothAdd = nSmoothAdd + numNNvert;
end

for i = 1:NUM_TL_HAND;
    NNvertSMPoint = NN_VERT_TL{TL_HAND_IDX(1,i),1};
    numNNvert = size(NNvertSMPoint,2);
    for j = 1:numNNvert
        E_SMOOTH_Hands = E_SMOOTH_Hands + norm(A(:,:,TL_HAND_IDX(1,i)) - A(:,:,NNvertSMPoint(j)),'fro')^2;
    end
    nSmoothAdd = nSmoothAdd + numNNvert;
end

for i = 1:NUM_TL_ALL;
    NNvertSMPoint = NN_VERT_TL{TL_ALL_IDX(1,i),1};
    numNNvert = size(NNvertSMPoint,2);
    for j = 1:numNNvert
        E_SMOOTH_Gen = E_SMOOTH_Gen + norm(A(:,:,TL_ALL_IDX(1,i)) - A(:,:,NNvertSMPoint(j)),'fro')^2;
    end
    nSmoothAdd = nSmoothAdd + numNNvert;
end

E_data_total = 0;

for i=1:NUM_TL_POINTS
    E_data_total = E_data_total + (IS_VALID_NN(i,1)*CONF_SC_POINTS(i,1)*(norm(A_DATA(:,:,i)*TL_POINTS(:,:,i) - [NN_SC(i,:) 1]'))^2);
end

E_DATA   = 1/sum(IS_VALID_NN) * (W_DATA * E_data_total);
E_SMOOTH = 1/nSmoothAdd * (W_SMOOTH_HEAD * E_SMOOTH_Head + W_SMOOTH_HANDS * E_SMOOTH_Hands + W_SMOOTH_GEN * E_SMOOTH_Gen);
E_LM = 0;
if (NUM_LM > 0)
    E_LM = 1/NUM_LM * W_LM * E_LM_Gen;
end
E_TOTAL = E_DATA + E_SMOOTH + E_LM;

E = E_TOTAL;
end
