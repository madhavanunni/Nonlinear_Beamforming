%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--Get the transmit element index for Plane wave transmit                %
%                                                                        %
%                                                                        %
% Author: Madhavanunni A N                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [txElements, vectors] = getTxElements_PW(bf_points,sysPara,steerAngles)
if nargin<3
    error('Not enough input arguments')
end
% Initializing the Tranmit element matrix
txElements = single(zeros(size(steerAngles,1),size(bf_points,1)));

    for etype = 1:length(steerAngles)
        vec_elem = reshape(repmat(sysPara.elem_pos,[size(bf_points,1) 1]),[size(sysPara.elem_pos,1) size(bf_points,1) 3]);
        vec_bf = permute(reshape(repmat(bf_points,[size(sysPara.elem_pos,1) 1]),[size(bf_points,1) size(sysPara.elem_pos,1)  3]),[2 1 3]);
        vectors.elem_to_bf = vec_elem-vec_bf;
        vectors.rx_distance = sqrt(sum((vectors.elem_to_bf).^2,3));
        if steerAngles(etype)==0
            [vectors.tx_distance(etype,:), txElements(etype,:)] = min(vectors.rx_distance,[],1);
        else
            vectors.tx_distance(etype,:) = bf_points(:,3).*(tand(steerAngles(etype)));
            [txElements(etype,:),~] = min(abs(repmat(vectors.tx_distance(etype,:),sysPara.nElements,1)-vectors.rx_distance),[],2);
        end
    end
if  ~isempty(find(~txElements,1))
  Warning('One or more invalid transmit elements')
end
end

