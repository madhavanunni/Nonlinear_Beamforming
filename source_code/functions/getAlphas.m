%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alpha estimation                                                       %
%                                                                        %
% Estimates all possible transmit to recieve angles using Law of Cosines %
% Also estimates the alpha value pairs where alpha left and alpha right  %
% is equal.                                                              %
%                                                                        %
% Author: Madhavanunni A N                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alphas, vectors, rxElements, txRxAngles] = getAlphas(sysPara,bf_points,alphaSet,txElements,vectors)

n = size(sysPara.elem_pos,1);
nAlphas = size(alphaSet,2);
txRxAngles.alpha_left = single(NaN(sysPara.no_emission_types,size(bf_points,1),nAlphas));
txRxAngles.alpha_right = single(NaN(sysPara.no_emission_types,size(bf_points,1),nAlphas));

rxElements.left_element_idx = single(NaN(sysPara.no_emission_types,size(bf_points,1),nAlphas));
rxElements.right_element_idx = single(NaN(sysPara.no_emission_types,size(bf_points,1),nAlphas));
alphas = single(NaN(sysPara.no_emission_types,n,size(bf_points,1)));

for etype=1:sysPara.no_emission_types
    
    for bfp=1:size(bf_points,1)
        tx_idx = txElements(etype,bfp);
        for i=1:n
            vectors.tx_to_rx(etype,i) = abs(sysPara.elem_pos(tx_idx)-sysPara.elem_pos(i));
            alphas(etype,i,bfp) = acos(((vectors.rx_distance(i,bfp))^2+(vectors.tx_distance(etype,bfp))^2-(vectors.tx_to_rx(etype,i))^2)/(2*vectors.rx_distance(i,bfp)*vectors.tx_distance(etype,bfp)));
        end
        
        for deg=1:nAlphas
            refAlpha = single(alphaSet(deg)*pi/180);
            [diff,left_ind] = min(abs(refAlpha-alphas(etype,1:tx_idx-1,bfp)));
            if (diff*180/pi) < 03
                alpha_left = alphas(etype,left_ind,bfp);
                [diff, right_ind] = min(abs(alpha_left-alphas(etype,tx_idx+1:sysPara.nElements,bfp)));
                right_ind=right_ind+tx_idx;
                if (diff*180/pi) < 03
                    valid_alphas(etype,bfp,deg)=refAlpha;
                    txRxAngles.alpha_left(etype,bfp,deg) = alphas(etype,left_ind,bfp);
                    rxElements.left_element_idx(etype,bfp,deg)=left_ind;
                    txRxAngles.alpha_right(etype,bfp,deg) = alphas(etype,right_ind,bfp);
                    rxElements.right_element_idx(etype,bfp,deg) = right_ind;
%                 else
%                     pause()
                end
%             else
%                     pause()
            end
        end
    end
    
end
end
