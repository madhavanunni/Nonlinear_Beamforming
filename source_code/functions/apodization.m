%-- Function which assigns different apodization to a set of elements
%-- Function prototype:
%-- [apod_weights_left,apod_weights_right] = apodization(sysPara,bf_points,txElements,rxElements,angleInd,window)
%-- Windows: 'singleElement', 'gaussian8', 'gaussian12', 'box'
%-- Authors: Madhavanunni A N

%-- Date:

function [apod_weights_left,apod_weights_right] = apodization(sysPara,~,bf_points,txElements,rxElements,angleInd,window)

nElements = sysPara.nElements;
elem_idx = 1:nElements;
elemIdx = (ones(size(bf_points,1),1)*elem_idx).';
elemPosX = sysPara.elem_pos(:,1);
elemPos = (sysPara.elem_pos(:,3)*ones(1,size(bf_points,1)));

switch(window)
    case 'singleElement'
        rxLeftElemPos = elemPosX(rxElements.left_element_idx(:,:,angleInd))*ones(1,nElements);
        rxRightElemPos = elemPosX(rxElements.right_element_idx(:,:,angleInd))*ones(1,nElements);
        apod_weights_left = single(elemPos==rxLeftElemPos.');
        apod_weights_right = single(elemPos==rxRightElemPos.');
    case 'gaussian8'
        txElemPos = (elemPosX(txElements)*ones(1,nElements)).';
        rxLeft = ones(nElements,1)*rxElements.left_element_idx(:,:,angleInd);
        rxRight = ones(nElements,1)*rxElements.right_element_idx(:,:,angleInd);
        leftApert = double(elemPosX<txElemPos);
        rightApert = double(elemPosX>txElemPos);
        apod_weights_left = single(exp(-(((elemIdx-rxLeft).^2)/(2*16))).*leftApert);
        apod_weights_right = single(exp(-(((elemIdx-rxRight).^2)/(2*16))).*rightApert);
    case 'gaussian12'
        txElemPos = (elemPosX(txElements)*ones(1,nElements)).';
        rxLeft = ones(nElements,1)*rxElements.left_element_idx(:,:,angleInd);
        rxRight = ones(nElements,1)*rxElements.right_element_idx(:,:,angleInd);
        leftApert = double(elemPosX<txElemPos);
        rightApert = double(elemPosX>txElemPos);
        apod_weights_left = single(exp(-(((elemIdx-rxLeft).^2)/(2*24))).*leftApert);
        apod_weights_right = single(exp(-(((elemIdx-rxRight).^2)/(2*24))).*rightApert);
    case 'box'
        txElemPos = (elemPosX(txElements)*ones(1,nElements)).';
        leftApert = double(elemPosX<txElemPos);
        rightApert = double(elemPosX>txElemPos);
        apod_weights_left = ones(size(leftApert)).*leftApert;
        apod_weights_right = ones(size(rightApert)).*rightApert;
        
    otherwise
        error('Unknown window type. Known types are: boxcar, Gaussian, single element');
end

end
