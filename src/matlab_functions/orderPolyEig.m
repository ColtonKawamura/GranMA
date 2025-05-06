function [outData] = orderPolyEig(inData)
% orderPolyEig - Sort eigenvalues and eigenvectors from polyeig

    outData = inData;

    for k=1:size(inData.pressure, 1)
        [~, idx] = sort(abs(imag(outData.eigenValues{k})));
        outData.eigenValues{k} = outData.eigenValues{k}(idx);
        outData.eigenVectors{k} = outData.eigenVectors{k}(:, idx);
    end

end