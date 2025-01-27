function write_array_to_file(array, filename)
    % Write elements of an array of arbitrary dimensions to a text file.
    % Each element is written on a new line.
    %
    % Args:
    %     array: Input array of arbitrary dimensions.
    %     filename: Name of the output text file.

    % Open the file for writing
    filename = ['/Users/dpn/', filename, '.m.txt']
    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not open the file for writing.');
    end

    % Flatten the array and write each element to the file
    tmpArray = array;

    numDims = ndims(tmpArray);
    arraySize = size(tmpArray);
    sizeStr = mat2str(arraySize);
    if numDims == 2
        tmpArray = transpose(tmpArray);
    elseif numDims == 3
        tmpArray = permute(tmpArray, [2 1 3]);
        #fprintf('skipping...\n');
    end
    fprintf(fid, '(');
    fprintf(fid, '%d,', arraySize);
    fprintf(fid, ')\n');

    tmpArray = tmpArray(:); % Flatten the array into a column vector
    cnt = 0;
    for i = 1:numel(tmpArray)
        fprintf(fid, '%12.10e\n', tmpArray(i)); % Write each element on a new line
        cnt = cnt+1;
        if cnt > 1000
            break;
        end
    end

    % Close the file
    fclose(fid);
end
