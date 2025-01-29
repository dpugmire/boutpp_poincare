function write_array_to_file(array, filename)
    % Write elements of an array of arbitrary dimensions to a text file.
    % Each element is written on a new line.
    %
    % Args:
    %     array: Input array of arbitrary dimensions.
    %     filename: Name of the output text file.

    % Open the file for writing
    filename = ['/Users/dpn/', filename, '.m.txt'];
    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not open the file for writing.');
    end

    maxCnt = 5000;
    maxCnt = -1;
    numDims = ndims(array);
    arraySize = size(array);
    sizeStr = mat2str(arraySize);

    cnt = 0;
    if numDims == 1
      nx = arraySize(1);
      fprintf(fid, '(%d)\n', nx);
      for i = 1 : nx
        if maxCnt > 0 && cnt > maxCnt
          break;
        endif
        fprintf(fid, '%d, %12.10e\n', i-1, array(i));
        cnt = cnt+1;
      end
    elseif numDims == 2
      [nx, ny] = size(array);
      fprintf(fid, '(%d, %d)\n', nx, ny);
      x0 = 1; x1 = nx; y0 = 1; y1 = ny;

      for i = x0 : x1
        for j = y0 : y1
          if maxCnt > 0 && cnt > maxCnt
            break;
          endif
          fprintf(fid, '%d, %d, %12.10e\n', i-1,j-1, array(i,j));
          cnt = cnt+1;
        end
      end
    elseif numDims == 3
     [nx, ny, nz] = size(array);
     fprintf(fid, '(%d, %d, %d)\n', nx, ny, nz);
     x0 = 1; x1 = nx; y0 = 1; y1 = ny; z0 = 1; z1 = nz;
     x1 = x1 / 2;
     y1 = y1 / 2;
     z1 = z1 / 2;

     for i = x0 : x1
       for j = y0 : y1
          for k = z0 : z1
            if maxCnt > 0 && cnt > maxCnt
              break;
            endif
            fprintf(fid, '%d, %d, %d, %12.10e\n', i-1,j-1,k-1,array(i,j,k));
            cnt = cnt+1;
          end
        end
     end
    end

    fclose(fid);
end
