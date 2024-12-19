def write_array_to_file(array, filename):
    filename = filename + '.p.txt'
    #print('array.shape= ', array.shape)
    shp = f"({', '.join(map(str, array.shape))})"
    with open(filename, "w") as f:
        # Flatten the array and write each element to the file
        f.write('%s\n' % shp)
        cnt = 0
        for element in array.flatten():
            f.write('%12.10e\n' % element)
            cnt = cnt+1
            if cnt > 10000 : break

            
