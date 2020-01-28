

inputf = 'dmetric_dx.txt'
outputf= 'dmetric_dx2fort.txt'

line_count = 1

fout = open(outputf,"w")

with open(inputf) as f:
    for iline in f:
        if line_count == 1:
            line_str = iline.split()
            fout.write("test_dmetric_dx_zetamax(:,:,%s,%s,%s) = reshape((/                                                                                 &\n" % (line_str[3].replace(',',''),line_str[5].replace(',',''),line_str[7].replace(',','')))     
            line_count += 1


        elif line_count == 4:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk  &\n" % (line_str[0],line_str[1],line_str[2]))     
            fout.write("                                /),(/3,3/))\n")
            fout.write("\n")
            line_count = 1

        else:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk, &\n" % (line_str[0],line_str[1],line_str[2]))     
            line_count += 1
              




fout.close()


