

inputf = 'dinterp_xcoords_dx.txt'
outputf= 'dinterp_xcoords_dx_dx2fort.txt'

line_count = 1

fout = open(outputf,"w")

with open(inputf) as f:
    for iline in f:
        if line_count == 1:
            line_str = iline.split()
            fout.write("test_dcoords_dx = reshape((/                                                                                 &\n" )     
            line_count += 1


        elif line_count == 9:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3],line_str[4],line_str[5],line_str[6],line_str[7],line_str[8],line_str[9],line_str[10],line_str[11],line_str[12],line_str[13],line_str[14],line_str[15],line_str[16],line_str[17],line_str[18],line_str[19],line_str[20],line_str[21],line_str[22],line_str[23]))     
            fout.write("                                /),(/24,8/))\n")
            fout.write("\n")
            line_count = 1

        else:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk,  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3],line_str[4],line_str[5],line_str[6],line_str[7],line_str[8],line_str[9],line_str[10],line_str[11],line_str[12],line_str[13],line_str[14],line_str[15],line_str[16],line_str[17],line_str[18],line_str[19],line_str[20],line_str[21],line_str[22],line_str[23]))     
            line_count += 1
              




fout.close()




inputf = 'dinterp_ycoords_dx.txt'
outputf= 'dinterp_ycoords_dx_dx2fort.txt'

line_count = 1

fout = open(outputf,"w")

with open(inputf) as f:
    for iline in f:
        if line_count == 1:
            line_str = iline.split()
            fout.write("test_dcoords_dx = reshape((/                                                                                 &\n" )     
            line_count += 1


        elif line_count == 9:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3],line_str[4],line_str[5],line_str[6],line_str[7],line_str[8],line_str[9],line_str[10],line_str[11],line_str[12],line_str[13],line_str[14],line_str[15],line_str[16],line_str[17],line_str[18],line_str[19],line_str[20],line_str[21],line_str[22],line_str[23]))     
            fout.write("                                /),(/24,8/))\n")
            fout.write("\n")
            line_count = 1

        else:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk,  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3],line_str[4],line_str[5],line_str[6],line_str[7],line_str[8],line_str[9],line_str[10],line_str[11],line_str[12],line_str[13],line_str[14],line_str[15],line_str[16],line_str[17],line_str[18],line_str[19],line_str[20],line_str[21],line_str[22],line_str[23]))     
            line_count += 1
              




fout.close()





inputf = 'dinterp_zcoords_dx.txt'
outputf= 'dinterp_zcoords_dx_dx2fort.txt'

line_count = 1

fout = open(outputf,"w")

with open(inputf) as f:
    for iline in f:
        if line_count == 1:
            line_str = iline.split()
            fout.write("test_dcoords_dx = reshape((/                                                                                 &\n" )     
            line_count += 1


        elif line_count == 9:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3],line_str[4],line_str[5],line_str[6],line_str[7],line_str[8],line_str[9],line_str[10],line_str[11],line_str[12],line_str[13],line_str[14],line_str[15],line_str[16],line_str[17],line_str[18],line_str[19],line_str[20],line_str[21],line_str[22],line_str[23]))     
            fout.write("                                /),(/24,8/))\n")
            fout.write("\n")
            line_count = 1

        else:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk, %s_rk, %s_rk,%s_rk,%s_rk,%s_rk,%s_rk,%s_rk,  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3],line_str[4],line_str[5],line_str[6],line_str[7],line_str[8],line_str[9],line_str[10],line_str[11],line_str[12],line_str[13],line_str[14],line_str[15],line_str[16],line_str[17],line_str[18],line_str[19],line_str[20],line_str[21],line_str[22],line_str[23]))     
            line_count += 1
              




fout.close()


