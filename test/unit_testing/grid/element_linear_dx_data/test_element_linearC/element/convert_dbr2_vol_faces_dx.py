

inputf = 'dbr2_vol_face1_dx.txt'
outputf= 'dbr2_vol_face1_dx2fort.txt'

line_count = 1

fout = open(outputf,"w")

with open(inputf) as f:
    for iline in f:
        if line_count == 1:
            line_str = iline.split()
            fout.write("test_dbr2_vol_face1_dx(:,:,%s,%s) = reshape((/                                                                                 &\n" % (line_str[3].replace(',',''),line_str[5].replace(',','')))     
            line_count += 1


        elif line_count == 9:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            fout.write("                                /),(/4,8/))\n")
            fout.write("\n")
            line_count = 1

        else:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk, &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            line_count += 1
              
fout.close()





inputf = 'dbr2_vol_face2_dx.txt'
outputf= 'dbr2_vol_face2_dx2fort.txt'

line_count = 1

fout = open(outputf,"w")

with open(inputf) as f:
    for iline in f:
        if line_count == 1:
            line_str = iline.split()
            fout.write("test_dbr2_vol_face2_dx(:,:,%s,%s) = reshape((/                                                                                 &\n" % (line_str[3].replace(',',''),line_str[5].replace(',','')))     
            line_count += 1


        elif line_count == 9:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            fout.write("                                /),(/4,8/))\n")
            fout.write("\n")
            line_count = 1

        else:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk, &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            line_count += 1

fout.close()





inputf = 'dbr2_vol_face3_dx.txt'
outputf= 'dbr2_vol_face3_dx2fort.txt'

line_count = 1

fout = open(outputf,"w")

with open(inputf) as f:
    for iline in f:
        if line_count == 1:
            line_str = iline.split()
            fout.write("test_dbr2_vol_face3_dx(:,:,%s,%s) = reshape((/                                                                                 &\n" % (line_str[3].replace(',',''),line_str[5].replace(',','')))     
            line_count += 1


        elif line_count == 9:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            fout.write("                                /),(/4,8/))\n")
            fout.write("\n")
            line_count = 1

        else:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk, &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            line_count += 1

fout.close()





inputf = 'dbr2_vol_face4_dx.txt'
outputf= 'dbr2_vol_face4_dx2fort.txt'

line_count = 1

fout = open(outputf,"w")

with open(inputf) as f:
    for iline in f:
        if line_count == 1:
            line_str = iline.split()
            fout.write("test_dbr2_vol_face4_dx(:,:,%s,%s) = reshape((/                                                                                 &\n" % (line_str[3].replace(',',''),line_str[5].replace(',','')))     
            line_count += 1


        elif line_count == 9:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            fout.write("                                /),(/4,8/))\n")
            fout.write("\n")
            line_count = 1

        else:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk, &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            line_count += 1

fout.close()





inputf = 'dbr2_vol_face5_dx.txt'
outputf= 'dbr2_vol_face5_dx2fort.txt'

line_count = 1

fout = open(outputf,"w")

with open(inputf) as f:
    for iline in f:
        if line_count == 1:
            line_str = iline.split()
            fout.write("test_dbr2_vol_face5_dx(:,:,%s,%s) = reshape((/                                                                                 &\n" % (line_str[3].replace(',',''),line_str[5].replace(',','')))     
            line_count += 1


        elif line_count == 9:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            fout.write("                                /),(/4,8/))\n")
            fout.write("\n")
            line_count = 1

        else:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk, &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            line_count += 1

fout.close()





inputf = 'dbr2_vol_face6_dx.txt'
outputf= 'dbr2_vol_face6_dx2fort.txt'

line_count = 1

fout = open(outputf,"w")

with open(inputf) as f:
    for iline in f:
        if line_count == 1:
            line_str = iline.split()
            fout.write("test_dbr2_vol_face6_dx(:,:,%s,%s) = reshape((/                                                                                 &\n" % (line_str[3].replace(',',''),line_str[5].replace(',','')))     
            line_count += 1


        elif line_count == 9:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            fout.write("                                /),(/4,8/))\n")
            fout.write("\n")
            line_count = 1

        else:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk, &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            line_count += 1

fout.close()
