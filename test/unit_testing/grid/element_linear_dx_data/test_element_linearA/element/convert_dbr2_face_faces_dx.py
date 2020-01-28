

inputf = 'dbr2_face_face1_dx.txt'
outputf= 'dbr2_face_face1_dx2fort.txt'

line_count = 1

fout = open(outputf,"w")

with open(inputf) as f:
    for iline in f:
        if line_count == 1:
            line_str = iline.split()
            fout.write("test_dbr2_face_face1_dx(:,:,%s,%s) = reshape((/                                                                                 &\n" % (line_str[3].replace(',',''),line_str[5].replace(',','')))     
            line_count += 1


        elif line_count == 5:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            fout.write("                                /),(/4,4/))\n")
            fout.write("\n")
            line_count = 1

        else:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk, &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            line_count += 1
              
fout.close()




inputf = 'dbr2_face_face2_dx.txt'
outputf= 'dbr2_face_face2_dx2fort.txt'

line_count = 1

fout = open(outputf,"w")

with open(inputf) as f:
    for iline in f:
        if line_count == 1:
            line_str = iline.split()
            fout.write("test_dbr2_face_face2_dx(:,:,%s,%s) = reshape((/                                                                                 &\n" % (line_str[3].replace(',',''),line_str[5].replace(',','')))     
            line_count += 1


        elif line_count == 5:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            fout.write("                                /),(/4,4/))\n")
            fout.write("\n")
            line_count = 1

        else:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk, &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            line_count += 1
              
fout.close()




inputf = 'dbr2_face_face3_dx.txt'
outputf= 'dbr2_face_face3_dx2fort.txt'

line_count = 1

fout = open(outputf,"w")

with open(inputf) as f:
    for iline in f:
        if line_count == 1:
            line_str = iline.split()
            fout.write("test_dbr2_face_face3_dx(:,:,%s,%s) = reshape((/                                                                                 &\n" % (line_str[3].replace(',',''),line_str[5].replace(',','')))     
            line_count += 1


        elif line_count == 5:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            fout.write("                                /),(/4,4/))\n")
            fout.write("\n")
            line_count = 1

        else:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk, &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            line_count += 1
              
fout.close()




inputf = 'dbr2_face_face4_dx.txt'
outputf= 'dbr2_face_face4_dx2fort.txt'

line_count = 1

fout = open(outputf,"w")

with open(inputf) as f:
    for iline in f:
        if line_count == 1:
            line_str = iline.split()
            fout.write("test_dbr2_face_face4_dx(:,:,%s,%s) = reshape((/                                                                                 &\n" % (line_str[3].replace(',',''),line_str[5].replace(',','')))     
            line_count += 1


        elif line_count == 5:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            fout.write("                                /),(/4,4/))\n")
            fout.write("\n")
            line_count = 1

        else:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk, &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            line_count += 1
              
fout.close()




inputf = 'dbr2_face_face5_dx.txt'
outputf= 'dbr2_face_face5_dx2fort.txt'

line_count = 1

fout = open(outputf,"w")

with open(inputf) as f:
    for iline in f:
        if line_count == 1:
            line_str = iline.split()
            fout.write("test_dbr2_face_face5_dx(:,:,%s,%s) = reshape((/                                                                                 &\n" % (line_str[3].replace(',',''),line_str[5].replace(',','')))     
            line_count += 1


        elif line_count == 5:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            fout.write("                                /),(/4,4/))\n")
            fout.write("\n")
            line_count = 1

        else:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk, &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            line_count += 1
              
fout.close()




inputf = 'dbr2_face_face6_dx.txt'
outputf= 'dbr2_face_face6_dx2fort.txt'

line_count = 1

fout = open(outputf,"w")

with open(inputf) as f:
    for iline in f:
        if line_count == 1:
            line_str = iline.split()
            fout.write("test_dbr2_face_face6_dx(:,:,%s,%s) = reshape((/                                                                                 &\n" % (line_str[3].replace(',',''),line_str[5].replace(',','')))     
            line_count += 1


        elif line_count == 5:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk  &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            fout.write("                                /),(/4,4/))\n")
            fout.write("\n")
            line_count = 1

        else:
            line_str = iline.split()
            fout.write("                                %s_rk, %s_rk, %s_rk,%s_rk, &\n" % (line_str[0],line_str[1],line_str[2],line_str[3]))     
            line_count += 1
              
fout.close()
