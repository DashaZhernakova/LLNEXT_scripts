import sys

f = open(sys.argv[1], mode = 'rt')

f.readline()

clumped = {}
for l in f:
    spl = l.rstrip().split("\t")
    feature = spl[0] + ":" + spl[1]
    chr = spl[2]
    spl[4] = int(spl[4])
    pos = spl[4]
    spl[11] = float(spl[11])
    pval = spl[11]
    
    saved_lines = clumped.get(feature, [])
    toappend = 1
    if len(saved_lines) > 0:
        for i in range(len(saved_lines)):
            saved_line = saved_lines[i]
            saved_chr = saved_line[2]
            saved_pos = saved_line[4]
            saved_pval = saved_line[11]
            if chr == saved_chr and abs(pos - saved_pos) < 250000:
                if pval < saved_pval:
                    saved_lines[i] = spl
                    #print("replacing")
                toappend = 0
                break
        if toappend:
            saved_lines.append(spl)
            #print("appending")
    else:
        saved_lines = [spl]
        clumped[feature] = saved_lines
        #print("new item")

for lines in clumped.values():
    for res_line in lines:
        print ("\t".join(map(str,res_line)))

