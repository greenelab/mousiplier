import sys

def main(argv):
    (infile, outfile) = argv
    outfh = open(outfile, 'w')
    outfh.write("LV_ID\tday\ttreatment\tregion\tlv_value\n")
    with open(infile) as f:
        header = f.readline()
        header = header.strip().split()
        LV_count = len(header) - 1
        #print(LV_count)
        #print(header)
        for line in f:
            line = line.strip().split()
            (day, treatment, region, rep) = line[0].split('_')
            #if treatment == "food":
            #    continue
            #print(region)
            for i in range(1, LV_count + 1):
                lv = "LV" + str(i)
                outfh.write(lv + '\t' + day + '\t' + treatment + '\t' + region + '\t' + str(line[i]) + '\n')



if __name__ == "__main__":
    main(sys.argv[1:])