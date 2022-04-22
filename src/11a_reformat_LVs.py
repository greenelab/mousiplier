import argparse

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='the transformed values of latent variables')
    parser.add_argument('outfile', help='the reformated values of latent variables')
    args = parser.parse_args()
    
    outfh = open(args.outfile, 'w')
    ### the header
    outfh.write("LV_ID\tday\ttreatment\tregion\tlv_value\n")
    
    with open(args.infile) as f:
        header = f.readline()
        header = header.strip().split()
        LV_count = len(header) - 1
        #print(LV_count)
        #print(header)
        for line in f:
            line = line.strip().split()
            (day, treatment, region, rep) = line[0].split('_')

            for i in range(1, LV_count + 1):
                lv = "LV" + str(i)
                outfh.write(lv + '\t' + day + '\t' + treatment + '\t' + region + '\t' + str(line[i]) + '\n')

    outfh.close()

if __name__ == "__main__":
    main(sys.argv[1:])
