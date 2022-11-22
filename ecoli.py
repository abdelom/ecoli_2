import sys
import re


def parse_sam(file_sam, recipient, donor, out):
    with open(f'sam/{file_sam}') as filin, open(f'position/{recipient}_{out}_pos.txt', "w")  as rec_file, \
    open(f'position/{donor}_{out}_pos.txt', "w") as don_file:
        d = {recipient: rec_file, donor: don_file}
        for line in filin:
            if line.startswith("@"):
                continue
            list_line = line[:-1].split()
            score = int(re.compile(r'AS:i:(\d+)').findall(line)[0])
            score_alternatif = int(re.compile(r'XS:i:(\d+)').findall(line)[0])
            if list_line[2] == "Plasmid" or list_line[2] == "*" or score == 0 or \
            score == score_alternatif:
                continue
            if re.compile(r'XA:Z:.+,.\d+,.*,\d+;.+,.\d+,.*,\d+;').findall(line):
                continue
            second_pos = -1
            second = re.compile(r'XA:Z:(.+),.(\d+),.*,\d+;').findall(line)
            if second:
                if second[0][0] == list_line[2]:
                    continue
                if score - score_alternatif < 30: 
                    second_pos = second[0][1]
            d[list_line[2]].write(f"{list_line[3]}\t{second_pos}\n")


def main():
    parse_sam(sys.argv[1] + ".sam", sys.argv[2], sys.argv[3], sys.argv[1])


if __name__=="__main__":
    main()

