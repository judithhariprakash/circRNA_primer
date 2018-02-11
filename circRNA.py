import sys
import requests
import os
import uuid
import json

SPECIES = 'homo_sapiens'
ASSEMBLY = 'GRCh38'
OFFSET = 200

def get_input():
    try:
        coord = sys.argv[1].rstrip('\n')
        strand = sys.argv[2].rstrip('\n')
        outfile = sys.argv[3].rstrip('\n')
    except IndexError:
        print "USAGE:\npython circRNA.py coordinates strand outfile\n"
        print "\t'coordinates': chrom:start-end format (exp: 10:25000000-25001000)"
        print "\t'strand'     : can be '+' or '-'"
        print "\t'outfile'    : Give a name where the output will be saved in CSV format"
        sys.exit()
    coord = coord.rstrip('\n').rstrip('\r')
    strand = strand.rstrip('\n').rstrip('\r')
    try:
        chrom, start, end = parse_coords(coord)
    except:
        print "ERROR: Please provide coordinate in chrom:start-end format (exp: 10:25000000-25001000)"
        sys.exit()
    if strand not in ['+', '-']:
        print "ERROR: Strand should be either '+' or '-'"
        sys.exit()
    return chrom.upper().replace('CHR', ''), start, end, strand, outfile

def parse_coords(coord):
    chrom = coord.split(":")[0]
    start = coord.split(":")[1].split("-")[0]
    end = coord.split(":")[1].split("-")[1]
    print "INFO: Input coordinates successfully parsed as following:"
    print "INFO: Chrom:%s Start:%s End:%s" % (chrom, start, end)
    return chrom, start, end

def get_circ_coordinates(start, end, strand, offset=OFFSET):
    # if its positive strand the offset is subtracted from end and add to the start sequence
    if strand == "+":
        start1 = int(end) - offset
        end1 = int(end)
        start2 = int(start)
        end2 = int(start) + offset
    # if its negative strand, reverse complement of the obtained sequence is joined together
    elif strand == "-":
        start1 = int(start)
        end1 = int(start) + offset
        start2 = int(end) - offset
        end2 = int(end)
    else:
        return False
    print "INFO: Calculation of backspliced coordinates with offset of %dbp complete" % offset
    return ((start1, end1), (start2, end2))

def fetch_ensembl(chrom, start, end, strand, species=SPECIES, assembly=ASSEMBLY):
    print "INFO: Will connect to ENSEMBL for species: %s and assembly:%s" % (species, assembly)
    server = "https://rest.ensembl.org"
    s = '1' if strand == '+' else '-1'
    ext = "/sequence/region/%s/%s:%s..%s:%s?coord_system_version=%s" % (
        species, chrom, str(start), str(end), s, assembly)
    print "INFO:", server+ext
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
    if not r.ok:
        print "ERROR: Unable to contact Ensembl to fetch sequences. Exiting now.."
        sys.exit()
    print "INFO: Ensembl query successful"
    return r.text

def get_backspliced_seq(circ_cords):
    print "INFO: Connecting to Ensembl for fetching 5' sequence"
    seq1 = fetch_ensembl(chrom, circ_cords[0][0], circ_cords[0][1], strand)
    print "INFO: Connecting to Ensembl for fetching 3' sequence"
    seq2 = fetch_ensembl(chrom, circ_cords[1][0], circ_cords[1][1], strand)
    # Strand info is allow implicit in seq1 and seq2
    return seq1+seq2

def run_primer3(seq):
    print "INFO: Creating template file for Primer3 with template length of: %d" % len(seq)
    temp_fn = str(uuid.uuid4())
    with open(temp_fn, 'w') as OUT:
        OUT.write("SEQUENCE_ID=sequence\nSEQUENCE_TEMPLATE=%s\n=\n" % seq)
    script_dir = os.path.dirname(os.path.realpath(__file__))
    exe = "%s/primer3_bin/primer3_core" % script_dir
    settings_fn = "%s/primer3_bin/circrna_primers_settings.p3" % script_dir
    print "INFO: Using Primer3 settings file: %s" % settings_fn
    cmd = "%s -p3_settings_file=%s -echo_settings_file %s" % (exe, settings_fn, temp_fn)
    print "INFO: Running Primer3..."
    print "INFO: ", cmd
    output = os.popen(cmd).read()
    os.system('rm %s' % temp_fn)
    print "INFO: Removed temporary input file for Primer3"
    return output

def prep_output(raw_output):
    output_dict = {}
    for line in raw_output.split('\n'):
        cols = line.rstrip('\n').split('=')
        try:
            output_dict[cols[0]] = cols[1]
        except IndexError:
            pass
    num_primers =  int(output_dict['PRIMER_PAIR_NUM_RETURNED'])
    valid_primers = 0
    if num_primers > 0:
        parsed = [','.join([
            "Primer ID", "Product size", "Left primer", "Right primer",
            "Left GC", "Right GC", "Left TM", "Right TM",
            "Left pos", "Right pos", "Left size", "Right size"
        ])]
        for i in range(num_primers):
            left_pos = output_dict['PRIMER_LEFT_%d' % i].split(',')[0]
            right_pos = output_dict['PRIMER_RIGHT_%d' % i].split(',')[0]
            if int(left_pos) < 195 and int(right_pos) > 205:
                valid_primers += 1
                parsed.append(",".join([
                    str(i), output_dict['PRIMER_PAIR_%d_PRODUCT_SIZE' % i],
                    output_dict['PRIMER_LEFT_%d_SEQUENCE' % i], output_dict['PRIMER_RIGHT_%d_SEQUENCE' % i],
                    output_dict['PRIMER_LEFT_%d_GC_PERCENT' % i], output_dict['PRIMER_RIGHT_%d_GC_PERCENT' % i],
                    output_dict['PRIMER_LEFT_%d_TM' % i], output_dict['PRIMER_RIGHT_%d_TM' % i],
                    left_pos, right_pos,
                    output_dict['PRIMER_LEFT_%d' % i].split(',')[1], output_dict['PRIMER_RIGHT_%d' % i].split(',')[1]
                ]))
    else:
        print "WARNING: No primers found"
        return "No primers found"
    print "INFO: %d primers found" % valid_primers
    return '\n'.join(parsed)

if __name__ == '__main__':
    print "Welcome to CircPrimer"    
    chrom, start, end, strand, outfile = get_input()
    cc = get_circ_coordinates(start, end, strand)
    back_seq = get_backspliced_seq(cc)
    primer3_output = run_primer3(back_seq)
    with open("%s_raw_primer3.csv" % outfile, 'w') as OUT:
        OUT.write(primer3_output)
    primers_table = prep_output(primer3_output)
    with open("%s.csv" % outfile, 'w') as OUT:
        OUT.write(primers_table)
    print "SUCCESS: Analysis complete.. Primers are saved in %s.csv" % outfile
