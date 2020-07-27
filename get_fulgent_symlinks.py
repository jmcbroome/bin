#helper script that automates process of creating symlinks for each data file in a targeted Fulgent sequencing output directory and places the symlinks in the current directory.
import glob
import subprocess
import argparse
import hashlib

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('target_directory', help = 'Path to directory to read barcodes and create symlinks from.')
    parser.add_argument('-c', '--cut', type = int, help = 'Number of split sections (delineated by -) to remove from the end of sample names when processing. Default 1', default = 1)
    parser.add_argument('-a', '--hash', action = 'store_true', help = 'Calculate md5 hashes for each file and compare to the correct hashes in the md5sum text file, as well as creating symlinks.')
    args = parser.parse_args()
    return args

def get_hash(filename):
    md5_hash = hashlib.md5()
    with open(filename,"rb") as f:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: f.read(4096),b""):
            md5_hash.update(byte_block)
        return md5_hash.hexdigest()

def check_md5sum(md5path, dirpath):
    cont = True
    with open(md5path) as inf:
        for entry in inf:
            code, fpath = entry.strip().split('  ')
            print("Checking hash of entry {}".format(fpath))
            #add the general path to the specific path here.
            full_path = dirpath + fpath[1:]
            fresh_code = get_hash(full_path)
            if code != fresh_code:
                print('Warning: md5sums for this entry do not match!')
                print(code, fresh_code)
                if full_path.split('.')[-2:] != ['fastq', 'gz']:
                    print("Allowing linking, as this file is not a core data file.")
                else:
                    cont = False
    return cont

def read_barcode(path, cut = 1):
    bcd = {}
    with open(path) as inf:
        for entry in inf:
            spent = entry.strip().split(',')
            if spent[0] != 'AccessionID':
                bcd[spent[0][3:]] = '-'.join(spent[1].split('-')[:-cut])
    return bcd

def identify_bar(files):
    bar = None
    for f in files:
        sf = f.split('/')[-1]
        if sf[:7] == 'barcode':
            bar = f
    if bar is None:
        print("Can't identify barcode file, quitting")
    return bar

def link_files(files, bcd):
    for f in files:
        #the ones we care about here are actually directories
        #check for an ending, and if there isn't one, re-glob that path to get its contents.
        if '.' not in f.split('/')[-1]:
            #symlink to the files with the appropriate endings in SD.
            sd_r1 = glob.glob(f + '/*R1*fastq.gz')[0]
            sd_r2 = glob.glob(f + '/*R2*fastq.gz')[0]
            #choose a new name based on bcd.
            try:
                name_r1 = './' + bcd[sd_r1.split('/')[-1].split('_')[1]] + "_R1.fq.gz"
                name_r2 = './' + bcd[sd_r2.split('/')[-1].split('_')[1]] + "_R2.fq.gz"
            except KeyError:
                print("Error matching file names")
                print(sd_r1, sd_r2)
                continue
            subprocess.check_call(['ln', '-s', sd_r1, name_r1])
            subprocess.check_call(['ln', '-s', sd_r2, name_r2])
    print("Linking complete.")

def main():
    args = argparser()
    files = glob.glob(args.target_directory + '/*')
    #before you do anything else, verify that the md5 hashes match up with the listings in md5sum.
    if args.hash:
        print("Checking hashes.")
        proceed = check_md5sum([f for f in files if f.split('.')[-1] == 'md5sum'][0], args.target_directory)
    else:
        print("Skipping hash checking (default behavior)")
        proceed = True
    if proceed:
        #identify the barcode file among files.
        bcp = identify_bar(files)
        if bcp != None:
            bcd = read_barcode(bcp)
            #print(bcd.keys())
            #now go through the files and do symlinking.
            link_files(files, bcd)

if __name__ == "__main__":
    main()